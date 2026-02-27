#!/bin/bash
# sine_consensus.sh — Bootstrapped SINE consensus builder
#
# Usage:  bash sine_consensus.sh [options] input.fasta
#
# Options:
#   -n SIZE    Subsample size per iteration (default: 100)
#   -m ITERS   Maximum iterations (default: 50)
#   -s THRESH  Convergence threshold for Hamming distance (default: 0.01)
#   -t PCT     Frequency threshold 0-100: base_count/total_seqs >= t to call (default: 50)
#              Gaps are in the denominator — this is what suppresses ragged flanks.
#   -c PCT     Min coverage 0-100: non_gap/total_seqs >= c to even try calling (default: 30)
#   -k         Keep intermediate files (trace alignment, log)
#   -h         Show this help
#
# Consensus algorithm (matches Toki-bio MSA-viewer exactly):
#   For each alignment column:
#     1. non_gap_count / nseq < min_cov  → '-'
#     2. best_base_count / nseq >= threshold  → call base  (gaps IN denominator)
#     3. Otherwise → '-'
#   Using gaps in the denominator is critical: flank columns present in only
#   20-30% of sequences score 0.2-0.3 frequency and are silently trimmed.
#   This also prevents flanks from accumulating into mini-consensuses during
#   iterative subsampling, which was the root cause of garbage output.
#
# Output:
#   ${BASENAME}_consensus.fasta  — ungapped consensus
# With -k:
#   ${BASENAME}_trace.aln.fasta  — mini-consensuses per iteration + final
#   ${BASENAME}_consensus.log    — convergence log
#
# Requires: mafft, shuf, awk, bc

set -euo pipefail

# ── Defaults ───────────────────────────────────────────────────────────────
SUBSAMPLE_SIZE=100
MAX_ITERS=50
STABILITY_THRESH=0.01
CONS_THRESH=50      # frequency threshold (gaps in denominator)
MIN_COVERAGE=30     # min non-gap coverage to attempt a call
KEEP=0
MIN_ITERS=5

# ── Parse arguments ────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,1\}//'
    exit 0
}

while getopts "n:m:s:t:c:kh" opt; do
    case $opt in
        n) SUBSAMPLE_SIZE=$OPTARG ;;
        m) MAX_ITERS=$OPTARG ;;
        s) STABILITY_THRESH=$OPTARG ;;
        t) CONS_THRESH=$OPTARG ;;
        c) MIN_COVERAGE=$OPTARG ;;
        k) KEEP=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

INPUT_FASTA=${1:?Error: input FASTA required}
[ -f "$INPUT_FASTA" ] || { echo "Error: $INPUT_FASTA not found" >&2; exit 1; }

BASENAME=$(basename "$INPUT_FASTA" | sed 's/\.\(fasta\|fa\|fas\)$//')
WORKDIR=$(mktemp -d "${BASENAME}_cons_XXXXXX")
LOG="${BASENAME}_consensus.log"

# ── Cleanup on exit ────────────────────────────────────────────────────────
cleanup() {
    rm -rf "$WORKDIR"
    [ "$KEEP" -eq 0 ] && rm -f "$LOG"
}
trap cleanup EXIT

# ── Progress display ───────────────────────────────────────────────────────
progress() {
    printf '\r\033[K%s' "$*" >&2
}

# ── Random subsample (shuf/awk, no seqkit) ────────────────────────────────
random_subsample() {
    local INFILE=$1 N=$2 OUTFILE=$3
    awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next }
         { printf("%s",$0) }
         END { printf("\n") }' "$INFILE" \
        | shuf > "$WORKDIR/shuffled.tmp"
    head -n "$N" "$WORKDIR/shuffled.tmp" \
        | awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' > "$OUTFILE"
    rm -f "$WORKDIR/shuffled.tmp"
}

# ── Consensus from alignment (MSA-viewer algorithm) ────────────────────────
# Gaps are IN the frequency denominator. A base must appear in >=threshold%
# of ALL sequences (not just those covering this column) to be called.
# This suppresses ragged flanks at every stage, including mini-consensuses,
# so garbage cannot accumulate into the master alignment.
compute_consensus() {
    local ALN=$1 CONS=$2
    awk -v thresh="$CONS_THRESH" -v min_cov="$MIN_COVERAGE" '
    /^>/ {
        if (seq != "") {
            nseq++
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i][c]++
            }
        }
        seq = ""; next
    }
    { seq = seq $0 }
    END {
        if (seq != "") {
            nseq++
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i][c]++
            }
        }
        thr = thresh / 100.0
        cov = min_cov / 100.0
        printf ">consensus\n"
        for (i = 1; i <= maxcol; i++) {
            gap    = (count[i]["-"] + 0) + (count[i]["."] + 0)
            nongap = nseq - gap

            # 1. Min coverage check
            if (nongap < nseq * cov) { printf "-"; continue }

            # 2. Most frequent base among ACGT
            split("A C G T", bases, " ")
            best = "-"; mx = 0; tie = 0
            for (b = 1; b <= 4; b++) {
                base = bases[b]
                c = count[i][base] + 0
                if (c > mx)             { mx = c; best = base; tie = 0 }
                else if (c == mx && c > 0) { tie = 1 }
            }

            # 3. Frequency with gaps in denominator — matches MSA-viewer
            if (mx == 0 || mx / nseq < thr) { printf "-"; continue }

            printf (tie ? "-" : best)
        }
        print ""
    }
    ' "$ALN" > "$CONS"
}

# ── Hamming distance between two single-sequence FASTA files ──────────────
approx_distance() {
    awk '
    NR==FNR { if (/^>/) next; seq1 = seq1 $0; next }
    /^>/ { next }
    { seq2 = seq2 $0 }
    END {
        len1 = length(seq1); len2 = length(seq2)
        if (len1 == 0 || len2 == 0) { print 1.0; exit }
        lendiff = (len1 > len2 ? len1 - len2 : len2 - len1)
        if (lendiff > len1 * 0.1) { print 1.0; exit }
        minlen = (len1 < len2 ? len1 : len2)
        diff = 0
        for (i = 1; i <= minlen; i++)
            if (substr(seq1, i, 1) != substr(seq2, i, 1)) diff++
        print diff / len1
    }
    ' "$1" "$2"
}

# ── Init ───────────────────────────────────────────────────────────────────
SEQ_COUNT=$(grep -c '^>' "$INPUT_FASTA")
[ "$SEQ_COUNT" -lt 3 ] && { echo "Error: need at least 3 sequences, got $SEQ_COUNT" >&2; exit 1; }
[ "$SUBSAMPLE_SIZE" -gt "$SEQ_COUNT" ] && SUBSAMPLE_SIZE=$SEQ_COUNT

{
    echo "=== sine_consensus.sh ==="
    echo "Input:  $INPUT_FASTA  ($SEQ_COUNT sequences)"
    echo "Params: subsample=$SUBSAMPLE_SIZE  max_iters=$MAX_ITERS  threshold=$CONS_THRESH%  min_cov=$MIN_COVERAGE%  stability=$STABILITY_THRESH"
    echo "Algorithm: MSA-viewer (gaps in frequency denominator)"
    echo "---"
} > "$LOG"

MASTER_RAW="$WORKDIR/master_raw.fasta"
MASTER_ALN="$WORKDIR/master.aln.fasta"
PREV_CONS="$WORKDIR/prev_cons.fasta"
CURR_CONS="$WORKDIR/curr_cons.fasta"

# ── Iteration 1 ────────────────────────────────────────────────────────────
progress "Iter 1/$MAX_ITERS: subsampling..."
random_subsample "$INPUT_FASTA" "$SUBSAMPLE_SIZE" "$WORKDIR/sub_1.fasta"

progress "Iter 1/$MAX_ITERS: aligning..."
mafft --auto --quiet "$WORKDIR/sub_1.fasta" > "$WORKDIR/aln_1.fasta"

compute_consensus "$WORKDIR/aln_1.fasta" "$WORKDIR/mini_1.fasta"

cp "$WORKDIR/mini_1.fasta" "$MASTER_RAW"
cp "$WORKDIR/mini_1.fasta" "$MASTER_ALN"
cp "$WORKDIR/mini_1.fasta" "$PREV_CONS"

echo "Iter  1: init" >> "$LOG"
rm -f "$WORKDIR/sub_1.fasta" "$WORKDIR/aln_1.fasta"

# ── Main loop ──────────────────────────────────────────────────────────────
CONVERGED=0
ITER=1
for ITER in $(seq 2 "$MAX_ITERS"); do
    progress "Iter $ITER/$MAX_ITERS: subsampling..."
    random_subsample "$INPUT_FASTA" "$SUBSAMPLE_SIZE" "$WORKDIR/sub_${ITER}.fasta"

    progress "Iter $ITER/$MAX_ITERS: aligning subsample..."
    mafft --auto --quiet "$WORKDIR/sub_${ITER}.fasta" > "$WORKDIR/aln_${ITER}.fasta"

    compute_consensus "$WORKDIR/aln_${ITER}.fasta" "$WORKDIR/mini_${ITER}.fasta"

    cat "$WORKDIR/mini_${ITER}.fasta" >> "$MASTER_RAW"

    progress "Iter $ITER/$MAX_ITERS: aligning master ($ITER mini-consensuses)..."
    mafft --auto --quiet "$MASTER_RAW" > "$MASTER_ALN"

    compute_consensus "$MASTER_ALN" "$CURR_CONS"

    CHANGE=$(approx_distance "$PREV_CONS" "$CURR_CONS")
    echo "Iter $(printf '%2d' "$ITER"): change=$CHANGE" >> "$LOG"
    progress "Iter $ITER/$MAX_ITERS: change=$CHANGE"

    if [ "$ITER" -ge "$MIN_ITERS" ] && [ "$(echo "$CHANGE < $STABILITY_THRESH" | bc -l)" -eq 1 ]; then
        echo "Converged at iter $ITER (change=$CHANGE < $STABILITY_THRESH)" >> "$LOG"
        CONVERGED=1
        break
    fi

    cp "$CURR_CONS" "$PREV_CONS"
    rm -f "$WORKDIR/sub_${ITER}.fasta" "$WORKDIR/aln_${ITER}.fasta"
done

progress ""
echo "" >&2

[ "$CONVERGED" -eq 0 ] && echo "Warning: did not converge within $MAX_ITERS iterations" >> "$LOG"
[ ! -f "$CURR_CONS" ] && cp "$PREV_CONS" "$CURR_CONS"

# ── Output ─────────────────────────────────────────────────────────────────
OUTFILE="${BASENAME}_consensus.fasta"
printf ">%s\n" "$BASENAME" > "$OUTFILE"
grep -v '^>' "$CURR_CONS" | tr -d '-' >> "$OUTFILE"
CONS_LEN=$(grep -v '^>' "$OUTFILE" | tr -d '\n' | wc -c)
echo "Output: $OUTFILE (${CONS_LEN} bp)" >> "$LOG"

# ── Trace alignment (if -k) ───────────────────────────────────────────────
if [ "$KEEP" -eq 1 ]; then
    TRACE="${BASENAME}_trace.aln.fasta"
    > "$TRACE"
    for f in $(ls "$WORKDIR"/mini_*.fasta 2>/dev/null | sort -t_ -k2 -n); do
        ITER_NUM=$(basename "$f" | grep -oP '\d+')
        awk -v n="$ITER_NUM" '/^>/{print ">iter_"n; next}{print}' "$f" >> "$TRACE"
    done
    awk -v name="$BASENAME" '/^>/{print ">FINAL_"name; next}{print}' "$CURR_CONS" >> "$TRACE"
    echo "Trace: $TRACE" >> "$LOG"
    trap 'rm -rf "$WORKDIR"' EXIT
fi

# ── Summary ───────────────────────────────────────────────────────────────
if [ "$CONVERGED" -eq 1 ]; then
    echo "$BASENAME: converged at iter $ITER → ${CONS_LEN} bp" >&2
else
    echo "$BASENAME: NOT converged after $MAX_ITERS iters → ${CONS_LEN} bp" >&2
fi
