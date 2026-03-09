#!/bin/bash
# sine_consensus.sh ΓÇö Bootstrapped SINE consensus builder
#
# Usage:  bash sine_consensus.sh [options] input.fasta
#
# Options:
#   -n SIZE    Subsample size per iteration (default: 100)
#   -m ITERS   Maximum iterations (default: 50)
#   -s THRESH  Convergence threshold for Hamming distance (default: 0.01)
#   -g PCT     Max gap fraction per column to call a base, 0-100 (default: 50)
#   -c PCT     Min fraction of all sequences that must be non-gap to call a base, 0-100 (default: 50)
#              Increase this to trim ragged flanks more aggressively
#   -k         Keep intermediate files (trace alignment, log)
#   -h         Show this help
#
# Consensus algorithm (matches original/preferred script):
#   For each alignment column:
#     1. If gap fraction > max_gap_pct ΓåÆ emit '-'  (gap-dominant column)
#     2. Among non-gap characters {A,C,G,T}, find the plurality base
#     3. If a plurality base exists with no tie ΓåÆ call it
#     4. If tied ΓåÆ emit 'N'
#   Gaps are NOT in the frequency denominator (only non-gap chars compete).
#   This produces denser, longer consensuses than MSA-viewer mode.
#
# Output (always):
#   ${BASENAME}_consensus.fasta         ΓÇö ungapped, ready for pairwise use
# With -k:
#   ${BASENAME}_trace.aln.fasta         ΓÇö mini-consensus per iteration + final
#   ${BASENAME}_consensus.log           ΓÇö convergence log
#
# Requires: mafft, shuf, awk, bc

set -euo pipefail

# ΓöÇΓöÇ Defaults ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
SUBSAMPLE_SIZE=100
MAX_ITERS=50
STABILITY_THRESH=0.01
MAX_GAP_PCT=50      # columns with >this% gaps ΓåÆ '-'
MIN_COVERAGE=50     # columns where <this% of sequences are non-gap ΓåÆ '-' (trims ragged flanks)
KEEP=0
MIN_ITERS=5         # minimum iterations before convergence check

# ΓöÇΓöÇ Parse arguments ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
usage() {
    grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,1\}//'
    exit 0
}

while getopts "n:m:s:g:c:kh" opt; do
    case $opt in
        n) SUBSAMPLE_SIZE=$OPTARG ;;
        m) MAX_ITERS=$OPTARG ;;
        s) STABILITY_THRESH=$OPTARG ;;
        g) MAX_GAP_PCT=$OPTARG ;;
        c) MIN_COVERAGE=$OPTARG ;;
        k) KEEP=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

INPUT_FASTA=${1:?Error: input FASTA required}
[ -f "$INPUT_FASTA" ] || { echo "Error: $INPUT_FASTA not found" >&2; exit 1; }

BASENAME=$(basename "$INPUT_FASTA" | sed 's/\.\(fasta\|fa\|fas\|bnk\)$//')
WORKDIR=$(mktemp -d "${BASENAME}_cons_XXXXXX")
LOG="${BASENAME}_consensus.log"

# ΓöÇΓöÇ Cleanup on exit ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
cleanup() {
    rm -rf "$WORKDIR"
    [ "$KEEP" -eq 0 ] && rm -f "$LOG"
}
trap cleanup EXIT

# ΓöÇΓöÇ Progress display ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
progress() {
    printf '\r\033[K%s' "$*" >&2
}

# ΓöÇΓöÇ Random subsample (shuf/awk, no seqkit) ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
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

# ΓöÇΓöÇ Consensus from alignment ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
# Original majority-vote algorithm: gaps NOT in the frequency denominator.
#   - Columns with >MAX_GAP_PCT% gaps ΓåÆ '-'
#   - Among non-gap bases, plurality wins; ties ΓåÆ 'N'
# This is intentionally more aggressive than MSA-viewer (gaps-in-denominator)
# mode, producing denser consensuses that better represent the repeat family.
compute_consensus() {
    local ALN=$1 CONS=$2
    awk -v max_gap_pct="$MAX_GAP_PCT" -v min_cov="$MIN_COVERAGE" '
    /^>/ {
        if (seq != "") {
            nseq++
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i, c]++
                total[i]++
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
                count[i, c]++
                total[i]++
            }
        }
        max_gap_frac = max_gap_pct / 100.0
        min_cov_frac = min_cov / 100.0
        printf ">consensus\n"
        for (i = 1; i <= maxcol; i++) {
            if (!total[i]) continue
            gap    = (count[i, "-"] + 0) + (count[i, "."] + 0)
            nongap = total[i] - gap

            # Coverage check: require enough sequences to be non-gap here
            if (nongap < nseq * min_cov_frac) { printf "-"; continue }

            # Gap-dominant column among those present
            if (nongap < total[i] * (1 - max_gap_frac)) { printf "-"; continue }

            # Plurality among non-gap bases (gaps NOT in denominator)
            best = ""; mx = 0; tie = 0
            split("A C G T", bases, " ")
            for (b = 1; b <= 4; b++) {
                ct = count[i, bases[b]] + 0
                if (ct > mx)       { mx = ct; best = bases[b]; tie = 0 }
                else if (ct == mx && ct > 0) { tie = 1 }
            }
            printf (tie ? "-" : best)
        }
        print ""
    }
    ' "$ALN" > "$CONS"
}

# ΓöÇΓöÇ Hamming distance between two single-sequence FASTA files ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
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

# ΓöÇΓöÇ Init ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
SEQ_COUNT=$(grep -c '^>' "$INPUT_FASTA")
[ "$SEQ_COUNT" -lt 3 ] && { echo "Error: need at least 3 sequences, got $SEQ_COUNT" >&2; exit 1; }
[ "$SUBSAMPLE_SIZE" -gt "$SEQ_COUNT" ] && SUBSAMPLE_SIZE=$SEQ_COUNT

{
    echo "=== sine_consensus.sh ==="
    echo "Input:  $INPUT_FASTA  ($SEQ_COUNT sequences)"
    echo "Params: subsample=$SUBSAMPLE_SIZE  max_iters=$MAX_ITERS  max_gap=$MAX_GAP_PCT%  min_cov=$MIN_COVERAGE%  stability=$STABILITY_THRESH"
    echo "Algorithm: majority-vote, gaps excluded from base-frequency denominator"
    echo "---"
} > "$LOG"

MASTER_RAW="$WORKDIR/master_raw.fasta"
MASTER_ALN="$WORKDIR/master.aln.fasta"
PREV_CONS="$WORKDIR/prev_cons.fasta"
CURR_CONS="$WORKDIR/curr_cons.fasta"

# ΓöÇΓöÇ Iteration 1 ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
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

# ΓöÇΓöÇ Main loop ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
CONVERGED=0
ITER=1   # ensure defined even if MAX_ITERS=1
for ITER in $(seq 2 "$MAX_ITERS"); do
    progress "Iter $ITER/$MAX_ITERS: subsampling..."
    random_subsample "$INPUT_FASTA" "$SUBSAMPLE_SIZE" "$WORKDIR/sub_${ITER}.fasta"

    progress "Iter $ITER/$MAX_ITERS: aligning subsample..."
    mafft --auto --quiet "$WORKDIR/sub_${ITER}.fasta" > "$WORKDIR/aln_${ITER}.fasta"

    compute_consensus "$WORKDIR/aln_${ITER}.fasta" "$WORKDIR/mini_${ITER}.fasta"

    # Append mini-consensus; re-align all accumulated mini-consensuses
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

[ ! -f "$CURR_CONS" ] && cp "$PREV_CONS" "$CURR_CONS"

# ΓöÇΓöÇ Fallback: simple consensus when bootstrap does not converge ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
#  Take min(100, total) random copies ΓåÆ MAFFT MSA ΓåÆ call consensus using
#  plurality ΓëÑ35% with gaps in denominator (MSA-viewer semantics), min_cov=30%.
FALLBACK_SUBSAMPLE=100
FALLBACK_THRESH=35
FALLBACK_MIN_COV=30

if [ "$CONVERGED" -eq 0 ]; then
    echo "Warning: did not converge within $MAX_ITERS iterations ΓÇö using fallback consensus" >> "$LOG"
    echo "Warning: bootstrap did not converge after $MAX_ITERS iterations. Using fallback (simple) consensus." >&2

    FB_N=$FALLBACK_SUBSAMPLE
    [ "$FB_N" -gt "$SEQ_COUNT" ] && FB_N=$SEQ_COUNT

    progress "Fallback: subsampling $FB_N copies..."
    random_subsample "$INPUT_FASTA" "$FB_N" "$WORKDIR/fb_sub.fasta"

    progress "Fallback: aligning $FB_N copies with MAFFT..."
    mafft --auto --quiet "$WORKDIR/fb_sub.fasta" > "$WORKDIR/fb_aln.fasta"

    progress "Fallback: calling consensus (plurality >= ${FALLBACK_THRESH}%)..."
    # Consensus with gaps in denominator (MSA-viewer / EMBOSS cons style)
    awk -v thresh="$FALLBACK_THRESH" -v min_cov="$FALLBACK_MIN_COV" '
    /^>/ {
        if (seq != "") {
            nseq++
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i, c]++
                total[i]++
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
                count[i, c]++
                total[i]++
            }
        }
        thr = thresh / 100.0
        mcov = min_cov / 100.0
        printf ">consensus\n"
        for (i = 1; i <= maxcol; i++) {
            if (!total[i]) continue
            gap    = (count[i, "-"] + 0) + (count[i, "."] + 0)
            nongap = total[i] - gap

            # Coverage check: skip if too few non-gap
            if (nongap < nseq * mcov) { printf "-"; continue }

            # Find plurality base
            best = ""; mx = 0; tie = 0
            split("A C G T", bases, " ")
            for (b = 1; b <= 4; b++) {
                ct = count[i, bases[b]] + 0
                if (ct > mx)            { mx = ct; best = bases[b]; tie = 0 }
                else if (ct == mx && ct > 0) { tie = 1 }
            }

            # Threshold check: denominator = ALL sequences (gaps included)
            freq = (nseq > 0) ? mx / nseq : 0
            if (mx > 0 && freq >= thr) {
                printf (tie ? "N" : best)
            } else {
                printf "-"
            }
        }
        print ""
    }
    ' "$WORKDIR/fb_aln.fasta" > "$CURR_CONS"

    echo "Fallback: $FB_N copies, threshold=${FALLBACK_THRESH}%, min_cov=${FALLBACK_MIN_COV}%" >> "$LOG"
    progress "Fallback: done"
    echo "" >&2
fi

# ΓöÇΓöÇ Output ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
# Ungapped
OUTFILE="${BASENAME}_consensus.fasta"
printf ">%s\n" "$BASENAME" > "$OUTFILE"
grep -v '^>' "$CURR_CONS" | tr -d '-' >> "$OUTFILE"
CONS_LEN=$(grep -v '^>' "$OUTFILE" | tr -d '\n' | wc -c)

echo "Output: $OUTFILE (${CONS_LEN} bp)" >> "$LOG"

# ΓöÇΓöÇ Trace alignment (if -k) ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
if [ "$KEEP" -eq 1 ]; then
    TRACE="${BASENAME}_trace.aln.fasta"
    > "$TRACE"
    for f in $(ls "$WORKDIR"/mini_*.fasta 2>/dev/null | sort -t_ -k2 -n); do
        ITER_NUM=$(basename "$f" | grep -o '[0-9]*')
        awk -v n="$ITER_NUM" '/^>/{print ">iter_"n; next}{print}' "$f" >> "$TRACE"
    done
    awk -v name="$BASENAME" '/^>/{print ">FINAL_"name; next}{print}' "$CURR_CONS" >> "$TRACE"
    echo "Trace: $TRACE" >> "$LOG"
fi

if [ "$KEEP" -eq 1 ]; then
    # Prevent cleanup from deleting log when -k is set
    trap 'rm -rf "$WORKDIR"' EXIT
fi

# ΓöÇΓöÇ Summary ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
if [ "$CONVERGED" -eq 1 ]; then
    echo "$BASENAME: converged at iter $ITER ΓåÆ ${CONS_LEN} bp" >&2
else
    echo "$BASENAME: fallback consensus (plurality >= ${FALLBACK_THRESH}%, $FB_N copies) ΓåÆ ${CONS_LEN} bp" >&2
fi
