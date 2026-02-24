#!/bin/bash
# sine_consensus.sh — Bootstrapped SINE consensus builder
# Usage: bash sine_consensus.sh input.fasta [subsample_size=100] [max_iters=50] [stability_thresh=0.01]
#
# Strategy: iterative subsampling + progressive consensus accumulation.
# Each iteration draws a random subsample, aligns with mafft --auto, then
# computes a majority-vote consensus via awk. Mini-consensuses are appended
# to master_raw.fasta and re-aligned with mafft --auto each iteration.
# Convergence checked via Hamming distance between successive master consensuses.
#
# NOTE: Master alignment contains mini-consensuses, not raw copies. Deliberate
# tradeoff — aligning all N copies simultaneously is intractable.
#
# Requires: mafft, shuf, awk (no seqkit/emboss dependency)
#
# Output:
#   final_consensus.fasta         — ungapped, ready for pairwise alignment use
#   final_consensus_gapped.fasta  — gap columns retained, shows alignment structure
#   master.aln.fasta              — accumulated master alignment (kept for inspection)
#   consensus_log.txt             — iteration log
#
# No automatic trimming. Trim manually after visual inspection.

INPUT_FASTA=${1:-input.fasta}
SUBSAMPLE_SIZE=${2:-100}
MAX_ITERS=${3:-50}
STABILITY_THRESH=${4:-0.01}

MASTER_RAW="master_raw.fasta"   # plain concatenation of all mini-consensuses
MASTER_ALN="master.aln.fasta"   # re-aligned each iteration with mafft --auto
PREV_CONS="prev_cons.fasta"     # consensus from previous iteration (for convergence check)
CURR_CONS="curr_cons.fasta"     # consensus from current iteration
LOG="consensus_log.txt"
BASENAME=$(basename "$INPUT_FASTA" .fasta)

# ── Consensus from alignment (majority-vote awk, single pass) ─────────────
# A base must appear at a column to be called; ties → N.
# Columns with >50% gaps → '-' (retained in gapped output, stripped in final).
# No external dependencies.
compute_consensus() {
    local ALN=$1
    local CONS=$2
    awk '
    /^>/ {
        if (seq != "") {
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i][c]++
                total[i]++
            }
        }
        seq = ""
        next
    }
    { seq = seq $0 }
    END {
        if (seq != "") {
            line = toupper(seq)
            if (length(line) > maxcol) maxcol = length(line)
            for (i = 1; i <= length(line); i++) {
                c = substr(line, i, 1)
                count[i][c]++
                total[i]++
            }
        }
        printf ">consensus\n"
        for (i = 1; i <= maxcol; i++) {
            if (!total[i]) continue
            gap    = (count[i]["-"] ? count[i]["-"] : 0)
            nongap = total[i] - gap
            if (nongap < total[i] / 2) { printf "-"; continue }
            best = ""; mx = 0; tie = 0
            for (b in count[i]) {
                if (b == "-") continue
                if (count[i][b] > mx)       { mx = count[i][b]; best = b; tie = 0 }
                else if (count[i][b] == mx) { tie = 1 }
            }
            printf (tie ? "N" : best)
        }
        print ""
    }
    ' "$ALN" > "$CONS"
}

# ── Hamming distance between two single-sequence FASTA files ──────────────
# Returns 1.0 if sequences differ in length by >10% (unconverged by definition).
approx_distance() {
    awk '
    NR==FNR { if (/^>/) next; seq1 = seq1 $0; next }
    /^>/ { next }
    { seq2 = seq2 $0 }
    END {
        len1 = length(seq1); len2 = length(seq2)
        if (len1 == 0 || len2 == 0) { print 1.0; exit }
        lendiff = (len1 > len2 ? len1 - len2 : len2 - len1)
        if (lendiff > len1 / 10) { print 1.0; exit }
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
echo "Starting bootstrapped consensus for $INPUT_FASTA" > "$LOG"
echo "Total sequences: $SEQ_COUNT  |  subsample=$SUBSAMPLE_SIZE  max_iters=$MAX_ITERS  thresh=$STABILITY_THRESH" >> "$LOG"

# ── Random subsample function (shuf/awk, no seqkit) ──────────────────────
# Linearize FASTA → shuf → head → re-format (from toki's sample script)
random_subsample() {
    local INFILE=$1 N=$2 OUTFILE=$3
    awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next } {printf("%s",$0)} END { printf("\n") }' "$INFILE" \
        | shuf | head -n "$N" \
        | awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' > "$OUTFILE"
}

# Iteration 1: bootstrap from first subsample
random_subsample "$INPUT_FASTA" "$SUBSAMPLE_SIZE" subsample_1.fasta
mafft --auto --quiet subsample_1.fasta > aligned_1.fasta
compute_consensus aligned_1.fasta "$PREV_CONS"
cp "$PREV_CONS" mini_cons_1.fasta
# Start raw collection with iter 1 mini-consensus
cp mini_cons_1.fasta "$MASTER_RAW"
# Align collection (1 seq — trivial, establishes file)
cp mini_cons_1.fasta "$MASTER_ALN"
echo "Iter 1: Initialized" >> "$LOG"

# ── Main loop ──────────────────────────────────────────────────────────────
CONVERGED=0
for ITER in $(seq 2 "$MAX_ITERS"); do
    random_subsample "$INPUT_FASTA" "$SUBSAMPLE_SIZE" subsample_${ITER}.fasta
    mafft --auto --quiet subsample_${ITER}.fasta > aligned_${ITER}.fasta
    compute_consensus aligned_${ITER}.fasta mini_cons_${ITER}.fasta

    # Append mini-consensus to raw collection, then re-align all with mafft --auto
    # This avoids mafft --add failing on single-sequence files
    cat mini_cons_${ITER}.fasta >> "$MASTER_RAW"
    mafft --auto --quiet "$MASTER_RAW" > "$MASTER_ALN"

    compute_consensus "$MASTER_ALN" "$CURR_CONS"
    CHANGE=$(approx_distance "$PREV_CONS" "$CURR_CONS")
    echo "Iter $ITER: change=$CHANGE" >> "$LOG"

    if [ "$ITER" -ge 5 ] && (( $(echo "$CHANGE < $STABILITY_THRESH" | bc -l) )); then
        echo "Converged at iter $ITER (change=$CHANGE < $STABILITY_THRESH)" >> "$LOG"
        CONVERGED=1
        break
    fi
    cp "$CURR_CONS" "$PREV_CONS"
done

[ "$CONVERGED" -eq 0 ] && echo "Warning: did not converge within $MAX_ITERS iterations" >> "$LOG"
[ ! -f "$CURR_CONS" ] && cp "$PREV_CONS" "$CURR_CONS"

# ── Output ─────────────────────────────────────────────────────────────────
cp "$CURR_CONS" final_consensus_gapped.fasta
sed -i "s/>consensus/>consensus_${BASENAME}_gapped/" final_consensus_gapped.fasta

printf ">consensus_%s\n" "$BASENAME" > final_consensus.fasta
grep -v '^>' "$CURR_CONS" | tr -d '-' >> final_consensus.fasta

echo "Output: final_consensus.fasta (ungapped), final_consensus_gapped.fasta" >> "$LOG"
echo "master.aln.fasta and master_raw.fasta kept for inspection." >> "$LOG"
echo "No trimming applied — trim manually after visual inspection." >> "$LOG"

# ── Convergence trace ──────────────────────────────────────────────────────
TRACE="convergence_trace.fasta"
> "$TRACE"
for f in $(ls mini_cons_*.fasta 2>/dev/null | sort -t_ -k3 -n); do
    ITER_NUM=$(echo "$f" | grep -oP '\d+')
    awk -v n="$ITER_NUM" '/^>/{print ">iter_"n; next}{print}' "$f" >> "$TRACE"
done
awk '/^>/{print ">FINAL_converged"; next}{print}' "$CURR_CONS" >> "$TRACE"
echo "Convergence trace: $TRACE (all mini-consensuses + final)" >> "$LOG"

# ── Cleanup (temp files only; keep master, trace and final output) ─────────
rm -f subsample_*.fasta aligned_*.fasta mini_cons_*.fasta \
      "$PREV_CONS" "$CURR_CONS" master_tmp.fasta 2>/dev/null

echo "Done. See $LOG for details."
