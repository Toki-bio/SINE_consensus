#!/bin/bash
# ............................................................................
# sine_pairwise_consensus.sh - Pairwise-alignment SINE consensus builder
# ............................................................................
#
# Approach:
#   Step 0: If no reference provided, auto-generate one using the
#           bootstrapped iterative consensus from sine_consensus.sh logic
#           (subsample -> MAFFT MSA -> mini-consensus -> accumulate -> converge)
#   Step 1: Each copy is pairwise-aligned (MAFFT --op 5) to that reference.
#           The alignment anchors the copy in the reference coordinate system.
#   Step 2: Per-position nucleotide frequencies are tallied across ALL copies.
#           The reference is used ONLY as a coordinate frame - its own bases
#           do NOT enter the frequency table.  Only the copies vote.
#
# Usage:  bash sine_pairwise_consensus.sh [options] input.fasta
#
# Options:
#   -r FILE    Reference consensus FASTA (skip Step 0 auto-generation)
#   -n SIZE    Max copies for pairwise step (default: 10000)
#   -t PCT     Base-call frequency threshold, 0-100 (default: 50)
#              A base must reach this % of ALL copies (gaps in denominator)
#   -c PCT     Min non-gap coverage to attempt a call, 0-100 (default: 30)
#   -j JOBS    Parallel pairwise-alignment jobs (default: 4)
#   -k         Keep intermediate files (freq table, insertion map, log)
#   -h         Show this help
#
# Step 0 (auto-reference) parameters:
#   -N SIZE    Subsample size per bootstrap iteration (default: 100)
#   -M ITERS   Max bootstrap iterations (default: 50)
#   -S THRESH  Bootstrap convergence threshold (default: 0.01)
#   -A ITERS   Anchor phase: activate after this many iters without convergence (default: 3)
#              Set to 0 to disable anchoring entirely.
#   -P PCT     Anchor phase: max anchor fraction of subsample, 0-100 (default: 50)
#
# Output:
#   ${BASENAME}_pw_consensus.fasta          - ungapped consensus
# With -k:
#   ${BASENAME}_pw_freqtable.tsv            - per-position frequencies
#   ${BASENAME}_pw_insertions.tsv           - insertion hotspot map
#   ${BASENAME}_pw_consensus.log            - run log
#   ${BASENAME}_pw_reference.fasta          - auto-generated reference (if no -r)
#
# Requires: mafft, seqkit, awk, shuf, bc

set -euo pipefail

# -- Defaults ---------------------------------------------------------------
REF_FILE=""
MAX_COPIES=10000
CONS_THRESH=30
MIN_COVERAGE=10
JOBS=4
KEEP=0
# Step 0 defaults (bootstrap)
BOOT_SUBSAMPLE=100
BOOT_MAX_ITERS=50
BOOT_STABILITY=0.01
BOOT_MIN_ITERS=5
BOOT_MAX_GAP=50
BOOT_MIN_COV=50
ANCHOR_AFTER=3
ANCHOR_MAX_PCT=50

# -- Usage ------------------------------------------------------------------
usage() {
    grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,1\}//'
    exit 0
}

# -- Parse arguments --------------------------------------------------------
while getopts "r:n:t:c:j:N:M:S:A:P:kh" opt; do
    case $opt in
        r) REF_FILE=$OPTARG ;;
        n) MAX_COPIES=$OPTARG ;;
        t) CONS_THRESH=$OPTARG ;;
        c) MIN_COVERAGE=$OPTARG ;;
        j) JOBS=$OPTARG ;;
        N) BOOT_SUBSAMPLE=$OPTARG ;;
        M) BOOT_MAX_ITERS=$OPTARG ;;
        S) BOOT_STABILITY=$OPTARG ;;
        A) ANCHOR_AFTER=$OPTARG ;;
        P) ANCHOR_MAX_PCT=$OPTARG ;;
        k) KEEP=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

INPUT_FASTA=${1:?Error: input FASTA required}
[ -f "$INPUT_FASTA" ] || { echo "Error: '$INPUT_FASTA' not found" >&2; exit 1; }

BASENAME=$(basename "$INPUT_FASTA" | sed 's/\.\(fasta\|fa\|fas\|bnk\)$//')
WORKDIR=$(mktemp -d "${BASENAME}_pw_XXXXXX")
LOG="${BASENAME}_pw_consensus.log"

# -- Cleanup ----------------------------------------------------------------
cleanup() {
    rm -rf "$WORKDIR"
    [ "$KEEP" -eq 0 ] && rm -f "$LOG"
}
trap cleanup EXIT

progress() { printf '\r\033[K%s' "$*" >&2; }

# ............................................................................
#  Step 0 - Auto-generate reference via bootstrap consensus
# ............................................................................

# -- Random subsample -------------------------------------------------------
random_subsample() {
    local INFILE=$1 N=$2 OUTFILE=$3
    seqkit seq -w 0 "$INFILE" \
      | awk '/^>/ { if (seq != "") print hdr "\t" seq; hdr = $0; seq = ""; next }
             { seq = seq $0 }
             END { if (seq != "") print hdr "\t" seq }' \
      | shuf -n "$N" \
      | awk -F'\t' '{ print $1 "\n" $2 }' \
      > "$OUTFILE"
}

# -- Replicate a consensus sequence N times --------------------------------
replicate_fasta() {
    local INFILE=$1 N=$2 OUTFILE=$3
    local header seq i
    header=$(grep '^>' "$INFILE" | head -1)
    seq=$(grep -v '^>' "$INFILE" | tr -d '\n')
    > "$OUTFILE"
    for i in $(seq 1 "$N"); do
        printf "%s_anchor%d\n%s\n" "$header" "$i" "$seq" >> "$OUTFILE"
    done
}

# -- Bootstrap consensus from alignment (gaps NOT in denominator) -----------
boot_consensus() {
    local ALN=$1 CONS=$2
    awk -v max_gap_pct="$BOOT_MAX_GAP" -v min_cov="$BOOT_MIN_COV" '
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
            if (nongap < nseq * min_cov_frac) { printf "-"; continue }
            if (nongap < total[i] * (1 - max_gap_frac)) { printf "-"; continue }
            best = ""; mx = 0; tie = 0
            split("A C G T", bases, " ")
            for (b = 1; b <= 4; b++) {
                ct = count[i, bases[b]] + 0
                if (ct > mx)             { mx = ct; best = bases[b]; tie = 0 }
                else if (ct == mx && ct > 0) { tie = 1 }
            }
            printf (tie ? "-" : best)
        }
        print ""
    }
    ' "$ALN" > "$CONS"
}

# -- Hamming distance -------------------------------------------------------
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

# -- Generate or load reference ---------------------------------------------
SEQ_COUNT=$(grep -c '^>' "$INPUT_FASTA")

if [ -n "$REF_FILE" ]; then
    seqkit seq -w 0 "$REF_FILE" \
      | awk '/^>/ { if (seq != "") { gsub(/[-.]/, "", seq); print hdr "\n" toupper(seq) }
                    hdr = $0; seq = ""; next }
             { seq = seq $0 }
             END { if (seq != "") { gsub(/[-.]/, "", seq); print hdr "\n" toupper(seq) } }' \
      > "$WORKDIR/reference.fasta"
    echo "Step 0: using provided reference $REF_FILE" >> "$LOG"
    progress "Step 0: using provided reference"
    echo "" >&2
else
    # -- Auto-generate via bootstrap consensus ------------------------------
    echo "=== Step 0: bootstrap consensus ===" > "$LOG"
    {
        echo "Input:  $INPUT_FASTA  ($SEQ_COUNT sequences)"
        echo "Params: subsample=$BOOT_SUBSAMPLE  max_iters=$BOOT_MAX_ITERS  stability=$BOOT_STABILITY"
        if [ "$ANCHOR_AFTER" -gt 0 ]; then
            echo "Anchor: enabled (after $ANCHOR_AFTER iters, ramp to ${ANCHOR_MAX_PCT}% of subsample)"
        else
            echo "Anchor: disabled"
        fi
        echo "---"
    } >> "$LOG"

    [ "$SEQ_COUNT" -lt 3 ] && { echo "Error: need at least 3 sequences, got $SEQ_COUNT" >&2; exit 1; }
    BSUB=$BOOT_SUBSAMPLE
    [ "$BSUB" -gt "$SEQ_COUNT" ] && BSUB=$SEQ_COUNT

    MASTER_RAW="$WORKDIR/boot_master_raw.fasta"
    MASTER_ALN="$WORKDIR/boot_master.aln.fasta"
    PREV_CONS="$WORKDIR/boot_prev.fasta"
    CURR_CONS="$WORKDIR/boot_curr.fasta"
    BEST_CONS="$WORKDIR/boot_best.fasta"
    BEST_CHANGE="1.0"

    # Iteration 1
    progress "Step 0 iter 1/$BOOT_MAX_ITERS: subsampling..."
    random_subsample "$INPUT_FASTA" "$BSUB" "$WORKDIR/boot_sub_1.fasta"
    progress "Step 0 iter 1/$BOOT_MAX_ITERS: aligning..."
    mafft --auto --quiet "$WORKDIR/boot_sub_1.fasta" > "$WORKDIR/boot_aln_1.fasta"
    boot_consensus "$WORKDIR/boot_aln_1.fasta" "$WORKDIR/boot_mini_1.fasta"
    cp "$WORKDIR/boot_mini_1.fasta" "$MASTER_RAW"
    cp "$WORKDIR/boot_mini_1.fasta" "$MASTER_ALN"
    cp "$WORKDIR/boot_mini_1.fasta" "$PREV_CONS"
    cp "$WORKDIR/boot_mini_1.fasta" "$BEST_CONS"
    echo "Step0 Iter  1: init" >> "$LOG"
    rm -f "$WORKDIR/boot_sub_1.fasta" "$WORKDIR/boot_aln_1.fasta"

    # Main loop
    CONVERGED=0
    BITER=1
    ANCHOR_STEP=0

    for BITER in $(seq 2 "$BOOT_MAX_ITERS"); do

        # Anchor count for this iteration
        ANCHOR_N=0
        if [ "$ANCHOR_AFTER" -gt 0 ] && [ "$BITER" -gt "$ANCHOR_AFTER" ]; then
            ANCHOR_STEP=$((ANCHOR_STEP + 1))
            ANCHOR_PCT=$(echo "scale=0; $ANCHOR_MAX_PCT * $ANCHOR_STEP / $ANCHOR_AFTER" | bc)
            [ "$ANCHOR_PCT" -gt "$ANCHOR_MAX_PCT" ] && ANCHOR_PCT=$ANCHOR_MAX_PCT
            ANCHOR_N=$(echo "scale=0; $BSUB * $ANCHOR_PCT / 100" | bc)
        fi
        RANDOM_N=$((BSUB - ANCHOR_N))
        [ "$RANDOM_N" -lt 1 ] && RANDOM_N=1

        # Build subsample
        BSUB_FILE="$WORKDIR/boot_sub_${BITER}.fasta"
        random_subsample "$INPUT_FASTA" "$RANDOM_N" "$BSUB_FILE"

        if [ "$ANCHOR_N" -gt 0 ]; then
            replicate_fasta "$BEST_CONS" "$ANCHOR_N" "$WORKDIR/boot_anc_${BITER}.fasta"
            cat "$WORKDIR/boot_anc_${BITER}.fasta" >> "$BSUB_FILE"
            rm -f "$WORKDIR/boot_anc_${BITER}.fasta"
            progress "Step 0 iter $BITER/$BOOT_MAX_ITERS: aligning (${RANDOM_N} + ${ANCHOR_N} anchors)..."
        else
            progress "Step 0 iter $BITER/$BOOT_MAX_ITERS: aligning subsample..."
        fi

        mafft --auto --quiet "$BSUB_FILE" > "$WORKDIR/boot_aln_${BITER}.fasta"
        boot_consensus "$WORKDIR/boot_aln_${BITER}.fasta" "$WORKDIR/boot_mini_${BITER}.fasta"
        # Only add to master if pure random subsample (no anchors)
        if [ "$ANCHOR_N" -eq 0 ]; then
            cat "$WORKDIR/boot_mini_${BITER}.fasta" >> "$MASTER_RAW"
            cp "$WORKDIR/boot_mini_${BITER}.fasta" "$WORKDIR/boot_mini_clean_${BITER}.fasta"
        fi

        progress "Step 0 iter $BITER/$BOOT_MAX_ITERS: aligning master ($BITER mini-consensuses)..."
        mafft --auto --quiet "$MASTER_RAW" > "$MASTER_ALN"
        boot_consensus "$MASTER_ALN" "$CURR_CONS"

        CHANGE=$(approx_distance "$PREV_CONS" "$CURR_CONS")

        # Track best consensus
        if [ "$(echo "$CHANGE < $BEST_CHANGE" | bc -l)" -eq 1 ]; then
            BEST_CHANGE=$CHANGE
            cp "$CURR_CONS" "$BEST_CONS"
        fi

        if [ "$ANCHOR_N" -gt 0 ]; then
            echo "Step0 Iter $(printf '%2d' "$BITER"): change=$CHANGE [anchor ${ANCHOR_N}/${BSUB}]" >> "$LOG"
        else
            echo "Step0 Iter $(printf '%2d' "$BITER"): change=$CHANGE" >> "$LOG"
        fi
        progress "Step 0 iter $BITER/$BOOT_MAX_ITERS: change=$CHANGE"

        if [ "$BITER" -ge "$BOOT_MIN_ITERS" ] && [ "$(echo "$CHANGE < $BOOT_STABILITY" | bc -l)" -eq 1 ]; then
            echo "Step0 converged at iter $BITER (change=$CHANGE < $BOOT_STABILITY)" >> "$LOG"
            CONVERGED=1
            break
        fi

        cp "$CURR_CONS" "$PREV_CONS"
        rm -f "$BSUB_FILE" "$WORKDIR/boot_aln_${BITER}.fasta"
    done

    # Use best consensus if not converged (avoids outputting a mid-spike result)
    if [ "$CONVERGED" -eq 0 ]; then
        cp "$BEST_CONS" "$CURR_CONS"
        echo "Step0 Warning: did not converge - using best consensus (change=$BEST_CHANGE)" >> "$LOG"
    fi
    [ ! -f "$CURR_CONS" ] && cp "$PREV_CONS" "$CURR_CONS"


    # -- Step 0 refinement: re-align all mini-consensuses pairwise ----------
    # against the converged consensus to stabilise coordinate-frame drift.
    REFINE_REF="$WORKDIR/refine_ref.fasta"
    printf ">refine_ref\n" > "$REFINE_REF"
    grep -v '^>' "$CURR_CONS" | tr -d '-' >> "$REFINE_REF"
    REFINE_REF_LEN=$(grep -v '^>' "$REFINE_REF" | tr -d '\n' | wc -c)
    REFINE_PROJ="$WORKDIR/refine_proj.txt"
    > "$REFINE_PROJ"
    MINI_N=0
    for mf in "$WORKDIR"/boot_mini_clean_*.fasta; do
        [ -f "$mf" ] || continue
        MINI_N=$((MINI_N + 1))
        RPAIR="$WORKDIR/rp_${MINI_N}.fasta"
        RALN="$WORKDIR/ra_${MINI_N}.fasta"
        cat "$REFINE_REF" "$mf" > "$RPAIR"
        if mafft --auto --op 5 --quiet "$RPAIR" > "$RALN" 2>/dev/null; then
            awk '/^>/{n++;next}{s[n]=s[n]$0}
            END{ref=toupper(s[1]);copy=toupper(s[2])
                for(i=1;i<=length(ref);i++){r=substr(ref,i,1)
                    c=(i<=length(copy))?toupper(substr(copy,i,1)):"-"
                    if(r!="-"&&r!=".")printf "%s",c}
                print ""}' "$RALN" >> "$REFINE_PROJ"
        fi
        rm -f "$RPAIR" "$RALN"
    done
    if [ "$MINI_N" -gt 0 ] && [ -s "$REFINE_PROJ" ]; then
        REFINE_CONS="$WORKDIR/boot_refined.fasta"
        awk -v mgp="$BOOT_MAX_GAP" -v mc="$BOOT_MIN_COV" '
        {nseq++;line=toupper($0);if(length(line)>maxcol)maxcol=length(line)
         for(i=1;i<=length(line);i++){c=substr(line,i,1);count[i,c]++;total[i]++}}
        END{printf ">consensus\n"
            for(i=1;i<=maxcol;i++){if(!total[i])continue
                gap=(count[i,"-"]+0)+(count[i,"."]+0);nongap=total[i]-gap
                if(nongap<nseq*(mc/100.0)){printf"-";continue}
                if(nongap<total[i]*(1-mgp/100.0)){printf"-";continue}
                best="";mx=0;tie=0;split("A C G T",b," ")
                for(j=1;j<=4;j++){ct=count[i,b[j]]+0
                    if(ct>mx){mx=ct;best=b[j];tie=0}else if(ct==mx&&ct>0)tie=1}
                printf(tie?"-":best)}
            print""}' "$REFINE_PROJ" > "$REFINE_CONS"
        REFINE_LEN=$(grep -v '^>' "$REFINE_CONS" | tr -d '-' | tr -d '\n' | wc -c)
        if [ "$REFINE_LEN" -ge "$REFINE_REF_LEN" ]; then
            cp "$REFINE_CONS" "$CURR_CONS"
            progress "Step 0 refinement: ${REFINE_LEN} bp from ${MINI_N} minis (was ${REFINE_REF_LEN} bp)"
            echo "" >&2
            echo "Step0 refinement: ${MINI_N} minis -> ${REFINE_LEN} bp (was ${REFINE_REF_LEN} bp)" >> "$LOG"
        fi
    fi
    # Fallback: if best consensus is empty/too short, do simple plurality consensus
    BOOT_LEN_CHECK=$(grep -v '^>' "$CURR_CONS" | tr -d '-' | tr -d '\n' | wc -c)
    if [ "$BOOT_LEN_CHECK" -lt 10 ]; then
        echo "Step0 Warning: bootstrap produced short sequence ($BOOT_LEN_CHECK bp) - using fallback" >> "$LOG"
        echo "Warning: bootstrap result too short, using fallback consensus" >&2
        FB_N=100; [ "$FB_N" -gt "$SEQ_COUNT" ] && FB_N=$SEQ_COUNT
        random_subsample "$INPUT_FASTA" "$FB_N" "$WORKDIR/fb_sub.fasta"
        mafft --auto --quiet "$WORKDIR/fb_sub.fasta" > "$WORKDIR/fb_aln.fasta"
        awk -v thresh=35 -v min_cov=30 '
        /^>/ { if (seq!="") { nseq++; line=toupper(seq); if(length(line)>maxcol) maxcol=length(line);
                for(i=1;i<=length(line);i++){c=substr(line,i,1);count[i,c]++;total[i]++} } seq="";next }
        { seq=seq $0 }
        END { if(seq!="") { nseq++; line=toupper(seq); if(length(line)>maxcol) maxcol=length(line);
                for(i=1;i<=length(line);i++){c=substr(line,i,1);count[i,c]++;total[i]++} }
              printf ">consensus\n"
              for(i=1;i<=maxcol;i++) { if(!total[i]) continue
                gap=(count[i,"-"]+0)+(count[i,"."]+0); nongap=total[i]-gap
                if(nongap<nseq*(min_cov/100.0)){printf"-";continue}
                best="";mx=0;tie=0; split("A C G T",b," ")
                for(j=1;j<=4;j++){ct=count[i,b[j]]+0; if(ct>mx){mx=ct;best=b[j];tie=0} else if(ct==mx&&ct>0)tie=1}
                if(mx>0&&mx/nseq>=(thresh/100.0)) printf(tie?"N":best); else printf"-" }
              print"" }
        ' "$WORKDIR/fb_aln.fasta" > "$CURR_CONS"
    fi

    # Write reference (ungapped)
    printf ">%s_bootstrap_ref\n" "$BASENAME" > "$WORKDIR/reference.fasta"
    grep -v '^>' "$CURR_CONS" | tr -d '-' >> "$WORKDIR/reference.fasta"

    BOOT_LEN=$(grep -v '^>' "$WORKDIR/reference.fasta" | tr -d '\n' | wc -c)
    if [ "$CONVERGED" -eq 1 ]; then
        progress "Step 0 done: bootstrap reference ${BOOT_LEN} bp (converged at iter $BITER)"
    else
        progress "Step 0 done: bootstrap reference ${BOOT_LEN} bp (best at iter, not converged)"
    fi
    echo "" >&2
    echo "Step0 output: bootstrap reference ${BOOT_LEN} bp" >> "$LOG"

    if [ "$KEEP" -eq 1 ]; then
        cp -f "$WORKDIR/reference.fasta" "${BASENAME}_pw_reference.fasta"
    fi

    rm -f "$WORKDIR"/boot_*.fasta "$WORKDIR"/boot_*.tmp 2>/dev/null || true
fi

# ............................................................................
#  Step 1 - Prepare copies and reference
# ............................................................................
progress "Step 1: preparing sequences..."

REF_HEADER=$(head -1 "$WORKDIR/reference.fasta")
REF_SEQ=$(tail -1 "$WORKDIR/reference.fasta" | tr -d '\n' | tr '[:lower:]' '[:upper:]')
REF_LEN=${#REF_SEQ}

[ "$REF_LEN" -lt 10 ] && { echo "Error: reference too short ($REF_LEN bp)" >&2; exit 1; }

seqkit seq -w 0 "$INPUT_FASTA" \
  | awk '/^>/ { if (seq != "") { gsub(/[-.]/, "", seq); if (length(seq) > 0) print hdr "\t" seq }
                hdr = $0; seq = ""; next }
         { seq = seq $0 }
         END { if (seq != "") { gsub(/[-.]/, "", seq); if (length(seq) > 0) print hdr "\t" seq } }' \
  > "$WORKDIR/copies_linear.tsv"

COPY_COUNT=$(wc -l < "$WORKDIR/copies_linear.tsv")
shuf -n "$MAX_COPIES" "$WORKDIR/copies_linear.tsv" > "$WORKDIR/copies_sampled.tsv"
SEQ_COUNT=$(wc -l < "$WORKDIR/copies_sampled.tsv")

[ "$SEQ_COUNT" -lt 1 ] && { echo "Error: no copy sequences found" >&2; exit 1; }

awk -F'\t' -v dir="$WORKDIR" '{
    file = dir "/copy_" NR ".fasta"
    print $1 "\n" toupper($2) > file
    close(file)
}' "$WORKDIR/copies_sampled.tsv"

{
    echo "=== Pairwise alignment phase ==="
    echo "Reference: $REF_HEADER ($REF_LEN bp)"
    echo "Copies:    $SEQ_COUNT (of $COPY_COUNT total)"
    echo "Params:    threshold=$CONS_THRESH%  min_cov=$MIN_COVERAGE%  jobs=$JOBS"
    echo "---"
} >> "$LOG"

# ............................................................................
#  Step 2 - Pairwise MAFFT alignment + projection onto reference coords
# ............................................................................
DONE=0
while [ "$DONE" -lt "$SEQ_COUNT" ]; do
    BATCH_END=$((DONE + JOBS))
    [ "$BATCH_END" -gt "$SEQ_COUNT" ] && BATCH_END=$SEQ_COUNT

    for i in $(seq $((DONE + 1)) "$BATCH_END"); do
        (
            pair="$WORKDIR/pair_${i}.fasta"
            aln="$WORKDIR/aln_${i}.fasta"
            cat "$WORKDIR/reference.fasta" "$WORKDIR/copy_${i}.fasta" > "$pair"
            if ! mafft --auto --op 5 --quiet "$pair" > "$aln" 2>/dev/null; then
                touch "$WORKDIR/fail_${i}"
                rm -f "$pair" "$aln"
                exit 0
            fi
            awk -v ins_file="$WORKDIR/ins_${i}.txt" '
            /^>/ { n++; next }
            { s[n] = s[n] $0 }
            END {
                ref  = toupper(s[1])
                copy = toupper(s[2])
                rpos = 0; ilen = 0; ins = ""
                for (i = 1; i <= length(ref); i++) {
                    r = substr(ref, i, 1)
                    c = (i <= length(copy)) ? toupper(substr(copy, i, 1)) : "-"
                    if (r != "-" && r != ".") {
                        if (rpos > 0 && ilen > 0) ins = ins rpos ":" ilen " "
                        rpos++; ilen = 0
                        printf "%s", c
                    } else {
                        if (c != "-" && c != ".") ilen++
                    }
                }
                if (rpos > 0 && ilen > 0) ins = ins rpos ":" ilen
                print ""
                print ins > ins_file
            }
            ' "$aln" > "$WORKDIR/proj_${i}.txt"
            rm -f "$pair" "$aln"
        ) &
    done
    wait
    DONE=$BATCH_END
    progress "Step 2: pairwise alignment $DONE/$SEQ_COUNT"
done
progress ""

PROJECTIONS="$WORKDIR/projections.txt"
FAILED=0; GOOD=0
for i in $(seq 1 "$SEQ_COUNT"); do
    if [ -f "$WORKDIR/fail_${i}" ]; then
        FAILED=$((FAILED + 1))
    elif [ -f "$WORKDIR/proj_${i}.txt" ]; then
        cat "$WORKDIR/proj_${i}.txt" >> "$PROJECTIONS"
        GOOD=$((GOOD + 1))
    fi
done

echo "Alignments: $GOOD successful, $FAILED failed" >> "$LOG"
[ "$GOOD" -eq 0 ] && { echo "Error: all alignments failed" >&2; exit 1; }

PROJ_LEN=$(awk 'NR==1 { print length($0); exit }' "$PROJECTIONS")
BAD_LEN=$(awk -v elen="$PROJ_LEN" 'length($0) != elen { n++ } END { print n+0 }' "$PROJECTIONS")
[ "$BAD_LEN" -gt 0 ] && echo "Warning: $BAD_LEN projections with unexpected length" >> "$LOG"

progress "Step 2 done: $GOOD copies projected onto $REF_LEN reference positions"
echo "" >&2

# ............................................................................
#  Step 3 - Count frequencies and call consensus
# ............................................................................
progress "Step 3: counting frequencies and calling consensus..."

FREQ_TABLE="$WORKDIR/freqtable.tsv"
GAPPED="$WORKDIR/consensus_gapped.txt"
UNGAPPED="$WORKDIR/consensus_ungapped.txt"

awk -v thresh="$CONS_THRESH" -v min_cov="$MIN_COVERAGE" \
    -v ref="$REF_SEQ" \
    -v freq_file="$FREQ_TABLE" \
    -v gapped_file="$GAPPED" \
    -v ungapped_file="$UNGAPPED" '
{
    n++
    line = toupper($0)
    len = length(line)
    if (len > maxcol) maxcol = len
    for (i = 1; i <= len; i++) {
        c = substr(line, i, 1)
        count[i, c]++
    }
}
END {
    thr = thresh / 100.0
    mcov = min_cov / 100.0
    printf "Pos\tRef\tA\tC\tG\tT\tGap\tNongap\tTotal\tCalled\tFreq\n" > freq_file
    gapped = ""; ungapped = ""; called_count = 0; changed_count = 0
    for (i = 1; i <= maxcol; i++) {
        cA=count[i,"A"]+0; cC=count[i,"C"]+0; cG=count[i,"G"]+0; cT=count[i,"T"]+0
        gap=count[i,"-"]+0; nongap=cA+cC+cG+cT
        ref_base=(i<=length(ref))?substr(ref,i,1):"?"
        called="-"; freq_val=0
        if (nongap >= n * mcov) {
            best="-"; mx=0; tie=0; split("A C G T",bases," ")
            for (b=1;b<=4;b++) { ct=count[i,bases[b]]+0
                if(ct>mx){mx=ct;best=bases[b];tie=0} else if(ct==mx&&ct>0)tie=1 }
            freq_val=(n>0)?mx/n:0
            if (mx>0 && freq_val>=thr) called=(tie?"N":best)
        }
        gapped=gapped called
        if (called!="-"&&called!="N") {
            ungapped=ungapped called; called_count++
            if (called!=ref_base) changed_count++
        }
        printf "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%.4f\n", \
            i,ref_base,cA,cC,cG,cT,gap,nongap,n,called,freq_val > freq_file
    }
    print gapped > gapped_file
    print ungapped > ungapped_file
    printf "%d %d %d %d\n", maxcol, called_count, changed_count, n
}
' "$PROJECTIONS" > "$WORKDIR/stats.txt"

read TOTAL_POS CALLED CHANGED NCOPY < "$WORKDIR/stats.txt"
TRIMMED=$((TOTAL_POS - CALLED))

# ............................................................................
#  Step 4 - Insertion hotspot analysis
# ............................................................................
progress "Step 4: analysing insertion hotspots..."

INS_SUMMARY="$WORKDIR/insertion_summary.tsv"
cat "$WORKDIR"/ins_*.txt 2>/dev/null | awk '
{
    for (i=1;i<=NF;i++) { split($i,a,":"); pos=a[1];len=a[2]
        if(pos>0&&len>0){ins_n[pos]++;if(len>ins_max[pos])ins_max[pos]=len;ins_total[pos]+=len} }
}
END {
    print "After_pos\tCopies_with_ins\tMax_ins_len\tMean_ins_len"
    for(pos in ins_n) printf "%d\t%d\t%d\t%.1f\n",pos,ins_n[pos],ins_max[pos],ins_total[pos]/ins_n[pos]
}
' | sort -t$'\t' -k1,1n > "$INS_SUMMARY"

NOTABLE_INS=$(awk -F'\t' -v n="$NCOPY" 'NR>1 && $2/n>=0.10' "$INS_SUMMARY" | wc -l)

# ............................................................................
#  Step 5 - Output
# ............................................................................
OUTFILE="${BASENAME}_pw_consensus.fasta"
CONS_SEQ=$(cat "$UNGAPPED")
CONS_LEN=${#CONS_SEQ}

printf ">%s_pw\n%s\n" "$BASENAME" "$CONS_SEQ" > "$OUTFILE"

{
    echo "---"
    echo "Reference:  $REF_LEN bp"
    echo "Consensus:  $CONS_LEN bp (ungapped)"
    echo "Positions:  $CALLED/$TOTAL_POS called, $TRIMMED trimmed"
    echo "Changed:    $CHANGED positions differ from reference"
    echo "Insertions: $NOTABLE_INS sites with >=10% copy frequency"
    echo "Output:     $OUTFILE"
} >> "$LOG"

if [ "$KEEP" -eq 1 ]; then
    cp -f "$FREQ_TABLE"  "${BASENAME}_pw_freqtable.tsv"
    cp -f "$INS_SUMMARY" "${BASENAME}_pw_insertions.tsv"
    GAPPED_FILE="${BASENAME}_pw_consensus_gapped.fasta"
    printf ">%s_pw_gapped\n%s\n" "$BASENAME" "$(cat "$GAPPED")" > "$GAPPED_FILE"
    echo "Kept:       ${BASENAME}_pw_freqtable.tsv, ${BASENAME}_pw_insertions.tsv, ${BASENAME}_pw_consensus_gapped.fasta" >> "$LOG"
    trap 'rm -rf "$WORKDIR"' EXIT
fi

# -- Summary ----------------------------------------------------------------
progress ""
echo "" >&2
echo ".. $BASENAME pairwise consensus .." >&2
echo "  Reference: $REF_LEN bp  .  Copies: $NCOPY" >&2
echo "  Consensus: $CONS_LEN bp  .  Changed: $CHANGED positions" >&2
echo "  Called:    $CALLED/$TOTAL_POS positions  .  Trimmed: $TRIMMED" >&2
if [ "$NOTABLE_INS" -gt 0 ]; then
    echo "  ! $NOTABLE_INS insertion hotspot(s) detected (>=10% of copies)" >&2
    echo "    Check ${BASENAME}_pw_insertions.tsv (-k) or re-run with longer reference" >&2
fi
echo "  Output:    $OUTFILE" >&2
