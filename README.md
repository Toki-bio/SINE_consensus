# SINE Consensus Builder

Tools for building consensus sequences of SINE (Short Interspersed Nuclear Element) families from collections of genomic copies. Designed for repeat annotation pipelines where thousands of diverged copies must be distilled into a single representative consensus.

## Scripts

### `sine_consensus.sh` — Bootstrap consensus

Iterative bootstrap approach. Repeatedly subsamples copies, aligns them with MAFFT, extracts a mini-consensus, then re-aligns all accumulated mini-consensuses until the result stabilizes (Hamming distance converges below a threshold).

```
bash sine_consensus.sh [options] input.fasta
```

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-n SIZE` | 100 | Subsample size per iteration |
| `-m ITERS` | 50 | Maximum iterations |
| `-s THRESH` | 0.01 | Convergence threshold (Hamming distance) |
| `-g PCT` | 50 | Max gap fraction per column to call a base |
| `-c PCT` | 50 | Min non-gap coverage to call a base |
| `-k` | off | Keep intermediate files (trace, log) |

**Output:**
- `{name}_consensus.fasta` — ungapped consensus

**Algorithm:**
- Gaps are **excluded** from the frequency denominator (only non-gap characters A/C/G/T compete)
- Produces denser, longer consensuses than gap-in-denominator methods
- If bootstrap does not converge: falls back to a simple consensus from 100 random copies using plurality ≥ 35% with gaps in denominator

---

### `sine_pairwise_consensus.sh` — Pairwise-alignment consensus

Each copy is individually aligned to a reference using MAFFT with high gap-opening penalty (`--op 5`), then projected onto the reference coordinate system. Per-position nucleotide frequencies are tallied across all copies. The reference is used only as a coordinate frame — its own bases do not vote.

```
bash sine_pairwise_consensus.sh [options] input.fasta
```

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-r FILE` | — | Reference consensus FASTA (skip auto-generation) |
| `-n SIZE` | 10000 | Max copies for pairwise step |
| `-t PCT` | 50 | Base-call frequency threshold (gaps in denominator) |
| `-c PCT` | 30 | Min non-gap coverage to attempt a call |
| `-j JOBS` | 4 | Parallel alignment jobs |
| `-N SIZE` | 100 | Bootstrap subsample size (Step 0) |
| `-M ITERS` | 50 | Max bootstrap iterations (Step 0) |
| `-S THRESH` | 0.01 | Bootstrap convergence threshold (Step 0) |
| `-k` | off | Keep intermediate files (freq table, insertions, log) |

**Output:**
- `{name}_pw_consensus.fasta` — ungapped consensus
- With `-k`: frequency table, insertion hotspot map, log, reference

**Steps:**
1. **Step 0** — Auto-generate reference via bootstrap consensus (same algorithm as `sine_consensus.sh`), or use `-r` to provide one
2. **Step 1** — Pairwise MAFFT alignment of each copy to the reference (`--op 5`)
3. **Step 2** — Project each aligned copy onto reference coordinates
4. **Step 3** — Count per-position frequencies, call consensus (gaps in denominator)
5. **Step 4** — Insertion hotspot analysis
6. **Step 5** — Output

If bootstrap (Step 0) does not converge: falls back to a simple consensus from 100 random copies using plurality ≥ 35% with gaps in denominator, and skips the pairwise phase.

---

## When to use which

| Scenario | Script | Why |
|----------|--------|-----|
| **Quick consensus** from a SINE family | `sine_consensus.sh` | Fast, no reference needed, works with any number of copies |
| **High-copy family** (thousands of copies) | `sine_pairwise_consensus.sh` | Pairwise alignment scales better than full MSA; avoids MSA artifacts from massive datasets |
| **Reference already available** | `sine_pairwise_consensus.sh -r ref.fa` | Skips bootstrap, directly uses your reference as coordinate frame |
| **Diverged or fragmented copies** | `sine_pairwise_consensus.sh` | High gap-opening penalty (`--op 5`) prevents spurious gap insertions in diverged copies |
| **Subfamily splitting** | `sine_consensus.sh` per subfamily | Run separately on each cluster of copies after subfamily assignment |

## Consensus algorithms

### Bootstrap consensus (both scripts, Step 0 / main loop)

- Per alignment column, count bases among **non-gap** characters only
- Plurality base wins; ties produce `N`
- Column emits gap if gap fraction exceeds threshold or coverage is below minimum
- **Gaps excluded from denominator** — this produces longer, denser consensuses

### Pairwise consensus (sine_pairwise_consensus.sh, Steps 3–5)

- Per reference position, count what each copy has at that position after pairwise alignment
- **Gaps included in denominator** — a base must appear in ≥ threshold % of ALL copies (not just those with a base)
- This is stricter and trims positions where many copies have deletions
- Insertion hotspots (positions where copies have insertions relative to the reference) are tracked separately

### Fallback consensus (both scripts, on non-convergence)

- Triggered when bootstrap does not converge within max iterations
- Takes min(100, total copies) random copies
- Aligns with MAFFT, calls consensus with **plurality ≥ 35%** and **gaps in denominator** (min coverage 30%)
- Exits with code 0 and a warning message

## Input format

Standard FASTA file containing SINE copies. Can be:
- Raw copies (no consensus sequence) — e.g., extracted from a genome with RepeatMasker or BLAST
- `.bnk` files from repeat annotation pipelines
- Files with an existing consensus as the first sequence (will be treated as just another copy)

Minimum 3 sequences required.

## Dependencies

- [MAFFT](https://mafft.cbrc.jp/alignment/software/) — multiple sequence alignment
- [seqkit](https://bioinf.shenwei.me/seqkit/) — FASTA manipulation (pairwise script only)
- Standard Unix tools: `awk`, `shuf`, `bc`

## Example

```bash
# Quick bootstrap consensus
bash sine_consensus.sh copies.fasta
# → copies_consensus.fasta

# Pairwise consensus with 8 parallel jobs, keep intermediates
bash sine_pairwise_consensus.sh -j 8 -k copies.fasta
# → copies_pw_consensus.fasta, copies_pw_freqtable.tsv, etc.

# Pairwise with existing reference, relaxed threshold
bash sine_pairwise_consensus.sh -r reference.fa -t 35 copies.fasta
# → copies_pw_consensus.fasta

# Bootstrap with aggressive flank trimming
bash sine_consensus.sh -c 70 -g 30 copies.fasta
# → shorter but higher-confidence consensus
```
