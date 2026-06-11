# SeqDef

**SeqDef** is a phylogenetically aware statistic, and a companion R implementation, for
prioritizing taxa for genome-scale sequencing. For each tip on a phylogeny it returns a
single **sequencing-deficiency** score that quantifies the *marginal* genomic information
gained by sequencing that taxon, given the sequence resources already available in its
relatives. A taxon is scored high when it is phylogenetically isolated from existing
genomes and low when close relatives are already sequenced; the score is therefore
*dynamic* and updates as relatives are sequenced.

Formally, SeqDef is the complement of phylogenetically weighted data availability, where
the weights come from a distance-decay kernel applied to the cophenetic distances on a
time-scaled tree. Scores can be combined with any user-supplied tip trait (extinction
risk, economic or vector/pest status, etc.) via `calc_priority` to produce a conservation-
or objective-aware priority ranking.

This repository contains the R implementation, the analysis pipeline, and all data,
results, and figures for the accompanying manuscript.

## Manuscript

*SeqDef: An R package for phylogenetically weighted genomic prioritization across the
tree of life* — manuscript **JEB-2026-00140** (*Journal of Evolutionary Biology*), final
pre-submission.

> **Citation:** TBD — to be added on acceptance.

## The statistic and R package

The package exposes two main functions, defined in [`function.R`](function.R):

- **`SeqDef(tree, df, lambda = "auto_max", kernel = "exponential", ...)`** — takes an `ape`
  `phylo` tree and a two-column data frame (tip label, availability score *S* in `[0, 1]`;
  e.g. binary genome presence/absence) and returns the pruned tree, per-tip SeqDef scores,
  the availability data, and the λ used.
- **`calc_priority(seqdef_res, trait_values, ...)`** — combines SeqDef scores with a tip
  trait. The default model mirrors EDGE: `Priority = SeqDef * base^trait` (e.g.
  `SeqDef × 2^GE`, where `GE` is the IUCN risk weight 0–4).

### Kernels

The framework admits any kernel; the package implements four:

| Kernel | Form | λ | Complexity |
| --- | --- | --- | --- |
| **exponential** (default) | `exp(-λ·d)` | tunable (half-life / "phylogenetic horizon") | exact **O(n)** tree traversal — no `n × n` matrix is ever formed |
| gaussian | `exp(-λ·d²)` | tunable | dense **O(n²)** |
| linear | `max(0, 1 - λ·d)` | tunable | dense **O(n²)** |
| brownian | Brownian-motion phylogenetic correlation | parameter-free (ignores λ) | **O(n)** |

The exponential kernel is the default because it is grounded in standard sequence
evolution, gives an interpretable tunable half-life, is standard in distance-decay
biodiversity measures, and uniquely combines an **exact O(n)** algorithm with a tunable
horizon. The parameter-free Brownian kernel is the flat, low-λ limit of the exponential
(≈ λ 0.5) and is therefore less discriminating.

### Choosing λ

`lambda` may be a numeric value or one of two automatic selectors:

- **`"auto_max"`** (default, recommended) — scans λ and picks the value that maximizes the
  variance of SeqDef scores (maximum discriminatory power).
- **`"by_genus"`** — anchors the phylogenetic half-life to the median pairwise distance
  between congeneric species (a biologically interpretable genus-level horizon).

## Repository layout

```
function.R        SeqDef + calc_priority (the package; start here)
analysis.R        end-to-end Chondrichthyes data pipeline (uses NCBI/IUCN APIs)
analyses/         revision analysis scripts (see "Reproducing", below)
results/          cached result tables (CSV) and intermediate objects (RDS)
figures/          manuscript figures (PDF) + figures/png/ (300-DPI PNGs of all 7)
data/             input data (phylogeny .nex, IUCN and NCBI assembly .rds)
revision/         response to reviewers and manuscript change-list
```

Other top-level files: `validation.R` (unit checks of the statistic), `context.md` and
`REVISION_ANALYSES.md` (working notes), and `results/SUMMARY.md` (headline numbers).

## Reproducing the analyses

From the project root:

```sh
Rscript analyses/run_all.R
```

`run_all.R` runs, in order: `figS1_lambda_stability.R` (Fig. S1 + λ/method tables),
`figS2_edge_benchmark.R` (Fig. S2 + EDGE tables), `figS3_runtime.R` (Fig. S3),
`kernel_comparison.R`, and `brownian_comparison.R` (the last two report values quoted
inline in the manuscript). It loads cached data from `results/*.rds` and never re-runs
`analysis.R`, which hits external APIs. `00_setup.R` documents the upstream data assembly.
All randomized steps use `set.seed(42)`.

**Fig. S3 dense points (HPC):** the O(n²) Gaussian/linear benchmark points that require a
large-memory machine are produced separately by `analyses/grace_dense_kernels.R` (submitted
via `analyses/grace_dense_kernels.slurm`) and read back from
`results/dense_kernel_grace.csv`. The O(n) exponential/Brownian points run on a laptop and
are produced by `figS3_runtime.R`. These dense runs are **not** part of `run_all.R`.

## Case study and key results

Applied to **Chondrichthyes** (sharks, rays, chimaeras): n = 877 species at the
intersection of the VertLife posterior phylogeny and the IUCN Red List (55 with an NCBI
assembly), availability scored as binary genome presence/absence, aggregated over 100
posterior trees.

- **Top priority:** *Centrophorus atromarginatus* (Dwarf Gulper Shark), recovered as the
  single top target in **85 / 100** posterior trees under the exponential kernel.
- **λ selection is stable:** `auto_max` median λ = 3.40; `by_genus` = 3.19; the two agree on
  the top target in 100/100 trees (Priority ρ = 0.998). Across λ ∈ [1, 25] the Priority
  ranking is highly stable (median ρ = 0.93, min 0.88 vs `auto_max`).
- **Kernels largely agree:** Priority ρ vs the exponential default = 0.94 (linear),
  0.88 (gaussian), 0.78 (brownian); the top target is recovered in 87, 83, and 23 of 100
  trees respectively.
- **SeqDef is a distinct axis:** near-orthogonal to EDGE (ρ = 0.07); SeqDef-Priority
  correlates with EDGE (ρ = 0.83) and EDGE2 (ρ = 0.79) through the shared IUCN term, with
  top-25 overlaps of 9 and 8.
- **Efficiency:** the exponential and Brownian kernels are O(n) (n = 10⁶ in < 1 s, ≈ 0.4 GB
  on an Apple M1 Pro laptop); the dense Gaussian/linear kernels reach ≈ 447 GB at n = 10⁵
  and ≈ 1 TB at n = 1.5 × 10⁵ on a large-memory HPC node. The 877-tip run is 0.02 s
  (fixed λ) / 0.21 s (`auto_max`); the full 100-tree posterior ≈ 40 s.

See [`results/SUMMARY.md`](results/SUMMARY.md) for the full set of reproduced numbers.

## Figures

Manuscript figures live in `figures/` (PDF) with 300-DPI PNGs in `figures/png/`.
Supplementary material is **three figures, no tables**:

- **Fig. S1** — λ rank-stability across the posterior (`figS1_lambda_stability.pdf`)
- **Fig. S2** — SeqDef-vs-EDGE orthogonality (`figS2_edge_vs_seqdef.pdf`)
- **Fig. S3** — computational cost across the four kernels (`figS3_runtime.pdf`)

Figure 1 (`fig1.pdf`) is an illustrative toy-tree panel ("Values shown are illustrative").

## Computing environment

- **Laptop:** Apple M1 Pro, 16 GB RAM, R 4.5.2 (O(n) kernels, primary analyses).
- **HPC:** TAMU Grace bigmem node, 3 TB RAM, 80 cores (single-threaded), R 4.4.2 (dense
  O(n²) Fig. S3 points).
