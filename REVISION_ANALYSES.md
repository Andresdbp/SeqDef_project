# SeqDef — revision analysis pipeline

**Repo:** `SeqDef_project/` (analysis + manuscript repo) and the sibling package repo `SeqDef/`.
**Manuscript:** *SeqDef: An R package for phylogenetically weighted genomic prioritization across the tree of life* — JEB-2026-00140 (Major Revision, now final pre-submission).

This file documents the analysis pipeline that backs the revised manuscript: one entry per current script, what each does, which reviewer/editor concern it resolves, its inputs and outputs, and how to run it. The pure-text manuscript edits are tracked separately (`revision/`). See `context.md` for repo orientation.

The supplementary materials are **three figures and no tables**:

| Figure | File | What it shows |
| --- | --- | --- |
| Fig S1 | `figures/figS1_lambda_stability.pdf` | Priority-ranking stability across λ (100 posterior trees) |
| Fig S2 | `figures/figS2_edge_vs_seqdef.pdf` | SeqDef-vs-EDGE orthogonality on the MCC tree (n = 877) |
| Fig S3 | `figures/figS3_runtime.pdf` | Computational cost across the four kernels (runtime + memory) |

The former Supplementary Tables (kernel agreement; Brownian comparison) were removed; those values are now reported **inline in the manuscript text**. The two scripts that produce them (`kernel_comparison.R`, `brownian_comparison.R`) write CSVs and print the in-text numbers, but no longer emit a figure or table.

---

## Ground rules

- **Reproducibility:** `set.seed(42)` is set in `00_setup.R` and in every script that draws random trees or samples. Hardware and `sessionInfo()` are written to `results/SESSION.txt` by the `write_session()` helper in `00_setup.R`.
- **No network, no credentials.** The pipeline never sources `analysis.R` (which historically held API keys). `00_setup.R` loads the cached data directly from `data/*.rds` and `data/chondrichthyes.nex`. There are no NCBI/IUCN calls in the analysis pipeline.
- **One script per output.** Scripts live in `analyses/` and are named to match the figure or value they produce. Heavy intermediates (the 100-tree × λ-grid loops, the MCC tree) cache to `results/*.rds`; delete those to force a clean recompute.
- **Backward compatibility.** `SeqDef()`'s default behaviour is unchanged. The only additive signature change is `kernel = c("exponential", "gaussian", "linear", "brownian")`, defaulting to `"exponential"`.

### Package surface (`SeqDef()` / `calc_priority()`)

- **`kernel`** — `"exponential"` (default; exact O(n) tree traversal, tunable horizon) / `"gaussian"` / `"linear"` (both dense O(n²)) / `"brownian"` (O(n), parameter-free, ignores λ; the flat low-λ limit of the exponential).
- **`lambda`** — `"auto_max"` (scans λ ∈ [1, 50] in 0.1 steps, maximizing SeqDef variance with a 10 % drop tolerance) / `"by_genus"` (λ = log 2 / median intra-genus normalized distance) / a numeric value.
- **`calc_priority()`** = SeqDef × `base ^ trait`, with `base = 2` and the Global Endangerment weights GE = {LC 0, NT 1, VU 2, EN 3, CR 4}.

The kernel is factored into one internal helper used by both the final calculation and the `auto_max` variance loop, mirrored in `SeqDef/R/SeqDef.R` and the local `function.R`. `function.R` also exposes the linear-time traversal helpers `.seqdef_exp_avail()` (exponential) and `.seqdef_bm_avail()` (Brownian) used by the runtime benchmark.

### Shared data and setup (`analyses/00_setup.R`)

Every analysis script begins with `source("analyses/00_setup.R")`. It loads the kernel-enabled `SeqDef()`/`calc_priority()` from `function.R`, reads the cached inputs, builds the IUCN table, and exposes the shared helpers:

- **Inputs:** `data/chondrichthyes.nex` (100 posterior trees, VertLife/Stein et al. 2018), `data/iucn_assessment_data.rds` (IUCN assessments), `data/ncbi_assembly_data.rds` (binary `assembly_availability`).
- **Derived:** `iucn_clean` (one row per species: name, IUCN category, order), `risk_index` (IUCN → GE weight), and `tree_mcc` — the maximum-clade-credibility tree, cached to `results/mcc_tree.rds`.
- **Helpers:** `build_input_binary()` (prune a tree to the IUCN intersection and attach binary availability), `ge_for_tree()` (GE vector aligned to a tree's tips), `clean_tree()`, and the `TARGET` constant (`Centrophorus_atromarginatus`, the case-study top target).
- **Case-study scope:** the phylogeny-∩-IUCN intersection is **n = 877 species**, of which **55 have an NCBI assembly**.

---

## Scripts

### `figS1_lambda_stability.R` — λ rank-stability + the two λ-selection methods
**Resolves:** Reviewer 1 #2a / Editor major #2a — *does λ change the prioritization, not just the variance?* and how do the two λ-selection methods compare on outcomes.

Runs `auto_max` and `by_genus` on all 100 posterior trees (cached), reproduces the topological-robustness distribution of top-1 winners, compares the two methods head-to-head, and sweeps a fixed-λ grid (λ ∈ [1, 25], step 0.5) on a 20-tree sample to measure ranking agreement against each tree's own `auto_max` solution.

- **Inputs:** `00_setup.R` (posterior trees, IUCN, NCBI).
- **Outputs:** `results/posterior_priority.rds` (cached per-tree runs), `results/lambda_rank_stability.csv`, `results/method_comparison.csv`, `figures/figS1_lambda_stability.pdf`.
- **Key results:** across λ ∈ [1, 25], Priority ρ vs `auto_max` = **0.93 median, 0.88 minimum**; *C. atromarginatus* is the top target across the sweep, peaking in the auto-selected band (λ ≈ 2.6–5). The two λ methods select similar values (`auto_max` median **3.40** [2.65, 5.00]; `by_genus` median **3.19** [2.87, 4.87]), their Priority outputs correlate at **ρ = 0.998**, and they agree on the single top target in **all 100/100 trees**.

### `figS2_edge_benchmark.R` — benchmark against EDGE / EDGE2
**Resolves:** Reviewer 1 #3 / Editor major #3 — position SeqDef relative to existing prioritization metrics. The framing is **complementarity**: EDGE ignores existing genomic data; SeqDef conditions on it.

Everything is computed on the MCC tree so ED, EDGE, EDGE2 and SeqDef share one topology. Computes Evolutionary Distinctiveness (`picante::evol.distinct`, fair-proportion and equal-splits), classic EDGE = log(1 + ED) + GE·log 2, and EDGE2 (expected-PD-loss formulation using the standard IUCN 50-yr extinction-probability map: LC 0.0009, NT 0.0071, VU 0.0513, EN 0.4276, CR 0.9688). It then correlates each against SeqDef-Priority, computes top-25 overlap, and isolates the divergence cases — high-EDGE / low-SeqDef species whose congener is already sequenced.

- **Inputs:** `00_setup.R` (MCC tree, IUCN, NCBI assembly presence/absence); `picante`.
- **Outputs:** `results/edge_benchmark.csv` (per species: ED, EDGE, EDGE2, SeqDef, Priority + each rank), `results/edge_divergence.csv`, `figures/figS2_edge_vs_seqdef.pdf`.
- **Key results:** ρ(EDGE, SeqDef) = **0.07** (near-orthogonal); ρ(EDGE, Priority) = **0.83**, ρ(EDGE2, Priority) = **0.79** (the shared IUCN threat term); top-25 overlap **9** (EDGE) and **8** (EDGE2). Worked divergence: *Sphyrna lewini* (CR) ranks EDGE 40 but SeqDef **779/877** because congener *S. mokarran* already has a genome (low marginal genomic gain); *C. atromarginatus* ranks high on both (Centrophoridae has no genome).

### `figS3_runtime.R` — computational cost across the four kernels
**Resolves:** Reviewer 1 #4 / Editor major #4 — scalability and the complexity statement.

Benchmarks all four kernels on `ape::rtree(n)` with random binary availability. The exponential and Brownian kernels use the exact O(n) traversal (timed on the laptop to **n = 10⁶**); the Gaussian and linear kernels build the full n × n distance matrix (O(n²), timed on the laptop to the memory limit and on an HPC node at larger n). Plots runtime and peak memory vs n on log-log axes with O(n) and O(n²) reference lines; the large-n O(n²) points are read back from the separately produced Grace CSV.

- **Inputs:** `function.R` (traversal helpers); `results/dense_kernel_grace.csv` (HPC points; produced by `grace_dense_kernels.R`, see below).
- **Outputs:** `results/traversal_benchmark.csv` (exponential + brownian, O(n)), `results/dense_kernel_benchmark.csv` (gaussian + linear, O(n²), laptop), `figures/figS3_runtime.pdf`.
- **Key results:** exponential and Brownian are **O(n)** — n = 10⁶ in **< 1 s using ≈ 0.4 GB** on an Apple M1 Pro laptop (16 GB); the case-study Chondrichthyes run (877 tips) is ~0.02 s at fixed λ / ~0.21 s with `auto_max`, and ~40 s for the full 100-tree posterior. Gaussian and linear scale as **O(n²)** (fitted log-log slopes ≈ 2.2 time / 2.0 memory), reaching **447 GB at n = 10⁵** and **~1 TB (1006 GB) at n = 1.5 × 10⁵** on the HPC node (runtime ≈ 726 s at n = 10⁵).

### `kernel_comparison.R` — kernel-agreement values (reported in text)
**Resolves:** Reviewer 1 #2b / Editor major #2b — justify defaulting to the exponential kernel by showing rankings are robust to the kernel choice.

Runs all four kernels (each with its own `auto_max` λ where applicable) on the MCC tree and all 100 posterior trees, then reports cross-kernel Priority correlations and whether the top target is preserved. **Produces a CSV only; the values are reported inline in the manuscript** (this replaced the former Supplementary Table S1).

- **Inputs:** `00_setup.R` (MCC + 100 posterior trees, IUCN, NCBI).
- **Outputs:** `results/kernel_comparison.csv`.
- **Key results (vs exponential, across the 100 posterior trees):** Priority ρ = **gaussian 0.88, linear 0.94, brownian 0.78**. *C. atromarginatus* recovered as the top target (of 100 trees) = exponential 85, gaussian 83, linear 87, **brownian 23**. Median `auto_max` λ = exponential 3.40, gaussian 1.80, linear 1.30 (Brownian is parameter-free).

### `brownian_comparison.R` — Brownian-kernel values (reported in text)
**Resolves:** Reviewer comments on alternative kernels — demonstrate empirically that the Brownian-motion kernel is the flat, low-discrimination limit of the exponential.

Compares the Brownian kernel against the exponential on the MCC tree: variance-vs-λ curve, agreement at small λ, and the Priority re-ranking. **Produces a CSV only; the values are reported inline in the manuscript** (this replaced the former Supplementary Table S2).

- **Inputs:** `00_setup.R` (MCC tree, IUCN).
- **Outputs:** `results/brownian_comparison.csv`.
- **Key results:** Brownian ≈ the flat low-λ limit of the exponential (ρ ≈ 0.98 vs exponential at λ ≈ 0.5). SeqDef variance 0.038 (Brownian) vs 0.045 at the variance-maximizing λ. The Brownian kernel re-ranks taxa (Priority ρ = 0.78 vs exponential), recovers the top target in only **23/100 trees**, and favours the deep-diverging chimaera *Callorhinchus*.

### `grace_dense_kernels.R` + `grace_dense_kernels.slurm` — HPC dense O(n²) benchmark
**Feeds:** Fig S3 (the large-n Gaussian/linear points).

The Gaussian and linear kernels build the full n × n distance matrix, so peak memory ~6·n²·8 bytes pushes past laptop RAM. This compute-only job runs them on a large-memory node and writes the CSV that `figS3_runtime.R` reads back. It runs single-threaded (BLAS threads pinned to 1) for a clean O(n²) scaling curve, over n ∈ {25 000, 50 000, 100 000, 150 000}.

- **Where:** a TAMU Grace bigmem node (3 TB, 80 cores, R 4.4.2). The SLURM script requests the `bigmem` partition with `--mem=1500G`; `function.R` must be present alongside the script in the run directory.
- **Inputs:** `function.R` (the kernel-enabled `SeqDef()`).
- **Output:** `results/dense_kernel_grace.csv` (n, time_s, peak_mb, kernel for gaussian + linear).
- **Run:** `sbatch grace_dense_kernels.slurm` (adjust the `cd` path and `R_LIBS_USER` to your scratch layout).

---

## How to run

From the project root (`SeqDef_project/`):

```sh
Rscript analyses/run_all.R
```

`run_all.R` sources, in order: `figS1_lambda_stability.R`, `figS2_edge_benchmark.R`, `figS3_runtime.R`, `kernel_comparison.R`, `brownian_comparison.R`. It never sources `analysis.R` and loads cached data directly. Heavy steps cache to `results/*.rds`; delete those to force a clean recompute.

The HPC dense-kernel points in Fig S3 are produced **separately** on a large-memory node (`grace_dense_kernels.R` via `grace_dense_kernels.slurm`) and read back from `results/dense_kernel_grace.csv`. `figS3_runtime.R` expects that CSV to be present; if it is missing, re-run the Grace job (or restore the committed CSV) before running the runtime script.

---

## Outputs at a glance

**Figures** (`figures/`): `figS1_lambda_stability.pdf`, `figS2_edge_vs_seqdef.pdf`, `figS3_runtime.pdf`. Manuscript Figures 1–4 (`fig1.pdf`–`fig4.pdf`) are produced elsewhere; `figures/png/` holds 300-DPI PNGs of all seven figures. **Note:** `fig1.pdf` is the original *illustrative* toy panel (caption: "Values shown are illustrative") — it is not regenerated by this pipeline.

**Results CSVs** (`results/`): `lambda_rank_stability.csv`, `method_comparison.csv`, `edge_benchmark.csv`, `edge_divergence.csv`, `traversal_benchmark.csv`, `dense_kernel_benchmark.csv`, `dense_kernel_grace.csv`, `kernel_comparison.csv`, `brownian_comparison.csv`. Plus `SESSION.txt` (hardware + `sessionInfo()`), `SUMMARY.md` (headline numbers), and cached `.rds` (`mcc_tree.rds`, `posterior_priority.rds`).

## Environment

- **Laptop:** Apple M1 Pro, 16 GB RAM, R 4.5.2 (all O(n) benchmarks and the figure builds).
- **HPC:** TAMU Grace bigmem node, 3 TB RAM, 80 cores, single-threaded, R 4.4.2 (the O(n²) dense-kernel benchmark only).

## Out of scope here

This file documents analyses only. Manuscript prose, the Data Availability Statement, CRediT contributions, a minted Zenodo DOI, and rotation of the API keys still present in git history are author-action items tracked elsewhere (`context.md`, `revision/`).
