# SeqDef revision — headline numbers (JEB-2026-00140)

Generated from `analyses/01–06`. Hardware: Apple Silicon (Darwin 25.5.0 arm64),
8 cores, **16 GB RAM**, R 4.5.2. All randomness seeded with `set.seed(42)`.
Tagged to the placeholders in `revision/SeqDef_response_to_reviewers.docx`.

> **Dataset note (flag for manuscript):** the reproducible pipeline prunes to
> **n = 877** species (MCC tips ∩ IUCN). The manuscript currently states **n = 850**.
> 55/877 species have an NCBI assembly. *Reconcile this number (likely an IUCN
> version difference) or state the filter that yields 850.*

---

## Baseline / robustness (reproduces Fig 3)
- 100 posterior trees; **C. atromarginatus is the top priority in 85/100 trees** (matches manuscript).
- Minority alternates: *Rhinobatos albomaculatus* (6), *R. jimbaranensis* (4),
  *Acroteriobatus variegatus* (2), *R. annandalei* (2), *R. rhinobatos* (2) — i.e. the
  CR guitarfish cluster, exactly as the manuscript narrative says.

## §1.1 / 2.3 — Continuous S (Analysis 6)  → `continuous_S_example.csv`, `figS_continuousS.pdf`
Illustrative continuous S from NCBI assembly count: **S = count/(count+1)** (no API call).
- Spearman ρ(SeqDef) binary vs continuous = **0.959**
- Spearman ρ(Priority) binary vs continuous = **0.981**; **top-10 overlap = 10/10**
- *C. atromarginatus* priority rank: binary **9 → continuous 9** (unchanged)
- λ: binary 4.80, continuous 4.20
- Redistribution toward species with only single-assembly relatives:
  *Carcharhinus falciformis* (738→389), *C. albimarginatus* (790→483), *Chimaera* spp. (~652→445)
- **Message:** high-priority region preserved; scores redistribute among species whose
  relatives are only partially covered. Confirms S need not be binary.

## §1.2 — λ sensitivity + method head-to-head (Analysis 1)  → `lambda_rank_stability.csv`, `method_comparison.csv`, `figS_lambda_stability.pdf`
- λ selected: **auto_max median 3.40 [2.65, 5.00]**; **by_genus median 3.19 [2.87, 4.87]**
- **Two methods are interchangeable:** Spearman ρ(Priority: auto vs genus) = **0.998 [0.985, 1.000]**;
  they agree on the top target in **100/100** trees; each selects *C. atromarginatus* in 85/100.
- **Ranking is λ-robust:** ρ(Priority vs auto_max) median **0.932** (min 0.878) across λ∈[1,25];
  ρ(SeqDef vs auto_max) median 0.818. Top-10 overlap plateaus ~8.3/10, peaks ~8.7 in the auto band.
- **Top target is the modal choice at every λ:** *C. atromarginatus* is top-1 in **65–80%** of
  sampled trees across the whole λ∈[1,25] range (never <50%), peaking at ~80% in the auto-selected
  band λ≈[2.6, 5]. (Honest framing: top-1 is the modal, auto-selected choice — not invariant at all λ.)

## §1.2 — Kernel comparison (Analysis 2)  → `kernel_comparison.csv`, `figS_kernel.pdf`
New `kernel=` arg (exponential/gaussian/linear), each with its own auto_max λ.
- ρ(Priority) vs exponential: **gaussian 0.882 [0.846, 0.961]**, **linear 0.943 [0.923, 0.958]**
- Target in top set (posterior, n=20): exp **16/20**, gaussian **15/20**, linear **16/20**
- **Message:** rankings are highly correlated across kernels and the top target is preserved;
  the exponential kernel is justified on principle (single interpretable half-life horizon).

## §1.3 — EDGE / EDGE2 benchmark (Analysis 3)  → `edge_benchmark.csv`, `edge_divergence.csv`, `figS_edge_vs_seqdef.pdf`
Computed on the MCC tree (ED via `picante` fair-proportion; classic EDGE = log(1+ED)+GE·log2;
EDGE2 = expected-PD-loss with 50-yr p_ext map LC=.0009/NT=.0071/VU=.0513/EN=.4276/CR=.9688, Gumbs et al. 2023).
- ρ(EDGE, Priority) = **0.832**; ρ(EDGE2, Priority) = **0.790**
- **ρ(EDGE, SeqDef) = 0.071** (essentially uncorrelated — they measure different things); ρ(ED, SeqDef) = 0.409
- Top-25 overlap: EDGE vs Priority **9/25**; EDGE2 vs Priority **8/25**
- *C. atromarginatus*: Priority rank 9, EDGE rank 43, EDGE2 rank 52, SeqDef rank 89 (CR) — **high on both**
- **Worked divergence (the value proposition):** *Sphyrna lewini* (CR, EDGE rank **40**) sits at
  SeqDef rank **779/877** because its congener **Sphyrna mokarran** already has an assembly — its
  marginal genomic info is already captured. Contrast: **Centrophorus (all 8 spp) and family
  Centrophoridae have ZERO assemblies**, so *C. atromarginatus* is high on both. (8 high-EDGE/
  low-SeqDef species total with a sequenced congener; also *Mobula mobular*←*M. hypostoma*,
  *Mustelus whitneyi*←*M. asterias*.)

## §1.4 — Runtime + memory + complexity (Analyses 4 + Grace) → `runtime_benchmark_grace.csv`, `figS_runtime.pdf`
Single λ (cost of ONE SeqDef call). rtree(n), random binary S, 100 reps where cheap (down to 3 at
the largest n). Benchmarked on a **TAMU Grace large-memory node (x86_64, 80 cores, 3 TB RAM),
single-threaded BLAS, R 4.4.2** — 36 points from **n = 100 to 100,000**.

| n | time (s) | peak memory |
|---|----------|-------------|
| 100 | 0.001 | 36 MB |
| 1,009 | 0.041 | 91 MB |
| 10,174 | 6.24 | 4.7 GB |
| 20,000 | 24.1 | 17.9 GB |
| 30,000 | 55.0 | 40 GB |
| 50,000 | 161 | 112 GB |
| 75,000 | 411 | 252 GB |
| 100,000 | **808** (13.5 min) | **447 GB** |
| 150,000 | — | OOM (exceeded the 800 GB allocation) |

- **Complexity: O(n²)** time & memory (cophenetic matrix + weight matrix + matrix–vector product);
  `auto_max` is O(k·n²), k ≤ 491 (λ grid `seq(1,50,0.1)` with early stop).
- **Empirical scaling confirms O(n²):** fitted log–log slope **2.09** (time) and **1.94** (memory)
  over n ≥ 1,000. Clean single-threaded curve across three orders of magnitude in n.
- **Memory is the binding constraint, not time.** Peak memory grows as O(n²): ~4.7 GB at n=10k,
  **447 GB at n=100k**; the run was pushed until memory exhaustion — n=150,000 exceeded the 800 GB
  allocation (OOM-killed), consistent with the fit. Extrapolating, n≈250k needs ≈2.8 TB (fits a 3 TB
  node) and whole-eukaryote scale (n~10⁶) would need tens of TB.
- **Real Chondrichthyes (877 tips):** single-λ **0.020 s**, auto_max **0.21 s**; 100-tree posterior ≈ **40 s**.
- **Revised prose:** efficient at clade-to-class scale and demonstrably tractable to n=100,000 on a
  large-memory node; the dense O(n²) matrix — not runtime — is the binding constraint at tree-of-life
  scale, motivating a sparse/truncated-kernel approximation (future work).
- *Note: the SLURM job writes the CSV only at the end; both 800 GB runs OOM-killed at n=150k before
  that write, so the data were recovered from the job's stdout log (`results/grace_runtime_18719664.out`)
  via `analyses/recover_grace_runtime.R`. A clean rerun would cap n at ≤120k or request ≥1.5 TB.*

## §1.4b — Exact O(n) algorithm for the exponential kernel (Analysis 7) → `traversal_benchmark.csv`, `figS_runtime.pdf`
The exponential kernel admits an **exact linear-time, linear-memory** algorithm. Because the
cophenetic distance is additive along tree paths and exp(−λd) = ∏ exp(−λ·edge/T), the weighted
availability A = W·s is computed by a two-pass sum-product tree traversal (Felsenstein-style) — **no
n×n matrix is built**. (Gaussian uses d² and linear uses 1−λd, neither of which factorizes, so they
keep the dense O(n²) path.)
- **Exact:** matches the dense result to ~1e-16 across `rtree`/`rcoal`, n=10–1000, λ=1–20, and the
  full `auto_max` pipeline; reproduces the Chondrichthyes MCC result identically (λ=4.80,
  *C. atromarginatus* rank 9).
- **Scaling (single λ, laptop):** n=10⁴ → 0.01 s / 49 MB; n=10⁵ → 0.06 s / 74 MB (dense = 447 GB);
  **n=10⁶ → 0.9 s / ~0.37 GB** (dense would need ~48 TB). ~13,000× faster and ~6,000× less memory at
  n=10⁵; reaches whole-tree-of-life scale on commodity hardware.
- **`auto_max` for the exponential kernel is now O(k·n)** (k = λ-grid length).
- **Implication for M4 / "why exponential":** the same memoryless/multiplicative property that makes
  the exponential the natural biological model (constant proportional decay of shared information) is
  exactly what makes it exactly computable in O(n). This turns the manuscript's "sparse-kernel
  approximation is future work" into "we provide an *exact* O(n) algorithm for the default kernel."
- Implemented on git branch `linear-time-exp` (merged into `revision-jeb`); helper `.seqdef_exp_avail`.

## Brownian-motion kernel (Analysis 8) → `brownian_comparison.csv`, `figS_brownian.pdf`
Added `kernel="brownian"`: a parameter-free phylogenetic correlation (shared ancestry under BM),
computed in O(n) via the three-point structure (exact vs dense `cov2cor(vcv())`, diff ~1e-15).
- **It is the flat limit of the exponential:** ρ(BM, exp λ=0.5) = 0.98; it sits at the low-λ,
  low-discrimination end of the variance-vs-λ curve (auto_max picks the peak, λ≈4.8).
- On Chondrichthyes it is **less discriminating** (var 0.038 vs 0.045 at auto_max) and selects a
  **different top target** (*Callorhinchus*, a deep-diverging chimaera) — a concrete illustration of
  why the tunable exponential is preferred.
- Efficiency nuance: BM is **also O(n)** (three-point), so the precise claim is "the exponential
  uniquely combines O(n) with a *tunable* horizon"; Gaussian/linear are O(n²). (Wording updated in
  SUMMARY §1.4b, the rebuttal, and the change-list.)

## §2.6 — Figure 1 regenerated with λ (Analysis 5)  → `figures/fig1.pdf`, `fig1_values.csv`
- Rebuilt from a REAL `SeqDef()` run on the 10-taxon toy tree (taxa10 = only sequenced tip; taxa9 = its sister).
- toy-tree auto_max **λ = 1.2**; figure shows three horizons: **λ = 1.2 (auto), 5, 15**.
- Horizon effect: *taxa9* SeqDef = **0.23 → 0.63 → 0.95** as λ = 1.2 → 5 → 15 (its sister taxa10 is
  sequenced, so it is "covered" only under a broad horizon; the signal vanishes as λ grows).

---

## New package capability
- `SeqDef(..., kernel = c("exponential","gaussian","linear"))` added (default exponential =
  identical to previous behaviour; verified). Used in both the final calculation and the `auto_max`
  variance loop. Mirrored in `function.R`; `man/SeqDef.Rd` re-documented.

## New / regenerated figures
- `figures/fig1.pdf` (regenerated, λ-annotated) — MAIN
- `figures/figS_lambda_stability.pdf`, `figS_kernel.pdf`, `figS_edge_vs_seqdef.pdf`,
  `figS_runtime.pdf`, `figS_continuousS.pdf` — SUPPLEMENTARY
