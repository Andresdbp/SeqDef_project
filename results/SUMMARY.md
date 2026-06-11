# SeqDef revision — headline numbers (JEB-2026-00140)

Reproduced with `set.seed(42)`. The O(n) benchmarks and all figure builds ran on an Apple M1 Pro
laptop (16 GB, R 4.5.2); the dense O(n²) kernel benchmark ran on a TAMU Grace bigmem node
(3 TB, 80 cores, single-threaded, R 4.4.2). Every number below is produced by a script in
`analyses/` and stored in `results/*.csv`.

> **Dataset:** the pipeline prunes to **n = 877** species (MCC tips ∩ IUCN; **55/877** have an NCBI
> assembly). The earlier "n = 850" was a transcription error, now corrected.

---

## Baseline / robustness (Fig 3) — `posterior_priority.rds`
- 100 posterior trees; **_C. atromarginatus_ is the top priority in 85/100 trees**.
- Minority alternates are the CR guitarfish cluster (*Rhinobatos* / *Acroteriobatus* spp.).

## λ sensitivity + the two λ-selection methods → Fig S1 — `lambda_rank_stability.csv`, `method_comparison.csv`
- λ selected: **auto_max median 3.40 [2.65, 5.00]**; **by_genus median 3.19 [2.87, 4.87]**.
- **Two methods interchangeable:** Priority ρ(auto vs genus) = **0.998 [0.985, 1.000]**; agree on the
  top target **100/100**; each selects *C. atromarginatus* in 85/100.
- **Ranking is λ-robust:** Priority ρ vs auto_max = **0.93 median (0.88 min)** across λ∈[1,25];
  SeqDef ρ median ≈ 0.82. *C. atromarginatus* is top-1 in ~60–85 % of trees at every λ, peaking in the
  auto-selected band (λ ≈ 2.6–5).

## Kernel comparison (4 kernels × all 100 posterior trees) — `kernel_comparison.csv` (reported inline)
- ρ(Priority vs exponential): **gaussian 0.88, linear 0.94, brownian 0.78**.
- *C. atromarginatus* top target (of 100): **exp 85, gaussian 83, linear 87, brownian 23**.
- Median auto_max λ: exp 3.40, gaussian 1.80, linear 1.30 (Brownian parameter-free).
- **Message:** the tunable dense kernels track the exponential; the parameter-free Brownian reranks and
  recovers the top target far less often — quantifying why the tunable exponential is the default.
  (Former Supplementary Table S1; now inline in the manuscript text.)

## EDGE / EDGE2 benchmark → Fig S2 — `edge_benchmark.csv`, `edge_divergence.csv`
ED via `picante` fair-proportion; classic EDGE = log(1+ED) + GE·log2; EDGE2 = expected-PD-loss with the
50-yr IUCN→p_ext map (LC .0009 / NT .0071 / VU .0513 / EN .4276 / CR .9688; Gumbs et al. 2023).
- **ρ(EDGE, SeqDef) = 0.07** (near-orthogonal). ρ(EDGE, Priority) = **0.83**; ρ(EDGE2, Priority) = **0.79**
  (the shared IUCN threat term). Top-25 overlap: **9** (EDGE), **8** (EDGE2).
- **Worked divergence:** *Sphyrna lewini* (CR) ranks EDGE **40** but SeqDef **779/877** — its congener
  *S. mokarran* is already sequenced. *C. atromarginatus* ranks high on both (Centrophoridae has no genome).

## Computational cost across the four kernels → Fig S3 — `traversal_benchmark.csv`, `dense_kernel_benchmark.csv`, `dense_kernel_grace.csv`
- **Exponential & Brownian: exact O(n)** tree traversal (no n×n matrix; identical to the dense result,
  ≤ 1e-15). On the M1 Pro laptop a single call at **n = 10⁶ runs in < 1 s using ≈ 0.4 GB** (0.66 s / 0.39 GB;
  the two kernels are within ~10 %). *Runtime is reproducible; R peak-memory is ~30 % noisy run-to-run.*
- **Gaussian & linear: dense O(n²)** in time and memory (cophenetic + weight matrix + matrix–vector
  product). On the Grace 3 TB node, real separate runs reach **447 GB at n = 10⁵** and **~1 TB (1006 GB)
  at n = 1.5×10⁵**; runtime ≈ **726 s** (gaussian) / 622 s (linear) at n = 10⁵. Fitted log–log slopes
  ≈ **2.2** (time) and **2.0** (memory). `auto_max` adds a constant factor (O(k·n²), k ≤ 491).
- **Real Chondrichthyes (877 tips):** **0.02 s** (fixed λ) / **0.21 s** (auto_max); 100-tree posterior ≈ **40 s**.

## Brownian kernel — `brownian_comparison.csv` (reported inline)
- The **flat low-λ limit of the exponential**: ρ(Brownian, exp at λ = 0.5) = **0.98**.
- Less discriminating: SeqDef variance **0.038 vs 0.045** at the variance-maximizing λ.
- Reranks taxa (Priority ρ = **0.78** vs exponential), recovers the top target in only **23/100** trees,
  and favours the deep-diverging chimaera *Callorhinchus*. (Former Supplementary Table S2; now inline.)

---

## Package capability
`SeqDef(..., kernel = c("exponential","gaussian","linear","brownian"))` — default `"exponential"` reproduces
the previous behaviour (verified). The kernel is used in both the final calculation and the `auto_max`
variance loop, mirrored in `function.R` and the package source.

## Figures
- **Fig 1** (`fig1.pdf`) — the original **illustrative** toy panel ("Values shown are illustrative"; no λ
  reported). Kept illustrative by decision; response §2.6 reworded to match (it no longer claims Fig 1
  reports λ, and points to Fig 2 / Supplementary Fig S1 for the quantitative λ analysis).
- **Supplementary = three figures, no tables:** **S1** `figS1_lambda_stability.pdf` · **S2**
  `figS2_edge_vs_seqdef.pdf` · **S3** `figS3_runtime.pdf`. 300-DPI PNGs of all seven figures are in
  `figures/png/`. (The former Supplementary Tables S1/S2 are now inline text; the continuous-S figure was dropped.)
