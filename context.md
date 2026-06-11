# SeqDef project — context (for a fresh chat)

## What this is
Analysis/manuscript repo for **SeqDef**, a phylogenetically weighted genomic-prioritization statistic
+ R package. Manuscript *"SeqDef: An R package for phylogenetically weighted genomic prioritization
across the tree of life"* — **Journal of Evolutionary Biology, JEB-2026-00140**. The revision is
**final and pre-submission**: the science, numbers, scripts, and figures are frozen; what remains is
author-side compliance paperwork and one response-letter consistency issue (see flags).

Two local repos (both on `main`, both pushed to GitHub, **public**):
- Package: `/Users/andresbarboza/Documents/GitHub/SeqDef` → github.com/Andresdbp/SeqDef
- Project: `/Users/andresbarboza/Documents/GitHub/SeqDef_project` → github.com/Andresdbp/SeqDef_project

Manuscript + rebuttal live on **Google Drive** (read via the Drive connector; can't edit in place):
- "manuscript" Doc id `1DF6OUwHqvZFZgv40ngG5BTNBHvt7WZP773rATRwoYzY`
- "Answer to reviewers" Doc id `11-lz-Sv1BZSPetPTl99MzyfvryxcNLXj_StHhHsldfI`

## Package (done, committed, pushed)
- `SeqDef(tree, df, data.col, invert, scale, lambda, kernel)` with **kernel = exponential (default) /
  gaussian / linear / brownian**.
  - **exponential** → exact **O(n)** two-pass tree traversal (`.seqdef_exp_avail`), no distance matrix.
  - **gaussian, linear** → dense O(n²) cophenetic path.
  - **brownian** → parameter-free phylogenetic correlation (`.seqdef_bm_avail`), also **O(n)**; ignores λ.
  - All verified identical to the dense computation (≤1e-15). Default behaviour unchanged.
- `lambda` = `"auto_max"` (variance-max scan 1–50 in 0.1 steps, early stop) / `"by_genus"`
  (half-life = median intra-genus distance) / numeric.
- `calc_priority(seqdef_res, trait, model="exponential", base=2)` = SeqDef × base^trait.
- DESCRIPTION author = **Andres Barboza <andresdbp00@gmail.com>** (renders "Barboza, A.").
  README credits **Dr. Andres Barboza** + a citation line.
- NOTE: DESCRIPTION says `License: MIT + file LICENSE` but **no LICENSE file exists** (held per user;
  R CMD check will warn). `function.R` (project) mirrors the package SeqDef + calc_priority.

## Analyses (project `analyses/`, run from project root)
Scripts are now **named for the figure or in-text result they produce** (no more 01–08 numbering).
`00_setup.R` is the shared setup; `run_all.R` drives the suite. **NEVER source `analysis.R`** — it is
the raw data-fetch script and held the API keys (load the cached `.rds`/`.nex` directly instead).

- `00_setup.R` — shared paths, seeds (`set.seed(42)`), helpers.
- `figS1_lambda_stability.R` — λ rank-stability + the two-method head-to-head → **Fig S1**.
- `figS2_edge_benchmark.R` — EDGE / EDGE2 benchmark and divergence → **Fig S2**.
- `figS3_runtime.R` — runtime/memory plot across the four kernels (dense O(n²) vs O(n) traversals) → **Fig S3**.
- `kernel_comparison.R` — 4-kernel comparison; values reported **inline** in the manuscript (no figure).
- `brownian_comparison.R` — Brownian kernel characterisation; values reported **inline** (no figure).
- `grace_dense_kernels.R` + `grace_dense_kernels.slurm` — HPC dense O(n²) gaussian+linear runs on
  TAMU Grace, feeding the O(n²) curve in Fig S3.
- `run_all.R` — runs the local analysis scripts in order.

**DELETED scripts** (no longer in the repo): `04_runtime_benchmark`, `05_fig1_regen`,
`06_continuous_S`, `plot_runtime`, `recover_grace_runtime`, `grace_runtime_compute`,
`grace_runtime.slurm`. (Older `01`–`08` names are gone; if referenced anywhere they are stale.)

## Results (`results/`; headline numbers also in `SUMMARY.md`)
CSVs kept: `brownian_comparison`, `dense_kernel_benchmark`, `dense_kernel_grace`, `edge_benchmark`,
`edge_divergence`, `kernel_comparison`, `lambda_rank_stability`, `method_comparison`,
`traversal_benchmark` (+ cached `mcc_tree.rds`, `posterior_priority.rds`, `SESSION.txt`).
**DELETED CSVs:** `continuous_S_example`, `fig1_values`, `runtime_benchmark`, `runtime_benchmark_grace`.

### Headline results (reproduced with `set.seed(42)`)
- **n = 877** species (phylo∩IUCN; **55** have an NCBI assembly). *C. atromarginatus* is the top
  priority in **85/100** posterior trees.
- **λ selection:** auto_max median **3.40 [2.65, 5.00]**; by_genus median **3.19 [2.87, 4.87]**. The two
  methods give Priority **ρ = 0.998** and agree on the single top target in **100/100** trees.
- **λ rank-stability:** across λ∈[1,25], Priority ρ vs auto_max = **0.93 median, 0.88 minimum**;
  *C. atromarginatus* is top-1 in 65–80 % of trees at every λ, peaking in the auto-selected band (λ≈2.6–5).
- **Kernels (4 kernels × all 100 posterior trees):** ρ(Priority vs exponential) = gaussian **0.88**,
  linear **0.94**, brownian **0.78**. *C. atromarginatus* recovered as top target (of 100): exponential
  **85**, gaussian **83**, linear **87**, brownian **23**. Median auto_max λ: exp **3.40**, gaussian
  **1.80**, linear **1.30** (Brownian is parameter-free).
- **EDGE:** ρ(EDGE, SeqDef) = **0.07** (near-orthogonal); ρ(EDGE, Priority) = **0.83**;
  ρ(EDGE2, Priority) = **0.79**; top-25 overlap **9** (EDGE) and **8** (EDGE2). Worked divergence:
  *Sphyrna lewini* sits at EDGE rank **40** but SeqDef rank **779/877** (congener *S. mokarran* already
  sequenced); *C. atromarginatus* ranks high on both (family Centrophoridae has no assembly).
- **Efficiency:** exponential & Brownian are **O(n)** (n = 10⁶ in < 1 s, ~0.4 GB, on an Apple M1 Pro
  laptop, 16 GB). Gaussian & linear are dense **O(n²)**; benchmarked on a TAMU Grace bigmem node (3 TB,
  80 cores, single-threaded, R 4.4.2), the real gaussian+linear runs reach **447 GB at n = 1×10⁵** and
  **~1 TB (1006 GB) at n = 1.5×10⁵**. The 877-tip Chondrichthyes run: **0.02 s** (fixed λ) /
  **0.21 s** (auto_max); full 100-tree posterior ≈ **40 s**.
- **Brownian = the flat low-λ limit of the exponential** (ρ ≈ 0.98 vs exponential at λ ≈ 0.5). On
  Chondrichthyes it is less discriminating (SeqDef variance **0.038 vs 0.045** at the variance-maximising
  λ), reranks taxa (Priority ρ = **0.78**), recovers the top target in only **23/100** trees, and instead
  favours the deep-diverging chimaera *Callorhinchus*.

## Figures
- **Fig 1** (`figures/fig1.pdf`) is the **original illustrative panel** — *not* λ-regenerated; the
  caption states **"Values shown are illustrative."** (See the response-letter flag below.)
- **Fig 2–4** (`fig2.pdf`, `fig3.pdf`, `fig4.pdf`): λ-variance sensitivity, topological robustness
  bar chart, and the circular Chondrichthyes priority tree, respectively.
- **Supplementary = exactly THREE figures, ZERO tables:**
  - **Fig S1** `figures/figS1_lambda_stability.pdf` — λ rank-stability across 100 posterior trees.
  - **Fig S2** `figures/figS2_edge_vs_seqdef.pdf` — SeqDef-vs-EDGE orthogonality on the MCC tree.
  - **Fig S3** `figures/figS3_runtime.pdf` — computational cost across the four kernels (O(n) vs O(n²)).
  - The former **Supplementary Tables S1 (kernel) and S2 (Brownian) were REMOVED**; their values are
    now reported **inline in the manuscript text**.
- `figures/png/` holds **300-DPI PNGs of all 7 figures** (fig1–4 + figS1–S3). `figures/poster_fig3.pdf`
  is an extra poster variant of Fig 3.

## Key decisions
- **n = 877** (the 850 in the original draft was an error).
- **Exponential is the default**, justified four ways: (i) grounded in sequence evolution (shared
  identity decays exponentially with divergence; Jukes & Cantor 1969); (ii) interpretable, tunable
  half-life horizon; (iii) standard in distance-decay/sequence-novelty work (Pavoine et al. 2005;
  Marini et al. 2022); (iv) it **uniquely combines exact O(n)** time/memory *with a tunable* horizon
  (Brownian is O(n) but parameter-free; gaussian/linear are O(n²)).
- **Continuous-S demonstration DROPPED** from the manuscript (binary presence/absence is justified on
  its own; script and CSV removed).
- **Brownian kept** but reported **inline** (no figure) — it sharpens the "why exponential" argument.
- **Tables S1/S2 converted to inline text** (supplementary is figures-only).
- **Fig 1 reverted** to the original illustrative panel (no λ annotation).

## What remains (author action, outside the manuscript body)
- **Data Availability Statement, CRediT, Funding** are on the **separate title page** (journal format),
  not the manuscript body. Still to do there: **mint the Zenodo DOI** to fill its placeholder.
- **Rotate the NCBI + IUCN API keys** still present in public git history (commit `126814d`).
- `revision/manuscript_changes.md` tracks the change-list; `revision/SeqDef_response_to_reviewers_FINAL`
  is the point-by-point response.
- New refs (Jukes & Cantor 1969, Gumbs et al. 2023, Kembel et al. 2010, Soares 2023) are all present in
  the reference list (verified).

## Figure 1 / response §2.6 — ✅ resolved (Option B)
**Decision:** keep the original **illustrative** Fig 1 (no λ; "Values shown are illustrative"). The formal
response §2.6 has been reworded to match — it states Fig 1 is a schematic, points to the reproducible λ
results (Fig 2, Supplementary Fig S1, in-text λ medians 3.40 / 3.19), and offers a λ-annotated panel if the
editor prefers. The earlier "regenerated λ = 1.2/5/15" claim was removed, so the response and the
manuscript figure are now consistent. No manuscript Fig 1 change needed.

## ⚠️ Security
`analysis.R` had real NCBI + IUCN API keys. The working copy is scrubbed (keys read from env vars),
**but the keys are still in the public git history** (commit `126814d`). **The user must rotate/revoke
them at NCBI and IUCN.** Never echo or re-commit keys; never source `analysis.R` in analyses (load the
cached `.rds`/`.nex` directly).

## Git
Both repos on `main` (pushed). Branches include `main`, `revision-jeb`, `linear-time-exp`,
`brownian-kernel`. End commit messages with:
`Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## Env
- **Laptop:** Apple M1 Pro, **16 GB RAM**, **R 4.5.2** (Apple Silicon; Darwin arm64).
- **HPC:** TAMU **Grace bigmem** node (**3 TB**, **80 cores**, single-threaded, **R 4.4.2**), reached via
  `ssh grace` (module `gfbf/2024a R/4.4.2`; user keeps an SSH ControlMaster open).
- Pandoc available for md→docx. `REVISION_ANALYSES.md` holds the original analysis build spec.
