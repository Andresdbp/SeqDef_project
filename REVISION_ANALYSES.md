# SeqDef — analyses to run for the JEB major revision

**For:** Claude Code, working in this repo (`SeqDef_project/`) and the sibling package repo (`SeqDef/`).
**Goal:** Produce the empirical results that the reviewer/editor comments require but that cannot be written as text until run. The pure-text manuscript edits are handled separately (see `revision/SeqDef_revision_strategy.docx` and the chat handoff). Do **not** rewrite the manuscript here — just generate analyses, figures, and a numbers summary.

Manuscript: *SeqDef: An R package for phylogenetically weighted genomic prioritization across the tree of life* — JEB-2026-00140, Major Revision. See `context.md` for repo orientation.

---

## 0. Ground rules

- **Reproducibility:** set `set.seed(42)` anywhere randomness is used (`rtree`, sampling, simulated S). Record `sessionInfo()` and hardware (`Sys.info()`, `parallel::detectCores()`, RAM) into `results/SESSION.txt`.
- **Outputs:** create `results/` (CSV + cached `.rds`) and write figures to `figures/`. Cache expensive intermediates (100-tree × λ-grid loops) as `.rds` so reruns are fast.
- **One script per analysis** in a new `analyses/` folder, plus `analyses/run_all.R`. Name them `01_…`.R through `06_…`.R as below.
- **Do NOT commit credentials.** `analysis.R` currently contains real-looking NCBI Entrez and IUCN keys (also in git history). Do not echo them, do not re-commit them. If an analysis needs to re-query NCBI/IUCN, read keys from `Sys.getenv()` only.
- **Path gotcha:** `analysis.R` reads `data/raw/iucn_assessment_data.rds` and `data/raw/ncbi_assembly_data.rds`, but those files actually live in `data/` (no `raw/` subfolder); `chondrichthyes.nex` is in `data/`. Reconcile before running (either move the two `.rds` into `data/raw/` or fix the paths). Confirm all three load.
- **Backward compatibility:** the default behaviour of `SeqDef()` must not change. The only signature change is an additive `kernel=` argument defaulting to `"exponential"` (Analysis 2).
- **End product:** `results/SUMMARY.md` containing the headline numbers, written to slot directly into the bracketed placeholders in `revision/SeqDef_response_to_reviewers.docx` (placeholder tags noted per analysis below).

## Repo assets you already have

- `function.R` — local copy of `SeqDef()` and `calc_priority()`.
- `SeqDef/R/SeqDef.R`, `SeqDef/R/calc_priority.R` — the package source (+ `man/`, `NAMESPACE`).
- `analysis.R` — Chondrichthyes pipeline (IUCN + NCBI + 100 posterior trees → SeqDef → `calc_priority` → Fig 4). Commented toy-tree block at the bottom is the basis for Fig 1.
- `validation.R` — already computes variance-vs-λ curves (Fig 2A/B), the λ distributions for both selection methods, and the topological-robustness test (Fig 3). **Extend this; don't duplicate it.**
- `data/chondrichthyes.nex` — 100 posterior trees (VertLife/Stein et al. 2018).
- `data/iucn_assessment_data.rds`, `data/ncbi_assembly_data.rds` — IUCN categories and NCBI assembly presence/absence (binary `assembly_availability`).
- `seqdef.csv` — toy 10-taxon availability vector (taxa1–10, last = 1).

Key `SeqDef()` internals to know: the kernel is hard-coded at `dist.prop <- exp(-final_lambda * dist.matrix / td)`; `lambda="auto_max"` scans `seq(1,50,0.1)` maximizing variance (10% drop tolerance) and its internal `calc_var()` helper **also** hard-codes the exponential kernel; `lambda="by_genus"` sets λ = log(2)/median intra-genus normalized distance. `calc_priority()` computes `s_scores * (base ^ trait)`, i.e. SeqDef × 2^GE with GE = {LC0, NT1, VU2, EN3, CR4}.

---

## Analysis 1 — λ sensitivity on *rankings* + head-to-head of the two λ methods
**Resolves:** Reviewer 1 #2 (first half), Editor major #2a. **Placeholders:** response §1.2.

`validation.R` already shows variance vs λ and the λ distributions. What's missing is the reviewer's actual question — *does λ change the prioritization, not just the variance?* Add:

1. **Rank stability vs λ.** On the MCC tree (and a sample of ≥20 posterior trees for CIs), compute SeqDef and the final Priority over `lambda = seq(1, 25, 0.5)`. For each λ:
   - Spearman ρ of the SeqDef ranking and of the Priority ranking vs the `auto_max` solution;
   - top-10 overlap (count or Jaccard) vs the `auto_max` top-10;
   - whether the top-1 target is *Centrophorus atromarginatus*.
   Plot ρ and top-10 overlap against λ (mean + 95% band across trees).
2. **Method head-to-head.** For each of the 100 posterior trees compute Priority under `auto_max` and under `by_genus`; report Spearman ρ between the two Priority vectors and the fraction of trees whose top target agrees. (Reuse the λ-value distributions already in `validation.R`.)

**Outputs:** `results/lambda_rank_stability.csv`, `results/method_comparison.csv`, `figures/figS_lambda_stability.pdf`.
**Report:** median ρ and range across λ∈[1,25]; the λ band over which top-1 is invariant; auto_max-vs-by_genus ρ and % top-target agreement.

## Analysis 2 — Kernel comparison (requires a small package change)
**Resolves:** Reviewer 1 #2 (second half), Editor major #2b. **Placeholders:** response §1.2 (kernel sentence).

**Code change (mirror in BOTH `SeqDef/R/SeqDef.R` and `function.R`):** factor the kernel into one internal helper and add a `kernel = c("exponential","gaussian","linear")` argument (default `"exponential"`, so existing behaviour is unchanged). The helper must be used in *both* the final calculation and the `auto_max` `calc_var()` loop so optimization matches the chosen kernel. Use normalized distance `x = d / T_depth`:

```r
kern <- function(x, lambda, kernel) {
  switch(kernel,
    exponential = exp(-lambda * x),
    gaussian    = exp(-lambda * x^2),
    linear      = pmax(0, 1 - lambda * x)   # clamp at 0
  )
}
```

Update roxygen for the new arg and re-run `devtools::document()`; add it to `man/SeqDef.Rd`.

**Analysis:** on the Chondrichthyes data (MCC + a posterior sample), compute SeqDef and Priority under each kernel (each with its own `auto_max` λ). Report Spearman ρ of rankings across kernels and whether the top target is preserved.

**Outputs:** `results/kernel_comparison.csv`, `figures/figS_kernel.pdf`.
**Report:** cross-kernel rank correlations; top-1 preserved (yes/no per kernel).

## Analysis 3 — Benchmark vs EDGE / EDGE2
**Resolves:** Reviewer 1 #3, Editor major #3. **Placeholders:** response §1.3.

Use the MCC tree + the IUCN table.

1. **ED:** `picante::evol.distinct(tree, type = "fair.proportion")` (also do `"equal.splits"` as a check), or `caper`.
2. **Classic EDGE (Isaac et al. 2007):** `EDGE = log(1 + ED) + GE * log(2)`, with `GE = {LC0, NT1, VU2, EN3, CR4}`.
3. **EDGE2 (Gumbs et al. 2023) — preferred, optional if impractical:** implement per the published protocol (ED2 accounts for relatives' extinction probabilities; score combines ED2 with the species' extinction probability). Use an IUCN→pₑₓₜ mapping and **cite the source you use**; a standard 50-year set is `LC=0.0009, NT=0.0071, VU=0.0513, EN=0.4276, CR=0.9688`. If a full EDGE2 implementation isn't feasible in time, do classic EDGE (required) and note EDGE2 as future work.
4. **Compare** EDGE / EDGE2 against SeqDef-Priority: Spearman ρ, overlap of the top-25 lists, and a scatter (EDGE vs Priority) with the top target highlighted.
5. **The key result — explained divergence.** Find ≥1 species that is high-EDGE but low-SeqDef *because a close relative is already sequenced* (verify against `ncbi_assembly_data.rds`), and confirm *C. atromarginatus* ranks high on both (Centrophoridae has no genome). This pair is the evidence of complementarity.

**Outputs:** `results/edge_benchmark.csv` (per species: ED, EDGE, EDGE2, SeqDef, Priority, and each rank), `figures/figS_edge_vs_seqdef.pdf`, plus a small divergence table.
**Report:** ρ(EDGE, Priority), ρ(EDGE2, Priority); top-25 overlap; the worked divergence example(s).

## Analysis 4 — Runtime + memory benchmark and complexity statement
**Resolves:** Reviewer 1 #4, Editor major #4. **Placeholders:** response §1.4.

- Benchmark `SeqDef()` on `ape::rtree(n)` with random binary S for `n ∈ {100, 250, 500, 1000, 2500, 5000, 10000}`, ≥3 reps, median reported. Time **single-λ** (e.g. `lambda = 10`) and **`auto_max`** separately. Capture peak memory (`bench::mark()` gives `mem_alloc`, or `peakRAM`/`gc()`).
- Time the **real Chondrichthyes run**: 850 tips × 100 trees × `auto_max` (wall-clock).
- Plot runtime vs n on log-log; fit the slope (expect ≈ 2, confirming ~O(n²)).
- Record hardware + R version in `results/SESSION.txt`.

**Outputs:** `results/runtime_benchmark.csv`, `figures/figS_runtime.pdf`.
**Report:** runtime/memory table; log-log slope; Chondrichthyes wall-clock; one line confirming the O(n²)-memory ceiling (≈7 GB at n≈30k).

## Analysis 5 — Regenerate Figure 1 with λ annotated
**Resolves:** Reviewer 2 (Fig 1), Editor seconded. **Placeholders:** response §2.6.

The current Fig 1 toy panel uses **hand-set** SeqDef values (commented block in `analysis.R`), so there is no real λ to report. Rebuild it from an actual `SeqDef()` run on the 10-taxon toy tree (`seqdef.csv` / the `tree_text` in `analysis.R`), and annotate the λ used. Show **2–3 panels** at different λ (e.g. a small λ, the `auto_max` λ, and a large λ) so readers see the "phylogenetic horizon" effect the reviewer asked about. Keep the existing visual style.

**Output:** overwrite `figures/fig1.pdf`; note the λ value(s) for the caption in `results/SUMMARY.md`.

## Analysis 6 — Continuous-S worked example (supplementary)
**Resolves:** Reviewer 1 #1 / Reviewer 2 (S examples). **Placeholders:** response §1.1.

Demonstrate that a continuous S behaves sensibly. `ncbi_assembly_data.rds` currently holds only binary `assembly_availability`, so either:
- **(preferred)** re-query NCBI for assembly *level* and/or contig N50 and map to [0,1] (e.g. Complete=1.0, Chromosome=0.9, Scaffold=0.6, Contig=0.4; needs `ENTREZ_KEY` via env var), or
- **(fallback)** construct an illustrative continuous S from a plausible quality proxy and label it clearly as illustrative.

Run SeqDef with continuous S on the Chondrichthyes data and compare to the binary-S ranking (Spearman ρ; show the high-priority region is preserved while scores redistribute among low-quality assemblies).

**Outputs:** `results/continuous_S_example.csv`, `figures/figS_continuousS.pdf`.
**Report:** ρ(continuous, binary) ranking; one or two species whose priority shifts and why.

---

## Deliverables checklist

- [ ] `analyses/01…06_*.R` + `analyses/run_all.R`
- [ ] `kernel=` argument added to `SeqDef()` in **both** `SeqDef/R/SeqDef.R` and `function.R`; roxygen + `man/` updated; `devtools::document()` run
- [ ] `results/` CSVs + cached `.rds`; `figures/` updated (fig1 regenerated; new figS_* added)
- [ ] `results/SESSION.txt` (hardware + sessionInfo)
- [ ] `results/SUMMARY.md` — every headline number, tagged to the response-doc placeholders (§1.1, §1.2, §1.3, §1.4, §2.6)
- [ ] Confirm no API keys are printed or committed; data paths reconciled

## Out of scope here
Manuscript prose edits, the notation/equation fixes, Data Availability Statement, CRediT, key rotation, and unit tests are tracked elsewhere (`context.md` "Still to do"). This file is analyses only.
