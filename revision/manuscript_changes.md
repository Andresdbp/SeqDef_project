# Manuscript change-list — JEB-2026-00140 (SeqDef)

**Status: FINAL pre-submission.** This file now tracks only what is *left to do*. Reconciled against the
live manuscript (`/tmp/manuscript.txt`, 2026-06-10) and the analysis outputs in `results/`. The large
revision (Intro reframe, S generalization, four-kernel justification, EDGE benchmark, Brownian rewrite,
notation/λ fixes, n = 877) is **applied in the manuscript** — those FIND/REPLACE blocks have been removed
to avoid confusion. The numbers below were re-verified against the committed CSVs.

---

## 1. Applied (summary — no action)

These edits are present in the current manuscript and verified:

- **Intro / motivation** — complementary-information framing, "marginal information value" definition,
  Marini et al. and Soares et al. citations.
- **Formal definition** — *S* generalized with worked examples; *λ* glyph and the displayed equations
  (`A_f = Σ w_fi S_i`, `SeqDef_f = 1 − A_{f,norm}`, `w_fi = exp(−λ · d_fi / T_depth)`) render correctly.
- **Four-kernel paragraph** — exponential / gaussian / linear / brownian described with the four reasons
  for the exponential default; Jukes & Cantor (1969) cited.
- **R implementation** — kernel/λ-selection and the O(n) vs O(n²) complexity paragraphs; the
  "We timed a single SeqDef calculation…" benchmark paragraph; `Priority = SeqDef_median × 2^GE`.
- **Practical workflow** — binary-coding justification (single clean copy; the duplicate paragraph and the
  continuous-S sentence were removed); the GE-doubling sentence (2⁴ = 16×); the 877-tip timing
  (0.02 s fixed-λ / 0.21 s auto_max per tree; ≈ 40 s for the 100-tree posterior).
- **Sensitivity analysis** — rank-stability across λ, the auto_max-vs-by_genus head-to-head, and the
  cross-kernel comparison.
- **EDGE benchmark** — methods paragraph and the "complementary to EDGE" results subsection
  (ρ = 0.07; *S. lewini* vs *C. atromarginatus* contrast); Gumbs et al. (2023) and Kembel et al. (2010)
  cited.
- **Discussion** — Brownian sentence rewritten to the flat low-λ-limit framing; grafting/insertion
  extension paragraph; calibrated efficiency claim in the Conclusions.
- **Supplementary figures renumbered to the final three** — S1 = λ rank-stability,
  S2 = SeqDef-vs-EDGE, S3 = computational cost. The former Supplementary Tables S1 (kernel) and S2
  (Brownian) were **removed**; their values are now reported inline in the text.
- **References** — Jukes & Cantor, Gumbs et al., Kembel et al. all appear in the reference list.

**Decisions locked in:** exponential is the default kernel; the continuous-S demonstration was dropped
(binary coding justified on its own); Brownian was kept but reported inline (no figure); the old
Tables S1/S2 became inline text. See §3 for the one figure decision that is still open.

---

## 2. Remaining manuscript fixes (apply in the Drive doc)

Line numbers refer to the body paragraphs in `/tmp/manuscript.txt`.

**Status (verified against the 2026-06-11 manuscript): items 1–6, 9, and 12 are APPLIED.** Still pending:
**7** (tip set → {1,…,n}), **8** (ggplot → ggplot2), **10** (Fig 4 caption "with an inner"), **11** (IUCN
2019/2024 year). ⚠️ **Verify:** item 2 was applied but ¶289 now shows "**3.40x**" — should be "3.40"
(looks like a stray keystroke).

1. **Duplicated kernel-comparison sentence** (Sensitivity analysis). The Gaussian/linear result is stated
   twice: once ending the paragraph at ¶283 ("Recomputing under the Gaussian and linear kernels gave
   Priority rankings highly correlated with the exponential default (ρ = 0.88 and 0.94 respectively) and
   preserved the top target.") and again, more completely, in the next paragraph at ¶294 ("Recomputing
   under the Gaussian, linear, and Brownian kernels … ρ = 0.88, 0.94, and 0.78, respectively …").
   **Delete the first (¶283) instance** and keep the fuller ¶294 paragraph.
2. **Stray space before comma** (¶288): "variance-maximization median at 3.40 , and genus-scale at 3.19"
   → "…median at 3.40, and genus-scale at 3.19".
3. **Conclusions — broken efficiency sentence** (¶423): "It is computationally efficient for the default
   exponential kernel, with an exact linear-time, linear-memory tree traversal **lets it** scale to
   phylogenies of millions of tips…" is ungrammatical. Recast, e.g.: "It is computationally efficient: for
   the default exponential kernel an exact linear-time, linear-memory tree traversal **lets** it scale to
   phylogenies of millions of tips, while the Gaussian and linear kernels use a dense O(n²) computation
   suited to clade-to-class scale."
4. **Subject–verb agreement** (¶394, Practical Considerations): "these data-driven approaches consistently
   **selects** parameters" → "**select**".
5. **Stray space before period** (¶122, four-kernel paragraph): "is flat at small distances then falls
   sharply ." → "…falls sharply.".

**Additional fixes found in the final consistency sweep (manuscript vs the committed CSVs):**

6. **Benchmark paragraph (¶173, "We timed a single SeqDef calculation…") — three stale numbers** vs the
   current Grace run (`dense_kernel_grace.csv`) and `traversal_benchmark.csv`:
   - dense runtime "**808 s** at n = 100,000" → **≈ 726 s** (gaussian; the 447 GB memory on the same line is correct);
   - exponential "**74 MB** at n = 100,000" → **≈ 105 MB** (current traversal benchmark);
   - fitted memory slope "**1.9**" → **≈ 2.0** (the time slope ≈ 2.1–2.2 is fine).
   (Optional: add the new n = 150,000 / ~1 TB point, already shown in the Fig S3 caption.)
7. **Tip-set notation** (¶68): "{*i*,…,*n*}" → "{1,…,*n*}".
8. **Abstract** (¶17): "ape and **ggplot**" → "**ggplot2**".
9. **Results / EDGE** (¶338): "SeqDef **focus** on…" → "**focuses**"; and the comma splice
   "different things, EDGE looks…" → "different things: EDGE looks…".
10. **Figure 4 caption** (¶347): "Circular Chondrichthyes tree**,** an inner…" → "…tree **with** an inner…".
11. **IUCN 2019 reference**: the lead year "2019" disagrees with the title's "…Threatened Species **2024**…";
    reconcile the year with the in-text "(IUCN, 2019)".
12. **Fig S2 — *C. borneensis* label removed** from the regenerated figure (it remains an unlabeled
    purple "divergent" point). In the manuscript, also drop *C. borneensis* from the **Fig S2 caption** so
    the figure, caption, and Results body (¶341–343) all name only *S. lewini*.

*Items 2–12 are surface fixes; none change a number or a scientific conclusion.*

**R2.4 (taxonomic-uncertainty weight) — leave as-is (decision, 2026-06-11).** No manuscript or response
change. Response §2.4 already frames this as a *reply-only* acknowledgement ("flagged as a planned option
… *rather than altering the core statistic*"), so it does not claim a manuscript sentence and there is no
contradiction. *Note for any follow-up:* the 100-posterior-tree design handles **topological** uncertainty,
which is distinct from the **taxonomic-placement** uncertainty the reviewer meant — so the rebuttal should
not lean on "multiple trees address it" (the current §2.4 wording does not).

---

## 3. Figure 1 — λ decision  ✅ RESOLVED (Option B)

**Decision:** keep the **original illustrative Fig. 1** (caption "Values shown are illustrative."; no λ).
Response §2.6 has been **reworded to match** (`SeqDef_response_to_reviewers_FINAL.md`): it now states
plainly that Fig. 1 is a schematic that does not report a λ, points the reviewer/editor to the
quantitative λ results that *are* reproducible (Fig. 2 variance-vs-λ + auto-selected band; Supplementary
Fig. S1 rank-stability; the in-text λ medians 3.40 / 3.19), and offers a λ-annotated panel if the editor
prefers. The earlier "regenerated λ = 1.2 / 5 / 15 Figure 1" claim was removed, so the response and the
manuscript figure are now consistent.

*Residual (accepted):* the editor's literal "report the λ on Fig. 1" is answered by reasoned explanation
plus an offer, not by the figure itself. **No manuscript change is required** — the Fig. 1 caption already
reads "Values shown are illustrative."

---

## 4. Compliance (author action, outside the manuscript body)

- **Data Availability Statement, CRediT, and Funding** live on the **separate title page** (per the
  journal format), not the manuscript body — ✅ handled there. Still to do within the DAS: **mint the
  Zenodo DOI** for the code snapshot (cite VertLife / Stein et al. 2018 for the trees; IUCN Red List data
  are not redistributable — give Red List version + query + access date; deposit the NCBI table + scripts).
- **Rotate the NCBI + IUCN API keys.** They were scrubbed from `analysis.R` at `HEAD` (commit `2f639a8`,
  env vars now) but **remain in public git history** (commit `126814d`, "Main push"). Rotating the keys is
  the only real remediation; rewriting history is optional.

---

## Appendix — verified key numbers (for reference)

Reproduced with `set.seed(42)`; sources in parentheses. Nothing here is an open item — included so the
inline manuscript values can be checked against the committed CSVs.

- **Dataset:** n = 877 species (phylo ∩ IUCN); 55 have an NCBI assembly. *C. atromarginatus* is the top
  priority in **85/100** posterior trees.
- **λ selection** (`results/method_comparison.csv`): auto_max median **3.40** [2.65, 5.00];
  by_genus **3.19** [2.87, 4.87]; Priority ρ = **0.998**; agree on top target **100/100**.
- **λ rank-stability** (`results/lambda_rank_stability.csv`): Priority ρ vs auto_max across λ ∈ [1, 25]
  ≈ **0.93** median, **0.88** minimum; *C. atromarginatus* top in ~65–85 % of trees, peaking near the
  auto-selected band.
- **Kernels, 100 posterior trees** (`results/kernel_comparison.csv`): ρ(Priority vs exponential) =
  gaussian **0.88**, linear **0.94**, brownian **0.78**; top-target recovery (of 100) = exp **85**,
  gaussian **83**, linear **87**, brownian **23**; median λ = exp **3.40**, gaussian **1.80**,
  linear **1.30** (Brownian parameter-free).
- **EDGE** (`results/edge_divergence.csv`, `results/edge_benchmark.csv`): ρ(EDGE, SeqDef) = **0.07**;
  ρ(EDGE, Priority) = **0.83**; ρ(EDGE2, Priority) = **0.79**; top-25 overlap **9** and **8**;
  *Sphyrna lewini* EDGE rank 40 but SeqDef rank 779/877 (congener *S. mokarran* already sequenced).
- **Brownian** (`results/brownian_comparison.csv`): flat low-λ limit of the exponential (ρ ≈ **0.98** vs
  exponential at λ ≈ 0.5); SeqDef variance **0.038** vs **0.045** at the variance-maximising λ; favours the
  deep-diverging chimaera *Callorhinchus*.
- **Efficiency.** Exponential and Brownian are O(n) (`results/traversal_benchmark.csv`, M1 Pro laptop,
  16 GB, R 4.5.2): n = 10⁶ in **0.66 s** / **0.62 s** using ≈ **0.4 GB**; n = 10⁵ in **0.06 s** /
  **0.05 s** using ≈ **0.1 GB**. The R peak-memory reading is noisy (~30 % run-to-run), so memory is
  reported as "≈ 0.1 / ≈ 0.4 GB" rather than a single value. Gaussian and linear are dense O(n²)
  (`results/dense_kernel_grace.csv`, TAMU Grace bigmem, 3 TB / 80 cores / single-threaded / R 4.4.2):
  **447 GB** at n = 10⁵, **~1 TB** (1006 GB) at n = 1.5 × 10⁵.

### Repository state (final)

- **Figures** (`figures/`): `fig1.pdf` (illustrative — see §3), `fig2–fig4.pdf`,
  `figS1_lambda_stability.pdf`, `figS2_edge_vs_seqdef.pdf`, `figS3_runtime.pdf`; 300-DPI PNGs of all seven
  in `figures/png/`.
- **Analysis scripts** (`analyses/`): `00_setup.R`, `figS1_lambda_stability.R`, `figS2_edge_benchmark.R`,
  `figS3_runtime.R`, `kernel_comparison.R` (in-text values), `brownian_comparison.R` (in-text values),
  `grace_dense_kernels.R` + `grace_dense_kernels.slurm` (HPC dense O(n²) run feeding Fig S3), `run_all.R`.
- **Results CSVs kept** (`results/`): `brownian_comparison`, `dense_kernel_benchmark`, `dense_kernel_grace`,
  `edge_benchmark`, `edge_divergence`, `kernel_comparison`, `lambda_rank_stability`, `method_comparison`,
  `traversal_benchmark`.
- **New references in the manuscript:** Jukes & Cantor (1969); Gumbs et al. (2023); Kembel et al. (2010).
  Soares et al. (2023) was already cited.
