# Manuscript change-list — JEB-2026-00140 (SeqDef) — REMAINING edits only

Reconciled against the live Google Doc ("manuscript", 2026-06-03). Edits you have **already applied**
are not repeated here (Intro reframe + Soares + "marginal information value"; S generalization with
worked examples; the four-kernel justification paragraph; the R-implementation kernel/λ + complexity
paragraphs; n = 877; the binary-coding justification). What remains is below, top-to-bottom.

**Decision applied here:** the continuous-S demonstration and its figure are **dropped** — the binary
coding is justified on its own (presence/absence is a legitimate, simple availability metric), and the
Methods "Defining S" examples already answer the reviewers' "examples beyond 0/1" request.

**Supplementary figures are renumbered S1–S5** (continuous-S removed):
S1 = `figS_lambda_stability`, S2 = `figS_kernel`, S3 = `figS_edge_vs_seqdef`, S4 = `figS_runtime`, S5 = `figS_brownian`.

Global: set **SeqDef / Priority** and label subscripts (norm, median) in upright (roman) type; italics
only for variables (*S, w, d, A, λ, n, i, f*). New references to add via Paperpile: **Jukes & Cantor 1969**
(now cited in the kernel paragraph), and **Gumbs et al. 2023** + **Kembel et al. 2010** (for the EDGE
benchmark below). Entries in the "New references" block.

---

## 1. Notation / dropped-λ fixes (Methods — Formal definition)

**1a. Kernel equation — restore λ.**
- FIND: "scaled by , a user-specified decay parameter that controls the phylogenetic scale of influence:"
  → "scaled by *λ*, a user-specified decay parameter that controls the phylogenetic scale of influence:"
- FIND (displayed eq): `wfi=exp(-dfi Tdepth)` → **w_fi = exp(−λ · d_fi / T_depth)**
- FIND: "based on the cophenetic distance dji between tips" → "…distance *d*_fi between tips" (index typo j→f).

**1b. Availability & SeqDef equations.**
- FIND `Af=i=1nwfisi` → **A_f = Σ_{i=1}^{n} w_fi S_i**
- FIND `SeqDeff=1-Af,norm` → **SeqDef_f = 1 − A_{f,norm}** ("SeqDef" roman, "norm" roman subscript).

**1c. Restore λ where the glyph was lost** (it renders as a blank "()"):
- FIND "the tunable decay parameter lambda ()," → "the tunable decay parameter *λ*,"
- In the four-kernel paragraph, restore λ in the formulas: `exp(-d)` → **exp(−λd)**; `exp(-d2)` → **exp(−λd²)**; `max(0,1-d)` → **max(0, 1−λd)**.
- (Discussion, later) FIND "the choice of the decay parameter () dictates" → "the choice of the decay parameter *λ* dictates"
- (Discussion, later) FIND "anchor  to the genus-level phylogenetic half-life" → "anchor *λ* to the genus-level phylogenetic half-life"

**1d. Figure-number fixes in the (already-applied) kernel paragraph** (renumbering S1–S5):
- FIND "(Supplementary Fig. S7)" → "(Supplementary Fig. S5)"  *(Brownian)*
- FIND "(Supplementary Fig. S3)" → "(Supplementary Fig. S2)"  *(kernel comparison)*

*Optional polish:* the "choice of the weight matrix W … we default to the exponential …" paragraph and
the new four-kernel paragraph now overlap on "why exponential." Consider trimming the older paragraph to
the W-flexibility point and letting the new paragraph carry the justification.

---

## 2. Methods — R implementation

**2a. Priority equation (Eq. 4).**
- FIND `Priority = SeqDefmedian2GE` → **Priority = SeqDef_median × 2^GE** ("Priority", "SeqDef" roman; "median" roman subscript; *GE* italic).

**2b. Add the "Computational benchmarks" subsection** (referenced by the kernel paragraph and the
complexity paragraph; numbers from `results/SUMMARY.md`). INSERT after the complexity paragraph:

> "**Computational benchmarks.** We timed a single SeqDef calculation on simulated trees (`ape::rtree`,
> random binary S). The dense path (Gaussian/linear kernels) scales as O(*n*²) in time and memory —
> fitted log–log slopes ≈ 2.1 and 1.9 — reaching 447 GB of peak memory at *n* = 100,000 on a 3 TB,
> 80-core node (R 4.4.2, single-threaded). For the default exponential kernel the exact O(*n*) traversal
> removes this wall: it returns identical scores (≤ 1 × 10⁻¹⁵) but runs *n* = 100,000 in 0.06 s using
> 74 MB, and *n* = 1,000,000 in < 1 s using ≈ 0.4 GB, on a laptop (Supplementary Fig. S4). The 877-tip
> Chondrichthyes analysis ran in 0.02 s (fixed λ) / 0.21 s (auto_max) per tree; the full 100-tree
> posterior completed in ≈ 40 s."

---

## 3. Methods — Practical Workflow (clean up the binary/continuous text)

The binary-coding text was **pasted twice** and still contains the continuous-S sentence. Fix:

1. **Keep** one clean sentence-pair: "*We defined availability (S) as binary (1 = presence of any assembly
   in NCBI GenBank; 0 = absence). For this demonstration, the binary coding allows us to have an objective
   metric in presence/absence, which is obtained uniformly for all species from a single NCBI query,
   avoids confounding the demonstration with heterogeneous and frequently missing quality metrics, and
   matches the Earth BioGenome framing of whether a reference genome yet exists.*"
2. **Delete** the trailing continuous-S sentence in that paragraph: "To show that the statistic behaves
   sensibly with a continuous input … minimal coverage." (and the stray `"`).
3. **Delete the entire duplicated paragraph** that begins "For this demonstration the binary coding is
   deliberate: assembly presence/absence is objective …" through "… minimal coverage."

**3a. Doubling sentence.**
- FIND: "This formula creates a multiplicative effect, assigning the highest priority to taxa that are both phylogenetically isolated from existing genomes and critically threatened."
- REPLACE: "Because GE enters as an exponent of 2, each one-category increase in extinction risk **doubles** a taxon's priority: a Critically Endangered taxon (GE = 4) receives 2⁴ = 16× the weight of a Least Concern taxon (GE = 0) with the same SeqDef. Priority therefore rises steeply for taxa that are simultaneously isolated from existing genomes and highly threatened."

---

## 4. Methods — Sensitivity analysis (add rank-stability + method head-to-head)

After the existing variance/λ paragraphs, INSERT:

> "Beyond the variance of the statistic, we asked whether λ changes the *prioritization*. Across
> λ ∈ [1, 25] the Priority ranking was highly stable (Spearman ρ vs the auto_max solution = 0.93 median,
> 0.88 minimum), and *Centrophorus atromarginatus* was the top target in 65–80 % of sampled trees at
> every λ, peaking in the auto-selected band (λ ≈ 2.6–5; Supplementary Fig. S1). We also compared the two
> λ-selection algorithms directly across the 100 posterior trees: they selected similar values
> (variance-maximization median 3.40 [2.65, 5.00]; genus-scale 3.19 [2.87, 4.87]), their Priority outputs
> correlated at ρ = 0.998, and they agreed on the single top target in all 100 trees. We recommend
> variance-maximization as the default and genus-scale where a biologically interpretable horizon is
> preferred. Recomputing under the Gaussian and linear kernels gave Priority rankings highly correlated
> with the exponential default (ρ = 0.88 and 0.94) and preserved the top target (Supplementary Fig. S2)."

**4a. (Optional) Figure 2 caption** — add: "Across the 100 posterior trees the two algorithms selected
similar λ (medians 3.40 and 3.19) and near-identical prioritizations (ρ = 0.998), agreeing on the top
taxon in all 100 trees."

---

## 5. Methods + Results — EDGE benchmark (new)

**5a. Methods** — INSERT a short paragraph (end of the workflow or sensitivity):

> "To position SeqDef against established metrics, we benchmarked it on the Chondrichthyes MCC tree
> against EDGE (Isaac et al., 2007) and EDGE2 (Gumbs et al., 2023). We computed Evolutionary
> Distinctiveness (fair-proportion; picante, Kembel et al., 2010), formed classic EDGE = log(1 + ED) +
> GE·log 2 and EDGE2 (expected PD loss with the standard IUCN→50-yr extinction-probability map), and
> compared each with SeqDef-Priority by Spearman correlation, top-25 overlap, and the species at which
> they most diverge."

**5b. Results** — INSERT a subsection after "Chondrichthyes Sequencing Recommendation":

> "**SeqDef is complementary to EDGE.** SeqDef and EDGE were near-orthogonal on the Chondrichthyes MCC
> tree (Spearman ρ = 0.07): they measure different things — evolutionary distinctiveness plus threat
> versus genomic-data deficiency. SeqDef-Priority correlated with EDGE (ρ = 0.83) and EDGE2 (ρ = 0.79)
> only through the shared IUCN threat term, and the top-25 lists overlapped by 9 and 8 species. The
> informative cases are the divergences: the scalloped hammerhead *Sphyrna lewini* (Critically
> Endangered) ranks high on EDGE (rank 40) but low on SeqDef (rank 779/877) because its congener
> *S. mokarran* already has a genome, so its marginal genomic information is largely captured. By
> contrast, *C. atromarginatus* ranks high on both because the family Centrophoridae has no assembly at
> all. SeqDef thus adds a distinct, dynamic axis — marginal genomic information given the data in hand —
> that updates as relatives are sequenced (Supplementary Fig. S3)."

---

## 6. Discussion

**6a. Fix the Brownian sentence** (it is now implemented — the current text is incorrect).
- FIND: "Additionally, we want to allow the users the ability to fit this tool to their preferred model of information decay. While not coded, the results of SeqDef can simulate what could be found using a Brownian motion model by using a very small value of lambda higher than zero."
- REPLACE: "The package also implements a parameter-free Brownian-motion kernel (the phylogenetic correlation under Brownian motion). This is exactly the flat, low-λ limit of the exponential (it matches the exponential at λ ≈ 0.5; ρ ≈ 0.98 on the Chondrichthyes tree) and, being parameter-free, cannot be calibrated to a phylogenetic horizon — on our data it is less discriminating and selects a different top target (Supplementary Fig. S5), which is precisely why we adopt the tunable exponential as the default."

**6b. Add an extension paragraph** (Practical Considerations):

> "As implemented, SeqDef prunes to the tree∩data intersection, so only placed taxa are scored. A natural
> extension — valuable for under-described tropical lineages — is to graft taxa lacking molecular data
> onto the phylogeny with taxonomic-insertion tools (e.g., rtrees, V.PhyloMaker, PASTIS/TACT) and then
> score the augmented tree; such taxa surface as high-deficiency targets. Soares et al. (2023) used
> exactly this insertion strategy for the fish tree of life, and our posterior-tree framework already
> integrates over uncertainty, so SeqDef could be averaged over stochastic insertions. A planned
> `prune = FALSE` option will make this a one-line user choice. Relatedly, a taxonomic-uncertainty weight
> could enter either as a trait passed to `calc_priority` or as a per-tip down-weight on S."

**6c. Conclusions — calibrate the efficiency claim.**
- FIND: "It is computationally efficient, interoperable with standard phylogenetic and bioinformatic R packages, and adaptable to context-dependent weighting schemes."
- REPLACE: "It is computationally efficient — for the default exponential kernel an exact linear-time, linear-memory tree traversal lets it scale to phylogenies of millions of tips, while the Gaussian and linear kernels use a dense O(n²) computation suited to clade-to-class scale. It is interoperable with standard phylogenetic and bioinformatic R packages and adaptable to context-dependent weighting schemes."

---

## 7. Figures

- **Fig. 1** — replace the image with the regenerated `figures/fig1.pdf` (λ-annotated) and update the caption:
  > "**Figure 1. Behaviour of the sequencing-deficiency statistic and the effect of the decay parameter λ.**
  > SeqDef on a 10-taxon toy tree in which only *taxa10* carries data (S = 1), shown for the linear kernel
  > and the exponential kernel at λ = 1.2 (auto_max) and λ = 5. As λ increases the phylogenetic horizon
  > shrinks: the sister of the sequenced tip (*taxa9*) is scored as well-covered only under a broad
  > horizon and becomes increasingly deficient as λ grows."
- **Supplementary figures** (files in `figures/`): add with captions —
  - **S1 `figS_lambda_stability.pdf`** — Priority-ranking stability vs λ (ρ and top-10 overlap vs the auto_max solution).
  - **S2 `figS_kernel.pdf`** — cross-kernel ranking agreement + MCC exponential-vs-Gaussian priority.
  - **S3 `figS_edge_vs_seqdef.pdf`** — SeqDef vs EDGE (near-orthogonal, ρ = 0.07); divergent species + *C. atromarginatus*.
  - **S4 `figS_runtime.pdf`** — runtime & memory vs n: dense O(n²) vs exponential traversal O(n) (to 10⁶ tips).
  - **S5 `figS_brownian.pdf`** — Brownian as the flat, low-λ limit of the exponential (variance-vs-λ; ranking scatter).

---

## 8. New references (add via Paperpile)

- Jukes, T.H. & Cantor, C.R. 1969. Evolution of protein molecules. In: *Mammalian Protein Metabolism* (H.N. Munro, ed.), pp. 21–132. Academic Press, New York.
- Gumbs, R., Gray, C.L., Böhm, M., Couchman, O.R., Faith, D.P., Forest, F., *et al.* 2023. The EDGE2 protocol: Advancing the prioritisation of Evolutionarily Distinct and Globally Endangered species. *PLoS Biology* 21: e3001991.
- Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D., *et al.* 2010. Picante: R tools for integrating phylogenies and ecology. *Bioinformatics* 26: 1463–1464.

*(Soares et al. 2023, Biological Conservation 285: 110235 — already cited.)*

---

## 9. Compliance (separate from the manuscript body)

- **Data Availability Statement:** Zenodo snapshot of the code (DOI: [mint]), the NCBI assembly table, and
  the analysis scripts; posterior trees cited to VertLife / Stein et al. (2018); IUCN data not
  redistributable — give Red List version + query + access date.
- **CRediT** author contributions [to complete].
- Outside the manuscript: **rotate the NCBI + IUCN keys** (still in public git history).
