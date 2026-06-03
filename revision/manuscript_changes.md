# Manuscript change-list — JEB-2026-00140 (SeqDef)

Apply these to the manuscript on Google Drive. Edits are ordered **top-to-bottom through the
paper** for easy application; each is tagged with the reviewer point it answers. Format:

- **FIND** = quote the existing text (search for it in the doc).
- **REPLACE WITH** / **INSERT** = the new text.
- Numbers come from `results/SUMMARY.md` (reproducible via `analyses/run_all.R`).

> **Reconciled against the live Google Doc ("manuscript", 2026-06-03).** Every FIND anchor below
> matches the document's *visible* text; the equations are present as the same plain text shown
> (e.g. `wfi=exp(-dfi Tdepth)`, `Priority = SeqDefmedian2GE`), so the find/replace works directly.
> Two formatting notes: (a) citations are **Paperpile-managed links** (italic *et al.* + hyperlink) —
> when adding the new references (Soares 2023, Gumbs 2023, Kembel 2010) insert them **through
> Paperpile** so they join the managed reference list, or add them manually and append the entries
> from the "New references" block; (b) the manuscript has lost several standalone **λ** glyphs — see
> the consolidated Edit M-3b.

Three things to do globally first:
1. **Set every occurrence of "SeqDef", "Priority", and label subscripts (norm, median, depth) in
   upright (roman) type**, reserving italics for true variables (*S*, *w*, *d*, *A*, *λ*, *n*, *i*, *f*). [R1 #5]
2. **Add the new references** (block at the end of this document) to the reference list.
3. **Replace "n = 850" with "n = 877"** (data-derived; see Methods edit M-7).

---

## ABSTRACT

**Edit A-1 [R1 #1].** The abstract implies a single binary notion of data. Add one clause so the
generality of S is visible up front.

- **FIND:** "we use an exponential kernel with a tunable decay parameter that controls the phylogenetic scale of influence."
- **REPLACE WITH:** "we use an exponential kernel with a tunable decay parameter that controls the phylogenetic scale of influence. The availability input is a general score on the unit interval — from a simple binary presence/absence to a continuous genome-quality measure — supplied by the user."

---

## INTRODUCTION

**Edit I-1 [R2 #2.1, m6].** Reframe the two traditions as complementary, not opposed.

- **FIND:** "While these approaches demonstrate that combining phylogeny and extinction threat level can produce targeted conservation priorities, simply overlapping phylogeny and threat level ignores the reality that we have already sequenced many lineages. Large sequencing efforts require a related but distinct perspective that explicitly accounts for the current distribution of sequence resources across a phylogeny so that sampling decisions favor taxa that contribute novel genomic information, which can then be flexibly weighted by any user-defined trait of interest (e.g., vector or pest status, extinction risk, economic impact)."
- **REPLACE WITH:** "These distinctiveness-and-threat frameworks are complemented by a second tradition — information-maximizing sampling schemes — that selects targets according to how much *new* information they add relative to what is already sampled (Marini et al., 2022). The two traditions are complementary rather than competing: the first asks which lineages are most valuable and imperilled, the second asks where sampling effort is most informative given the data already in hand. Large genome-sequencing efforts sit at their intersection, because the value of sequencing a taxon depends on the genomes already available in its relatives. SeqDef bridges the two: it imports the 'condition on what is already sampled' logic of sampling schemes into the conservation-prioritization tradition, producing a genomic-deficiency score that can then be weighted by any user-defined trait (e.g., extinction risk, vector or pest status, economic impact). Notably, the recent EDGE2 framework already conditions a species' score on the extinction risk of its relatives (Gumbs et al., 2023), so the field is converging on exactly the conditioning logic that SeqDef formalizes for genomic-data availability."

**Edit I-2 [R2 #2.5, m9 line 40].** Add the Soares et al. (2023) citation to the existing-frameworks sentence.

- **FIND:** "with extinction threat status (Faith, 1992; e.g., EDGE, Isaac et al., 2007)."
- **REPLACE WITH:** "with extinction threat status (Faith, 1992; e.g., EDGE, Isaac et al., 2007; EDGE2, Gumbs et al., 2023; Soares et al., 2023)."

**Edit I-3 [R2 #2.5, m7 line 48].** Define "marginal information value" at first use.

- **FIND:** "SeqDef quantifies (for any focal taxon) the degree to which existing sequence data in related taxa reduces the marginal information value of sequencing that focal taxon."
- **REPLACE WITH:** "SeqDef quantifies, for any focal taxon, its *marginal information value* — by which we mean the expected gain in genome-derived comparative information from adding that taxon's genome to the current set, given the genomes already available in its relatives. A taxon has high marginal value when little of its genomic information is already captured by sequenced close relatives; operationally this is the complement of phylogenetically weighted availability (1 − A_f, defined below)."

---

## METHODS — Formal definition

**Edit M-1 [R1 #1, R2 #2.3].** Generalize S and add worked continuous constructions.

- **FIND:** "define a scalar Si[0,1] quantifying the availability of sequence data (where Si=1 indicates high availability, e.g., a completed genome, and Si=0 indicates absence)."
- **REPLACE WITH:** "define a scalar *S*_i ∈ [0,1] quantifying the availability of sequence data for tip *i*. *S* is a general, user-supplied score: any monotone construction is admissible. Examples include (i) a binary indicator of whether any assembly exists (*S* ∈ {0,1}); (ii) genome completeness, e.g., the BUSCO complete fraction; (iii) a composite of rescaled BUSCO completeness, contiguity (log₁₀ N50) and base accuracy (QV); (iv) data-type breadth, the fraction of desired resources present (e.g., {genome, transcriptome, resequencing} → 0, ⅓, ⅔, 1); or (v) a saturating function of assembly count such as *S* = *n*/(*n*+1). The statistic requires only that *S*_i increases with data availability; the case study below uses the binary form, and a supplementary analysis demonstrates the continuous case."

**Edit M-2 [R1 #5 — the missing-λ bug].** Fix the kernel equation and its lead-in sentence.

- **FIND (sentence):** "The weights wfi are derived from a distance-decay kernel based on the cophenetic distance dji between tips, normalized by tree depth Tdepth, and scaled by , a user-specified decay parameter that controls the phylogenetic scale of influence:"
- **REPLACE WITH:** "The weights *w*_fi are derived from a distance-decay kernel based on the cophenetic distance *d*_fi between tips, normalized by tree depth *T*_depth and scaled by *λ*, a user-specified decay parameter that controls the phylogenetic scale of influence:"
- **FIND (displayed equation, currently missing λ):** `w_fi = exp(− d_fi / T_depth)`  *(rendered in the manuscript as "wfi=exp(-dfi Tdepth)")*
- **REPLACE WITH (Eq. 2):**  **w_fi = exp( − λ · d_fi / T_depth )**
  *(this is the equation the text and the code both use; the λ and the division had been dropped).*

**Edit M-3 [R1 #5].** Ensure the availability and SeqDef equations render with roman operators / proper subscripts.

- **A_f equation (Eq. 1):** render as **A_f = Σ_{i=1}^{n} w_fi · S_i** (sum index explicit).
- **SeqDef equation (Eq. 3):** render as **SeqDef_f = 1 − A_{f,norm}**, with "SeqDef" upright and "norm" an upright subscript.

**Edit M-3b [R1 #5 — dropped λ glyphs].** The live document has lost the standalone λ symbol in several
places (it renders as an empty "()" or a blank). Verified against the Drive doc — restore λ at each:
- Methods: **FIND** "scaled by , a user-specified decay parameter" → **"scaled by λ, a user-specified decay parameter"** (same as M-2).
- Methods (kernel choice): **FIND** "the tunable decay parameter lambda ()," → **"the tunable decay parameter λ,"**
- Discussion (Practical Considerations): **FIND** "the choice of the decay parameter () dictates" → **"the choice of the decay parameter λ dictates"**
- Discussion (Practical Considerations): **FIND** "anchor  to the genus-level phylogenetic half-life" → **"anchor λ to the genus-level phylogenetic half-life"**
  *(Tip: in Google Docs, searching the blank "()" is hard — search the surrounding words shown above.)*

**Edit M-4 [R1 #2b].** Add a principled justification of the kernel and mention the new argument.
Append to the paragraph that begins "The choice of the weight matrix W=[wfi] is flexible…":

- **INSERT (after that paragraph):** "Although the framework admits any kernel, the package implements four — `kernel = "exponential"` (default), `"gaussian"`, `"linear"`, and `"brownian"` — encoding different models of how informational redundancy relates to evolutionary divergence. The **exponential** kernel, exp(−λd), models a constant *proportional*, memoryless decay (redundancy halves every fixed distance — the phylogenetic half-life), mirroring how molecular similarity saturates with divergence under a clock-like process. The **Gaussian** kernel, exp(−λd²), is flat at small distances then falls sharply (a buffer of near-equivalence with a crisp outer boundary). The **linear** kernel, max(0, 1−λd), declines at a constant absolute rate to a hard cutoff (an explicit bounded horizon). The **Brownian** kernel weights by shared ancestry (the phylogenetic correlation under Brownian motion) and is *parameter-free*. We adopt the exponential as the default for four reasons: (i) it is the natural model of molecular-information decay; (ii) it provides a single interpretable, tunable horizon (the half-life); (iii) it is standard in distance-decay biodiversity measures (Pavoine et al., 2005) and sequence-novelty work (Marini et al., 2022); and (iv) because exp(−λd) factorizes into independent per-branch factors, it admits an *exact* O(n)-time, O(n)-memory algorithm (see Computational benchmarks). The exponential thus **uniquely combines O(n) scalability with a tunable horizon**: the Gaussian and linear kernels are O(n²), while the Brownian kernel — though also O(n) — is parameter-free and corresponds to the exponential at a fixed, small λ (≈0.5), i.e. its flat, low-discrimination limit (on the Chondrichthyes data it is less discriminating and selects a different top target; Supplementary Fig. S7). The same memoryless, per-branch-multiplicative property underlies both the exponential's biological naturalness and its computational tractability. A kernel-comparison analysis (Supplementary Fig. S3) confirms the prioritization is robust across the distance-decay kernels."

---

## METHODS — R implementation

**Edit M-5 [R1 #4].** Replace the unsupported efficiency claim with a complexity statement.

- **FIND:** "The deficiency calculation is fully vectorized using matrix algebra (WS), ensuring high computational efficiency even for large phylogenies."
- **REPLACE WITH:** "The deficiency calculation is fully vectorized as a matrix–vector product (**W S**). For *n* tips the cophenetic matrix, the weight matrix and the product are each O(*n*²) in time and memory, and automatic λ selection adds a constant factor (a scan of up to 491 λ values with early stopping), giving O(*k·n*²) with *k* ≤ 491. **For the default exponential kernel this cost is avoided entirely:** because exp(−λ*d*) factorizes along tree paths, the weighted availability is computed by an exact two-pass tree traversal in O(*n*) time and memory, forming no distance matrix at all; the O(*n*²) cost therefore applies only to the Gaussian and linear kernels. We quantify both regimes with empirical benchmarks below."

**Edit M-6 [R1 #2b].** Note the kernel argument in the function description.

- **FIND:** "the function includes two algorithms for automated lambda determination."
- **REPLACE WITH:** "the function includes two algorithms for automated lambda determination and a `kernel` argument selecting the weighting kernel — exponential (default), Gaussian, linear, or a parameter-free Brownian-motion correlation."

---

## METHODS — Practical Workflow (Chondrichthyes)

**Edit M-7 [admin — dataset count].** Reconcile n. The reproducible pipeline yields **877** species
(phylogeny ∩ IUCN; 55 with an NCBI assembly).

- **FIND:** "the intersection of species present in both the IUCN and phylogenetic datasets (n = 850)."
- **REPLACE WITH:** "the intersection of species present in both the IUCN and phylogenetic datasets (n = 877)."
  (877 is the data-derived count reproduced by the pipeline; the earlier 850 was an error.)

**Edit M-8 [R1 #1, R2 #2.3].** Justify the binary coding and point to the continuous analysis.

- **FIND:** "We defined availability (s) as binary (1 = presence of any assembly in NCBI GenBank; 0 = absence)."
- **REPLACE WITH:** "We defined availability (*S*) as binary (1 = presence of any assembly in NCBI GenBank; 0 = absence). For this demonstration the binary coding is deliberate: assembly presence/absence is objective, uniformly obtainable for all species from a single NCBI query, avoids confounding the demonstration with heterogeneous and frequently missing quality metrics, and matches the Earth BioGenome framing of whether a reference genome yet exists. To show that the statistic behaves sensibly with a continuous input, we repeated the analysis with an illustrative continuous score derived from assembly count, *S* = count/(count+1) (Supplementary Fig. S5); the binary and continuous prioritizations are highly concordant (Spearman ρ = 0.98 for Priority; identical top-10), while scores redistribute toward species whose relatives have only minimal coverage."

**Edit M-9 [R1 #5].** Fix the Priority equation typography (Eq. 4).

- **FIND (displayed):** `Priority = SeqDef_median 2^GE`  *(rendered as "Priority = SeqDefmedian2GE")*
- **REPLACE WITH (Eq. 4):**  **Priority = SeqDef_median × 2^GE**
  with "Priority" and "SeqDef" upright (roman) and "median" an upright subscript; *GE* italic.

**Edit M-10 [R2 #2.5, m9 lines 178–180].** Rewrite the multiplicative-effect sentence.

- **FIND:** "This formula creates a multiplicative effect, assigning the highest priority to taxa that are both phylogenetically isolated from existing genomes and critically threatened."
- **REPLACE WITH:** "Because GE enters as an exponent of 2, each one-category increase in extinction risk **doubles** a taxon's priority: a Critically Endangered taxon (GE = 4) receives 2⁴ = 16× the weight of a Least Concern taxon (GE = 0) with the same SeqDef. Priority therefore rises steeply for taxa that are simultaneously isolated from existing genomes and highly threatened."

---

## METHODS — Sensitivity analysis (expand this subsection) [R1 #2a, #2b]

**Edit M-11.** After the existing variance-vs-λ description, **INSERT** the following two paragraphs
(new analyses; numbers from Analysis 1 & 2):

> "Beyond the variance of the statistic, we asked whether λ changes the *prioritization itself*.
> Across λ ∈ [1, 25] we computed the Spearman correlation of the SeqDef and Priority rankings
> against each tree's automatically selected (auto_max) solution. The Priority ranking was highly
> stable (ρ = 0.93 median, 0.88 minimum across the interval), and *Centrophorus atromarginatus* was
> the single highest-priority taxon in 65–80 % of sampled trees at every λ, peaking near the
> auto-selected band (λ ≈ 2.6–5; Supplementary Fig. S2). We further compared the two automatic
> λ-selection algorithms directly: across the 100 posterior trees they selected similar values
> (variance-maximization median λ = 3.40 [2.65, 5.00]; genus-scale median λ = 3.19 [2.87, 4.87]),
> their Priority outputs correlated at ρ = 0.998, and they agreed on the single top target in
> 100 of 100 trees. We therefore recommend variance-maximization as the default for discriminatory
> power and genus-scale calibration where a biologically interpretable horizon is preferred.
>
> To assess the kernel choice, we recomputed SeqDef and Priority under exponential, Gaussian and
> linear kernels (each with its own auto_max λ). Priority rankings were highly correlated with the
> exponential default (Gaussian ρ = 0.88; linear ρ = 0.94) and the top target was preserved in the
> large majority of trees (Supplementary Fig. S3), indicating that the prioritization is governed by
> tree topology and the data distribution rather than by the specific kernel form."

**Edit M-12 [R1 #3].** Add the EDGE benchmark to Methods (brief) — **INSERT** a short paragraph:

> "To position SeqDef relative to established conservation metrics, we benchmarked it against EDGE
> (Isaac et al., 2007) and EDGE2 (Gumbs et al., 2023) on the Chondrichthyes MCC tree. We computed
> Evolutionary Distinctiveness (fair-proportion; picante, Kembel et al., 2010), formed classic EDGE
> = log(1 + ED) + GE·log 2, and computed EDGE2 as expected phylogenetic-diversity loss using the
> standard IUCN→50-year extinction-probability mapping (Gumbs et al., 2023). We compared each with
> SeqDef-Priority by Spearman correlation, top-25 list overlap, and the species at which the metrics
> most diverge."

---

## METHODS — Runtime benchmark (new short subsection) [R1 #4]

**Edit M-13.** **INSERT** a brief subsection (numbers from the Grace benchmark):

> "**Computational benchmarks.** We timed a single SeqDef calculation (fixed λ, single-threaded BLAS)
> on simulated trees (`ape::rtree`) with random binary S at 36 sizes from n = 100 to 100,000, on a
> large-memory compute node (80 cores, 3 TB RAM; R 4.4.2). Both runtime and peak memory scaled as
> O(n²) (fitted log–log slopes of 2.09 and 1.94, respectively, for n ≥ 1,000): a single call took
> 0.04 s at n = 1,000, 6.2 s at n = 10,000 and 808 s (≈13.5 min) at n = 100,000. Peak memory grew
> from ~4.7 GB at n = 10,000 to 447 GB at n = 100,000; the dense O(n²) distance and weight matrices —
> not runtime — are the binding constraint, and a benchmark run exceeded an 800 GB allocation at
> n = 150,000 (these O(n²) figures govern the Gaussian and linear kernels). **For the default
> exponential kernel, the exact O(n) traversal removes this constraint entirely:** it returns
> identical scores (agreement ≤ 1e-16) but a single call runs at n = 100,000 in 0.06 s using 74 MB
> (versus 808 s and 447 GB for the dense computation) and at n = 1,000,000 in under 1 s using ~0.4 GB,
> on a laptop — so prioritization scales exactly to whole-tree-of-life size for the default kernel.
> The 877-tip Chondrichthyes analysis ran in 0.02 s (fixed λ) / 0.21 s (auto_max) per tree, and the
> full 100-tree posterior completed in ≈40 s (Supplementary Fig. S6)."

---

## RESULTS — add EDGE benchmark subsection [R1 #3]

**Edit R-1.** **INSERT** a new subsection after "Chondrichthyes Sequencing Recommendation":

> "**SeqDef is complementary to EDGE.** On the Chondrichthyes MCC tree, SeqDef and EDGE were
> near-orthogonal (Spearman ρ = 0.07): they quantify different things — evolutionary distinctiveness
> plus threat versus genomic-data deficiency. SeqDef-Priority correlated with EDGE (ρ = 0.83) and
> EDGE2 (ρ = 0.79) only because all three share the IUCN threat term, and the top-25 lists overlapped
> by 9 and 8 species, respectively. The informative cases are where the metrics diverge. The
> scalloped hammerhead *Sphyrna lewini* (Critically Endangered) ranks high on EDGE (rank 40) but low
> on SeqDef (rank 779 of 877), because its congener *S. mokarran* already has a genome assembly, so
> its marginal genomic information is largely captured. In contrast, *C. atromarginatus* ranks high
> on both because the genus *Centrophorus* and the family Centrophoridae contain no assembly at all.
> SeqDef thus identifies a distinct, dynamic axis of priority — marginal genomic information given the
> data already in hand — that updates as relatives are sequenced and that static distinctiveness
> metrics do not capture (Supplementary Fig. S4)."

---

## DISCUSSION

**Edit D-1 [R1 #2b].** Soften the Brownian-motion claim.

- **FIND:** "Additionally, we want to allow the users the ability to fit this tool to their preferred model of information decay. While not coded, the results of SeqDef can simulate what could be found using a Brownian motion model by using a very small value of lambda higher than zero."
- **REPLACE WITH:** "Additionally, the framework accommodates the user's preferred model of information decay; the package implements a parameter-free Brownian-motion kernel (the phylogenetic correlation under Brownian motion). This kernel is exactly the flat, low-λ limit of the exponential (it matches the exponential at λ ≈ 0.5; Spearman ρ ≈ 0.98 on the Chondrichthyes tree), and being parameter-free it cannot be calibrated to a phylogenetic horizon — on our data it is less discriminating and selects a different top target (Supplementary Fig. S7). This makes concrete why we adopt the tunable exponential as the default."

**Edit D-2 [R2 #2.2, R2 #2.4].** Add an extension paragraph (rtrees + taxonomic uncertainty).

- **INSERT** (new Discussion paragraph, e.g., under "Practical Considerations"):

> "**Extending SeqDef to unsampled taxa.** As implemented, SeqDef prunes to the intersection of the
> tree and the data, so only taxa already placed on the phylogeny are scored. A natural extension —
> particularly valuable for under-described tropical lineages — is to first graft taxa that lack any
> molecular data onto the phylogeny using taxonomic-insertion tools (e.g., rtrees, V.PhyloMaker,
> `phytools::add.species.to.genus`, or PASTIS/TACT) and then score the augmented tree; such taxa,
> having no data and often sitting in undersampled clades, will surface as high-deficiency targets.
> Soares et al. (2023) used exactly this insertion-based strategy to build a synthesis phylogeny and
> prioritization for the fish tree of life, and our posterior-tree framework already integrates over
> phylogenetic uncertainty, so SeqDef could be averaged over many stochastic insertions to propagate
> placement error. We plan a `prune = FALSE` option so users can supply an augmented tree directly.
> Relatedly, where taxonomic-placement confidence is known, it can be incorporated either as an extra
> trait passed to `calc_priority` or as a per-tip down-weight on S, allowing users to discount
> confidently for uncertain tips; we leave this as a lightweight future option rather than a change to
> the core statistic."

**Edit D-3 [R1 #4].** In "Practical Considerations" / Conclusions, calibrate the efficiency wording.

- **FIND (Conclusions):** "It is computationally efficient, interoperable with standard phylogenetic and bioinformatic R packages, and adaptable to context-dependent weighting schemes."
- **REPLACE WITH:** "It is computationally efficient: for the default exponential kernel an exact linear-time, linear-memory tree traversal lets it scale to phylogenies of millions of tips, while the Gaussian and linear kernels use a dense O(n²) computation suited to clade-to-class scale. It is interoperable with standard phylogenetic and bioinformatic R packages and adaptable to context-dependent weighting schemes."

---

## FIGURES

**Fig. 1 (replace).** Use the regenerated `figures/fig1.pdf` (built from a real `SeqDef()` run on the
toy tree). New caption:

> "**Figure 1. Behaviour of the sequencing-deficiency statistic and the effect of the decay parameter
> λ.** SeqDef scores for a 10-taxon toy phylogeny in which only *taxa10* carries sequence data
> (availability S = 1). Columns show SeqDef computed at three decay parameters — λ = 1.2 (selected by
> the variance-maximization algorithm), 5 and 15. As λ increases the 'phylogenetic horizon' shrinks:
> the sister of the sequenced tip (*taxa9*) is scored as well-covered (low SeqDef) only under a broad
> horizon (0.23 at λ = 1.2) and becomes increasingly deficient as λ grows (0.95 at λ = 15), while
> distant tips are deficient throughout."

**Fig. 2 (caption update) [R1 #2a].** Add a sentence reporting the head-to-head result:

> "...the secondary genus-scale algorithm. Across the 100 posterior trees the two algorithms selected
> similar λ (medians 3.40 and 3.19) and produced almost identical prioritizations (Spearman ρ = 0.998),
> agreeing on the single highest-priority taxon in all 100 trees."

**New supplementary figures** (files in `figures/`): add with these captions.

- **Fig. S2 — `figS_lambda_stability.pdf` [R1 #2a].** "Stability of the prioritization to the decay
  parameter λ. (A) Spearman correlation of the Priority ranking with the auto_max solution across
  λ ∈ [1,25] (mean ± 95 % band over sampled posterior trees); the red band marks the 95 % interval of
  λ selected by auto_max. (B) Top-10 overlap with the auto_max ranking."
- **Fig. S3 — `figS_kernel.pdf` [R1 #2b].** "Robustness of the prioritization to the kernel. (A)
  Cross-kernel Spearman correlation of Priority rankings across posterior trees. (B) Per-species
  Priority on the MCC tree under the exponential vs Gaussian kernels; *C. atromarginatus* in red."
- **Fig. S4 — `figS_edge_vs_seqdef.pdf` [R1 #3].** "SeqDef versus EDGE on the Chondrichthyes MCC tree
  (near-orthogonal, ρ = 0.07). Purple: species that are high-EDGE but low-SeqDef because a congener is
  already sequenced (e.g., *Sphyrna lewini*); red: *C. atromarginatus*, high on both."
- **Fig. S5 — `figS_continuousS.pdf` [R1 #1].** "Binary versus continuous availability S. SeqDef (A)
  and Priority (B) computed with binary S vs an illustrative continuous score S = count/(count+1);
  rankings are highly concordant (ρ = 0.96 and 0.98)."
- **Fig. S6 — `figS_runtime.pdf` [R1 #4].** "Runtime (A) and peak memory (B) of a single SeqDef
  calculation versus tree size n (log–log; single-threaded). Orange: the dense O(n²) computation used
  by the Gaussian/linear kernels (n = 100 to 100,000; peak memory 447 GB at n = 100,000, beyond which
  it exhausts an 800 GB allocation). Green: the exact O(n) tree traversal used by the default
  exponential kernel (n = 100 to 1,000,000; ~0.4 GB and < 1 s at n = 1,000,000, on a laptop). Dotted
  lines are O(n²) and O(n) references. The two methods give identical scores (≤ 1e-16); the
  exponential kernel scales exactly to whole-tree-of-life size."
- **Fig. S7 — `figS_brownian.pdf` [R1 #2b].** "The Brownian-motion kernel as the flat limit of the
  exponential. (A) Variance of SeqDef scores on the Chondrichthyes MCC tree as a function of λ; the
  variance-maximizing `auto_max` value (λ ≈ 4.8) is far from the Brownian limit (λ ≈ 0.5), which sits
  in the low-discrimination region. (B) Per-species SeqDef under the exponential (`auto_max`) versus
  the parameter-free Brownian kernel; the two rank taxa differently (Spearman ρ ≈ 0.6), and Brownian
  selects a different top target — illustrating the value of a tunable horizon."

---

## DATA AVAILABILITY STATEMENT (new section) [admin]

> "**Data availability.** The SeqDef R package and all analysis scripts are archived at
> [Zenodo DOI — to be minted] and developed at https://github.com/Andresdbp/SeqDef. The NCBI assembly
> presence/absence table and the analysis outputs are included in the archive. The posterior
> phylogenies are from VertLife (Stein et al., 2018; [URL/DOI]). IUCN Red List assessments are not
> redistributable under IUCN's terms; we provide the Red List version, the exact query, and the access
> date so the assessment data can be reproduced [version __; accessed __]."

**CRediT (new section).** Add author contributions in CRediT taxonomy [to be completed by authors —
required at resubmission and not editable after acceptance].

---

## NEW REFERENCES TO ADD

- Gumbs, R., Gray, C.L., Böhm, M., Burfield, I.J., Couchman, O.R., Faith, D.P., Forest, F., Hoffmann, M., Isaac, N.J.B., Jetz, W., Mace, G.M., Mooers, A.O., Safi, K., Scott, O., Steel, M., Tucker, C.M., Pearse, W.D., Owen, N.R. & Rosindell, J. 2023. The EDGE2 protocol: Advancing the prioritisation of Evolutionarily Distinct and Globally Endangered species for practical conservation action. *PLoS Biology* 21: e3001991.
- Soares, B.E., Nakamura, G., Freitas, T.M.S., Richter, A. & Cadotte, M.W. 2023. Quantifying and overcoming Darwinian shortfalls to conserve the fish tree of life. *Biological Conservation* 285: 110235.
- Kembel, S.W., Cowan, P.D., Helmus, M.R., Cornwell, W.K., Morlon, H., Ackerly, D.D., Blomberg, S.P. & Webb, C.O. 2010. Picante: R tools for integrating phylogenies and ecology. *Bioinformatics* 26: 1463–1464.

*(The 50-year IUCN→extinction-probability values used for EDGE2 are those given by Gumbs et al., 2023 — cite that paper for the mapping.)*
