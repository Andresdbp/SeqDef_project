# Manuscript change-list — JEB-2026-00140 (SeqDef)

Apply these to the manuscript on Google Drive. Edits are ordered **top-to-bottom through the
paper** for easy application; each is tagged with the reviewer point it answers. Format:

- **FIND** = quote the existing text (search for it in the doc).
- **REPLACE WITH** / **INSERT** = the new text.
- Numbers come from `results/SUMMARY.md` (reproducible via `analyses/run_all.R`).

Three things to do globally first:
1. **Set every occurrence of "SeqDef", "Priority", and label subscripts (norm, median, depth) in
   upright (roman) type**, reserving italics for true variables (*S*, *w*, *d*, *A*, *λ*, *n*, *i*, *f*). [R1 #5]
2. **Add the new references** (block at the end of this document) to the reference list.
3. Decide on the **n = 850 vs 877** discrepancy (see Methods edit M-7).

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

**Edit M-4 [R1 #2b].** Add a principled justification of the kernel and mention the new argument.
Append to the paragraph that begins "The choice of the weight matrix W=[wfi] is flexible…":

- **INSERT (after that paragraph):** "Although the framework admits any kernel, the package implements three (`kernel = "exponential"`, `"gaussian"`, `"linear"`) and we adopt the exponential kernel as the default. It is the natural choice because it yields a single, interpretable tunable horizon — the phylogenetic half-life — and is standard in distance-decay biodiversity measures (Pavoine et al., 2005) and sequence-novelty work (Marini et al., 2022). A kernel-comparison analysis (below; Supplementary Fig. S3) confirms that the prioritization is robust to this choice."

---

## METHODS — R implementation

**Edit M-5 [R1 #4].** Replace the unsupported efficiency claim with a complexity statement.

- **FIND:** "The deficiency calculation is fully vectorized using matrix algebra (WS), ensuring high computational efficiency even for large phylogenies."
- **REPLACE WITH:** "The deficiency calculation is fully vectorized as a matrix–vector product (**W S**). For *n* tips the cophenetic matrix, the weight matrix and the product are each O(*n*²) in time and memory, and automatic λ selection adds a constant factor (a scan of up to 491 λ values with early stopping), giving O(*k·n*²) with *k* ≤ 491. The method is therefore efficient at clade-to-class scale and remains tractable to ~10⁵ tips on large-memory hardware, with the dense O(*n*²) matrix becoming the binding constraint only near tree-of-life scale; we quantify this with empirical benchmarks below."

**Edit M-6 [R1 #2b].** Note the kernel argument in the function description.

- **FIND:** "the function includes two algorithms for automated lambda determination."
- **REPLACE WITH:** "the function includes two algorithms for automated lambda determination and a `kernel` argument selecting the distance-decay kernel (exponential, Gaussian, or linear; exponential by default)."

---

## METHODS — Practical Workflow (Chondrichthyes)

**Edit M-7 [admin — dataset count].** Reconcile n. The reproducible pipeline yields **877** species
(phylogeny ∩ IUCN; 55 with an NCBI assembly).

- **FIND:** "the intersection of species present in both the IUCN and phylogenetic datasets (n = 850)."
- **REPLACE WITH:** either "(n = 877)" — matching the current pipeline — **or** add the filter that
  produces 850 (e.g., excluding additional IUCN categories) and state it. *Verify against your IUCN
  download version before choosing.*

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
> n = 150,000. By extrapolation, whole-eukaryote scale (n ~ 10⁶) would require tens of terabytes,
> motivating sparse-kernel approximations for the largest trees. The 877-tip Chondrichthyes analysis
> ran in 0.02 s (fixed λ) / 0.21 s (auto_max) per tree, and the full 100-tree posterior completed in
> ≈40 s (Supplementary Fig. S6)."

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
- **REPLACE WITH:** "Additionally, the framework can accommodate the user's preferred model of information decay. As a qualitative analogy, a small λ produces a broad, slowly decaying weighting reminiscent of a Brownian-motion covariance structure; a formal Brownian-covariance kernel could be added but is not implemented here, and we do not claim the exponential kernel reproduces it exactly."

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
- **REPLACE WITH:** "It is computationally efficient at clade scale (hundreds to a few thousand tips), with an explicit O(n²)-memory cost that makes naive whole-tree-of-life application memory-bound; because distant pairs contribute negligibly, a truncated or nearest-neighbour kernel is a natural route to larger trees. It is interoperable with standard phylogenetic and bioinformatic R packages and adaptable to context-dependent weighting schemes."

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
  calculation versus tree size n (log–log; single-threaded; n = 100 to 100,000 on a large-memory
  node). Points are medians over up to 100 replicate trees; dotted lines are O(n²) references.
  Fitted slopes are 2.09 (time) and 1.94 (memory), confirming quadratic scaling; peak memory reaches
  447 GB at n = 100,000."

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
