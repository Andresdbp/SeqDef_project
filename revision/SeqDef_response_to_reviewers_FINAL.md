# Response to Reviewers

**SeqDef: An R package for phylogenetically weighted genomic prioritization across the tree of life**
Manuscript JEB-2026-00140 · Journal of Evolutionary Biology

*Reviewer and editor comments are in italics. Our replies follow under "Response." All
numbers below are from the new analyses scripted in `analyses/01–06` and summarised in
`results/SUMMARY.md` (Apple Silicon, 8 cores, 16 GB RAM, R 4.5.2; `set.seed(42)`).
Manuscript line numbers will be inserted against the tracked-changes file at submission;
each edit is specified in the accompanying manuscript change-list.*

---

## Response to the Associate Editor

We thank the Associate Editor and both reviewers for a constructive assessment. We have
addressed every point. The five major concerns required clarification and evidence rather
than any change to the statistic itself; we ran four new analyses (all on data already in
the repository), regenerated Figure 1, added a `kernel`
argument to the package, and corrected the notation.

**Summary of major changes**

- **Variable S.** Clarified that S is a general availability score on [0,1]; the binary
  case-study coding is framed as a deliberate operational choice, with worked continuous
  constructions given in the Methods (no separate continuous-S run is needed).
- **λ and kernel.** Added a λ rank-stability analysis, a head-to-head of the two λ-selection
  methods (ρ = 0.998; identical top target in 100/100 trees), and a kernel comparison via a
  new `kernel` argument (exponential/Gaussian/linear; rankings ρ ≥ 0.88).
- **Benchmark.** Added a direct EDGE / EDGE2 comparison framed around complementarity
  (EDGE vs SeqDef ρ = 0.07 — they are near-orthogonal).
- **Efficiency.** Added an O(n²) complexity statement and empirical runtime/memory benchmarks
  with hardware specs; recalibrated the efficiency claims.
- **Notation.** Corrected Eq. 4 italics, **restored the missing λ in the kernel equation**,
  and standardised notation.
- **Framing/definitions.** Reframed the Introduction as complementary traditions, defined
  "marginal information value," added an insertion-based extension paragraph, and added the
  suggested citations.

---

## Reviewer 1

### 1.1 Inconsistency in the definition of Sᵢ
> *In line 62 "Si" is continuous on [0,1]; in line 164 it is treated as binary … If continuous,
> explain how it is computed (e.g., BUSCO); if binary, justify the simplification.*

**Response.** This was a presentation issue, not a methodological one. S is a general
continuous availability score on [0,1]: the function multiplies the phylogenetic weight
matrix by the numeric S vector (`dist.prop %*% as.numeric(df[, data.col])`) and never
requires binary input. We have separated the general definition from the case-study choice:

- In the **Methods** we now state S ∈ [0,1] is any user-supplied, monotone availability
  measure, and give worked continuous constructions: (i) genome completeness = BUSCO complete
  fraction; (ii) a composite of rescaled BUSCO completeness, contiguity (log₁₀ N50) and base
  accuracy (QV); (iii) data-type breadth (fraction of {genome, transcriptome, resequencing});
  and (iv) a saturating function of assembly count, S = n/(n+1).
- Where the **Chondrichthyes pipeline** is introduced we justify the binary coding: assembly
  presence/absence is objective, uniformly obtainable for all species from a single NCBI query,
  avoids confounding the demonstration with heterogeneous and often-missing quality metrics, and
  matches the Earth BioGenome framing of whether a reference yet exists.
- The binary case-study coding is retained and justified explicitly; the worked constructions
  above show how a continuous S is specified, so we did not add a separate continuous-S run.

### 1.2 Dependence on λ and the choice of kernel
> *Results appear highly dependent on λ and the exponential kernel … the two λ-selection methods
> are not compared, and sensitivity to λ is not assessed … the kernel choice should be justified
> by systematic comparison.*

**Response.** We added three analyses.

**(a) λ rank-stability.** Across λ ∈ [1, 25] we computed the Spearman correlation of the
SeqDef and Priority rankings against each tree's own `auto_max` solution. The **Priority ranking
is highly stable** (ρ median = **0.932**, minimum **0.878** across the whole interval; SeqDef
ρ median 0.818). *C. atromarginatus* is the top target in **65–80 % of sampled trees at every λ
in [1,25]** (never below 50 %), peaking at ~80 % within the band both automatic methods select
(λ ≈ 2.6–5). We therefore present this as the λ analogue of our topological-robustness test
(Fig. 3): the prioritization is robust across a broad λ band, and the top target is the modal,
auto-selected choice rather than an artifact of one λ (Supplementary Fig. S1).

**(b) The two λ-selection methods are essentially interchangeable.** Across the 100 posterior
trees, `auto_max` selects λ median **3.40** (95 % interval [2.65, 5.00]) and `by_genus` selects
**3.19** ([2.87, 4.87]); the two Priority vectors correlate at ρ = **0.998** ([0.985, 1.000]) and
**agree on the single top target in 100/100 trees** (each selects *C. atromarginatus* in 85/100).
We recommend `auto_max` as the default for discriminatory power and `by_genus` where a
biologically interpretable horizon is preferred, and we state this (revised Fig. 2).

**(c) Kernel comparison.** We added a `kernel` argument (`exponential`/`gaussian`/`linear`),
applied in both the final calculation and the `auto_max` optimisation. On the Chondrichthyes data
(each kernel with its own `auto_max` λ) the Priority rankings are highly correlated with the
exponential default (Gaussian ρ = **0.882** [0.846, 0.961]; linear ρ = **0.943** [0.923, 0.958]),
and the top target is preserved (in the top set for 16/20, 15/20 and 16/20 posterior trees,
respectively). We also justify the exponential kernel on principle — it yields a single,
interpretable "phylogenetic half-life" horizon and is standard in distance-decay biodiversity
measures (Pavoine et al. 2005) and sequence-novelty work (Marini et al. 2022) — and we
implemented a fourth, parameter-free **Brownian-motion kernel** (phylogenetic correlation) and use it
to make the earlier Brownian remark precise: it is exactly the flat, low-λ limit of the exponential
(Spearman ρ ≈ 0.98 with the exponential at λ ≈ 0.5), and on the Chondrichthyes data it is less
discriminating and selects a different top target — illustrating the value of the exponential's
tunable horizon (Supplementary Fig. S5). We additionally note a decisive computational argument for
the exponential default (developed under 1.4): among the distance-decay kernels, exp(−λd) is the only
one that factorizes along tree paths, so it alone admits an exact O(n)-time, O(n)-memory algorithm,
whereas the Gaussian (d²) and linear (1−λd) kernels are O(n²). (The Brownian kernel is also O(n) via
the three-point structure, but parameter-free.) The exponential therefore **uniquely combines O(n)
scalability with a tunable horizon**, and the same memoryless, per-branch-multiplicative property
underlies both its biological naturalness and its computational tractability.

### 1.3 Lack of comparison with existing methods
> *The manuscript mentions EDGE but does not compare against it.*

**Response.** We added a direct benchmark on the Chondrichthyes MCC tree (new Results
subsection). We computed Evolutionary Distinctiveness (fair-proportion), classic EDGE
(Isaac et al. 2007) and EDGE2 (Gumbs et al. 2023; expected-PD-loss formulation with the
standard 50-yr IUCN→extinction-probability map), and compared them with SeqDef-Priority.

The central result is **complementarity, not competition**: EDGE and SeqDef are
**near-orthogonal (ρ = 0.07)** — they answer different questions. Priority correlates with EDGE
(ρ = 0.832) and EDGE2 (ρ = 0.790) only because both carry the shared IUCN threat term; the
top-25 lists overlap 9/25 and 8/25. The divergence is the value proposition, and we give a worked
example: the scalloped hammerhead *Sphyrna lewini* (Critically Endangered) ranks high on EDGE
(rank 40) but **low on SeqDef (rank 779/877)** because its congener *Sphyrna mokarran* already has
a genome assembly — its marginal genomic information is largely captured. By contrast,
*C. atromarginatus* ranks high on **both** because the entire genus *Centrophorus* (8 species) and
family Centrophoridae have **no assembly at all**. We note that EDGE2 already conditions a
species' score on its relatives' extinction risk, so SeqDef extends the same conditioning logic
to genomic-data availability, and — unlike static EDGE — is **dynamic**, updating as relatives
are sequenced (Supplementary Fig. S3; divergence table in the supplement).

### 1.4 Claims about computational efficiency
> *Efficiency is asserted twice with no evidence; provide complexity analysis and/or empirical
> benchmarks with hardware and input-size scaling.*

**Response.** We replaced the unsupported assertions with both.

- **Complexity.** For n tips the cophenetic distance matrix, the weight matrix and the
  matrix–vector product are each **O(n²)** in time and memory; `auto_max` scans up to 491 λ
  values with early stopping, giving **O(k·n²)** with k ≤ 491. The binding constraint is the
  dense O(n²) matrix in memory.
- **Empirical (single λ; TAMU Grace large-memory node, x86_64, 80 cores, 3 TB RAM, single-threaded
  BLAS, R 4.4.2).** We benchmarked a single SeqDef call at 36 tree sizes from n = 100 to 100,000
  (100 replicates where inexpensive). Runtime and peak memory both scale as **O(n²)** — fitted
  log–log slopes of **2.09** (time) and **1.94** (memory) for n ≥ 1,000. A single call takes
  0.04 s at n = 1,000, 6.2 s at n = 10,000 and **808 s (13.5 min) at n = 100,000**.
- **Memory is the binding constraint, not runtime.** Peak memory grows quadratically — ~4.7 GB at
  n = 10,000 and **447 GB at n = 100,000**; we pushed the benchmark until memory exhaustion, and the
  n = 150,000 run exceeded an 800 GB allocation, consistent with the fit. (These O(n²) figures apply
  to the Gaussian and linear kernels, which use the dense distance matrix.)
- **For the default exponential kernel we now remove the O(n²) wall entirely.** Because exp(−λd)
  factorizes along tree paths, the phylogenetically weighted availability is computed by an exact
  O(n)-time, O(n)-memory tree traversal — no distance matrix is ever formed. It is identical to the
  dense result (≤ 1e-16) yet runs a single call at n = 100,000 in **0.06 s using 74 MB** (vs 808 s /
  447 GB dense) and at **n = 1,000,000 in < 1 s using ~0.4 GB**, on a laptop. The default kernel
  therefore scales *exactly* to whole-tree-of-life size, superseding the sparse-kernel approximation
  we had flagged as future work (Supplementary Fig. S4 now contrasts the two regimes).
- **Real data.** The 877-tip Chondrichthyes tree runs in **0.02 s** (single λ) / **0.21 s**
  (`auto_max`); the full 100-tree posterior with `auto_max` completes in ≈ **40 s**.
- We therefore state that SeqDef with the **default exponential kernel scales exactly to whole-tree-
  of-life size in O(n)**, while the Gaussian/linear kernels are O(n²) and efficient at clade-to-class
  scale (Supplementary Fig. S4). For those dense kernels, a truncated/sparse approximation remains a
  natural route to the very largest trees (future work).

### 1.5 Typographical problems in the mathematical notation
> *In Eq. 4, "Priority" (italic P) and "SeqDef" (italic f) collide with variable notation; set
> multi-letter labels upright.*

**Response.** Corrected. "Priority," "SeqDef" and all label subscripts (norm, median, depth) are
now upright (roman), with italics reserved for scalar variables (S, w, d, A, λ, n) and indices
(i, f). While doing this we caught and fixed a substantive omission: the displayed kernel equation
had **lost its λ** — it read *w*ᵢ = exp(−*d*ᵢ / *T*) — whereas the text and code use
*w*_fi = exp(−λ·*d*_fi / *T*_depth). The equation now matches the implementation, and the priority
formula is rendered unambiguously as Priority = SeqDef_median × 2^GE. We performed a full
consistency pass so every symbol is defined once and used identically in text, equations and code.

---

## Reviewer 2

### 2.1 Introduction conflates two distinct frameworks
> *The Introduction mixes distinctiveness/threat ranking and information-maximising sampling
> schemes as if opposed.*

**Response.** We rewrote the passage to present two **complementary** traditions —
distinctiveness/threat ranking (Faith's PD; EDGE/EDGE2; Soares et al. 2023 for fishes) and
information-maximising sampling schemes (Marini et al. 2022) — and position SeqDef as a bridge
that imports the "condition on what is already sampled" logic into the conservation-prioritization
tradition, weightable by any threat or value trait. We removed the oppositional phrasing and note
that EDGE2 itself conditions on relatives' risk, so the field is converging on the conditioning
logic SeqDef formalises for genomic data.

### 2.2 Extension to taxa lacking any molecular data (e.g., rtrees)
> *Use SeqDef to identify unsampled taxa that would most improve coverage, by inserting them
> (e.g., rtrees) and scoring the augmented tree.*

**Response.** Excellent suggestion; we added a Discussion paragraph. The current design prunes to
the tree∩data intersection, so only placed taxa are scored. Species lacking molecular data — common
for tropical lineages — can be grafted onto the phylogeny with taxonomic insertion tools (rtrees,
V.PhyloMaker, `phytools::add.species.to.genus`, PASTIS/TACT) and then scored; having no data and
often sitting in undersampled clades, they surface as high-deficiency targets. We note the caveat
(insertion adds placement/branch-length uncertainty) and that our posterior-tree framework already
integrates over uncertainty, so one could average SeqDef over stochastic insertions. We cite
Soares et al. (2023), who used exactly this insertion strategy for the fish tree of life, and flag
a planned `prune = FALSE` option so users can supply an augmented tree directly.

### 2.3 Definition of S and examples beyond 0/1
**Response.** Addressed together with Reviewer 1's 1.1: S is a general score on [0,1]; we give
worked continuous constructions in the Methods (data-type breadth and quality gradients). The binary
coding is retained for the main case study with its rationale stated.

### 2.4 Incorporating a taxonomic-uncertainty weight
**Response.** We note that a taxonomic-uncertainty weight slots naturally into the framework either
as an additional trait passed to `calc_priority` or as a per-tip down-weight on S, and we flag it
as a planned option (future work) rather than altering the core statistic. Combined with the
insertion-based extension (2.2), this would let users discount confidently for poorly placed tips.

### 2.5 Line-level comments
> *Line 40: see Soares et al., Conservation Biology.*

**Response.** Added. The intended reference is Soares et al. (2023), "Quantifying and overcoming
Darwinian shortfalls to conserve the fish tree of life," published in **Biological Conservation**
(the reviewer wrote "Conservation Biology" — we cite the correct journal; please confirm). It
supports both the Introduction framing and the new insertion-based Discussion paragraph.

> *Line 48: define "marginal information value."*

**Response.** Defined at first use: *by marginal information value we mean the expected gain in
genome-derived comparative information from adding a focal taxon's genome to the current set, given
the genomes already available in its relatives — high when little of the taxon's information is
already captured by sequenced close relatives.* We tie this to the formal object, the complement
of phylogenetically weighted availability (1 − A_f).

> *Lines 178–180: clarify the multiplicative ("doubling") effect.*

**Response.** Rewritten: "Because GE enters as an exponent of 2, each one-category increase in
extinction risk **doubles** a taxon's priority — a Critically Endangered taxon (GE = 4) receives
2⁴ = 16× the weight of a Least Concern taxon (GE = 0) with the same SeqDef — so priority rises
steeply for taxa that are simultaneously isolated from existing genomes and highly threatened."

### 2.6 Figure 1 should report the λ value
**Response.** Agreed (and seconded by the Associate Editor). The previous panel used illustrative
hand-set values; we **regenerated Figure 1 from an actual `SeqDef()` run** on the toy tree
(`auto_max` λ = 1.2) and now show **three λ values (1.2, 5, 15)** so the horizon effect is explicit:
*taxa9* (sister of the only sequenced tip) has SeqDef = 0.23, 0.63 and 0.95 as λ increases — its
"coverage" by a sequenced relative vanishes as the phylogenetic horizon shrinks.

---

## Editorial and compliance items

- **Data Availability.** We completed the Data Availability Statement: a citable Zenodo snapshot of
  the SeqDef code (DOI: [to be minted]), the NCBI assembly table and the analysis scripts; the
  posterior trees are cited to VertLife / Stein et al. (2018). IUCN Red List data are not
  redistributable, so we provide the Red List version, the exact query and the access date.
- **Author contributions.** CRediT contributions for all authors are provided with the revised
  submission [to be completed by the authors].
- **Reproducibility housekeeping.** Before public archiving we will remove and rotate the API keys
  embedded in the analysis script, complete the package `DESCRIPTION` (authors, ORCIDs, license)
  and citation, and add unit tests covering the three kernels and the λ-selection edge cases.
- **Dataset count.** We corrected the species count to **n = 877** (the value reproduced by the
  analysis pipeline; the previous "n = 850" was a transcription error).

*We thank the editor and reviewers for feedback that has measurably strengthened the paper.*
