# Response to Reviewers

**Manuscript JEB-2026-00140.R1** — *SeqDef: An R package for phylogenetically weighted prioritization of genome-sequencing targets* (previously "… genomic prioritization across the tree of life")

Reviewer comments are quoted in italics, followed by our response. Line numbers refer to the revised manuscript (clean version). All analyses were re-run for this revision with SeqDef v0.2.0 on all 100 posterior trees; scripts and outputs are in https://github.com/Andresdbp/SeqDef_project (folder `analyses/`, one script per result, listed in the README). Values in brackets are 2.5–97.5% quantiles across the 100 posterior trees.

---

## Response to the Handling Editor and the Associate Editor

We thank the Handling Editor, the Associate Editor and both reviewers. Reviewer 2's central point is correct, and we have revised the paper around it rather than defending our original wording.

**Does SeqDef provide novel information for prioritizing sequencing?** We now compare SeqDef directly with the two simpler measures Reviewer 2 proposed, in a new Results subsection (lines 408–430) and a new Supplementary Figure S4. In brief:

- On a time-calibrated tree, SeqDef assigns one score to every species in a clade without sequenced members. This holds for any measure built only from the tree and the list of sequenced species, including the two simpler measures. SeqDef therefore ranks *lineages*, and a second criterion (here IUCN status) chooses species within a lineage. We now say this in the Abstract, Methods, Results and Discussion.
- SeqDef agrees only moderately with distance to the nearest sequenced relative (ρ = 0.75 [0.63, 0.83] across species; 0.64 [0.52, 0.75] across clades) and weakly with whether the family already has an assembly (ρ = 0.43 [0.34, 0.53]). The top-ranked target would be its family's first assembly in 99 of 100 trees, but this is true of only 40% of the top 25 and 24% of the top 50; 68% of unsequenced species belong to families that already have an assembly, where a family rule gives no guidance.
- Nearest-relative distance is a special case of SeqDef: it is the large-λ limit of the exponential kernel (ρ = 0.97 at λ = 10 and 0.995 at λ = 25). At the selected λ, SeqDef also accounts for how many sequenced relatives lie within the phylogenetic horizon.

We have scaled our claims to match. SeqDef is now presented as a tunable generalization of nearest-relative distance, with an exact linear-time algorithm, support for posterior samples of trees, and any decision trait. Its distinctive contribution appears below the very top of the list, where it largely agrees with the simpler rules (lines 471–482).

**Other changes.** We corrected Figure 3 (ties) and now report robustness statistics for SeqDef alone, for Priority, and for Priority within IUCN categories. We corrected the Brownian-kernel statement (it is exactly the linear kernel with λ = 0.5) and made the *Sphyrna lewini* comparison like for like. We added a data-assembly paragraph and a sensitivity analysis on all 1,192 tips. Figure 1 now shows values computed by the package at three values of λ. We corrected the min–max normalization rationale, and the package now returns the unnormalized availability. We added the repository URL and unit tests comparing the O(n) and full-matrix calculations. The title is narrower, a paragraph on tree requirements has been added, and Reviewer 1's continuous-S example has been added. Figures 1, 3, S1 and S2 were replaced, and Figure S4 is new.

**Journal requirements.** The Data Availability Statement now gives the repository URLs and data sources. A CRediT statement for all authors is provided on the title page, and a tracked-changes manuscript is uploaded. Code and data are archived on GitHub (package: https://github.com/Andresdbp/SeqDef, release v0.2.0; analyses: https://github.com/Andresdbp/SeqDef_project). We are happy to deposit a versioned snapshot with a DOI if the editor prefers.

---

## Reviewer 1

> *… it would be useful to provide one illustrative example of how such information could be translated into a value between 0 and 1. For example, the authors could indicate whether measures such as BUSCO completeness, assembly level, or transcriptome availability could be used individually or combined.*

We thank the reviewer. We added one worked construction to the Methods (lines 75–80) and did not change the case study. In this construction, S is the product of the best available resource level (0, none; 0.25, transcriptome only; 0.5, contig-level; 0.75, scaffold-level; 1, chromosome-level assembly) and its BUSCO completeness. For example, a scaffold-level assembly with 92% BUSCO completeness gives S = 0.69, and a chromosome-level assembly with 96% gives 0.96. Either component can be used alone, and we note that a binary coding remains preferable when quality metrics are missing for many taxa.

---

## Reviewer 2 — major concerns

### 1. SeqDef cannot distinguish species within a clade lacking sequenced members

> *On an ultrametric (time-calibrated) tree, take any clade with no sequenced members … every tip in such a clade receives exactly the same A_f, and therefore the same SeqDef.*

The reviewer is correct, and we thank them for stating it so clearly. We now state the property in the Methods (lines 111–117) and in the package documentation. We also quantify it (lines 410–412): the 822 unsequenced species fall into 178 [168, 189] clades without a sequenced member, and within these clades SeqDef varies by at most 2 × 10⁻⁸, which is rounding in the tree file. The same property holds for nearest-relative distance and for the family rule, so it reflects the information available (the tree and the list of sequenced species) rather than the kernel.

> *Consequence of the case study … The abstract's claim that "SeqDef identifies" the top target should be reframed accordingly.*

Agreed. SeqDef identifies the lineage (Centrophoridae, which has no assembly), and the IUCN weight selects *C. atromarginatus*, the only Critically Endangered centrophorid, within it. We rewrote the Abstract, which now also states the relation to nearest-relative distance and stays within the 250-word limit (lines 14–19), Results (lines 385–387), Discussion (lines 442–444) and Conclusions (lines 536–538).

> *Figure 3: … the counts shown in the figure sum to 107 …*

Thank you for catching this. The discrepancy came from ties. In 85 trees *C. atromarginatus* is the unique top-ranked species. In the remaining 15 trees, the top position is an exact tie among all Critically Endangered members of one clade without sequenced members: the eight CR species of *Rhinobatos* and *Acroteriobatus* in 14 trees, and six CR *Squatina* species in one tree. This is exactly the behaviour the reviewer describes. The R1 figure counted every tied species, and because tied scores differed only by floating-point noise, the ties were also split unevenly. The revised Figure 3 counts one unit per tree (85 + 14 + 1 = 100) and labels tied sets. The text and caption explain how ties were detected (relative tolerance 10⁻⁶) and how many trees had them (lines 358–364; caption lines 367–371). The 85/100 result is unchanged.

> *Figure S2 shows many tips sharing identical SeqDef values …*

Yes, this is the same property. We now say so in the Results (lines 410–412) and in the caption of Figure S2 (lines 554–555).

> *Comparison with simpler measures … report the correlation between SeqDef and distance to the nearest sequenced relative; report the proportion of top-ranked targets that would be the first assembly in their family; show specific cases where SeqDef ranks species differently …*

We added a Methods description (lines 285–288), a Results subsection (lines 408–430) and Supplementary Figure S4.

| Requested statistic | Value |
| --- | --- |
| ρ(SeqDef, nearest-relative distance), 822 unsequenced species | 0.75 [0.63, 0.83] |
| Same, one value per unsequenced clade | 0.64 [0.52, 0.75] |
| ρ(SeqDef, family lacks an assembly) | 0.43 [0.34, 0.53] |
| Top-ranked target would be its family's first assembly | 99 of 100 trees |
| Top 10 / 25 / 50 targets that would be first in their family | 90% / 40% / 24% |
| Unsequenced species whose family already has an assembly | 68% |
| Top-25 overlap, Priority vs nearest-distance × 2^GE | 12 [11, 18] |
| Same top-ranked species | 68 of 100 trees |
| Pairs of unsequenced clades ordered differently (MCC tree) | 24% |

*How the measures relate.* Nearest-relative distance is the large-λ limit of SeqDef. As λ grows, the closest sequenced relative dominates A_f, and ρ rises from 0.46 at λ = 1 to 0.75 at the selected λ, 0.97 at λ = 10 and 0.995 at λ = 25 (Figure S4B). At the selected λ, SeqDef also counts how many sequenced relatives lie within the horizon.

*Where they differ, and why.* *Callorhinchus capensis* has a sequenced congener (nearest-relative distance 0.02 of tree depth) but SeqDef = 0.77, because *C. milii* is the only genome in the holocephalan neighbourhood. The carcharhinid *Lamiopsis temminckii* is farther from its nearest genome (0.21) but has seven sequenced relatives within half the tree depth, and its SeqDef is 0.33. *Narcine* species score 0.95 although Narcinidae has an assembly, because that assembly is distant. Gurgesiellidae has no assembly but scores 0.45–0.49, because it is nested among sequenced skates. The family rule cannot make either distinction.

*Which ranking is preferable* depends on the purpose, and we now say so (lines 471–482). If the goal is a close reference for one particular taxon (e.g. for read mapping or reference-guided assembly), nearest-relative distance, which SeqDef reproduces with a large λ, is the right criterion. If the goal is to fill the largest gaps in comparative coverage across a clade, a lineage whose only nearby reference is a single genome is less well covered than one surrounded by many genomes, and SeqDef at the selected λ captures this. Because the measures largely agree at the very top of the list, we no longer claim that SeqDef finds a different first target. We describe its value as a single, tunable and reproducible ranking of every taxon that updates as genomes accumulate, scales to millions of tips, and integrates phylogenetic uncertainty.

### 2. The *Sphyrna lewini* example compares unlike things

We agree. On the MCC tree, *S. lewini* ranks 40th on EDGE and 103rd on Priority (threat-weighted vs threat-weighted), and 269th on ED and 779th on SeqDef (threat-free vs threat-free). Most of the R1 "divergence" therefore disappears in a fair comparison, and we have removed that claim. (A small note: 0.365 × 16 = 5.84 rather than 0.58, but the reviewer's conclusion holds.)

The like-for-like comparison gives a clearer illustration of complementarity (lines 400–405). Species high on EDGE that already have genome assemblies drop sharply on Priority, e.g. *Sphyrna mokarran* (EDGE rank 20 → Priority 223) and *Galeorhinus galeus* (6 → 107). Species in unsequenced lineages rise, e.g. *Centrophorus granulosus* (182 → 70) and *Urolophus orarius* (191 → 76). The Methods now state which scores are compared with which (line 284).

The purple points in Figure S2 are the five high-EDGE species whose congener is sequenced and whose EDGE–SeqDef rank gap is largest: *Carcharhinus borneensis*, *S. lewini*, *S. corona*, *C. signatus* and *S. tudes*. All five are now labelled in the figure and named in the caption (lines 552–554).

### 3. The robustness analysis mainly reflects the IUCN weight

The reviewer is right that the R1 statistics were computed on Priority, where the IUCN term dominates. We now report every statistic for SeqDef alone, for Priority, and for Priority within each IUCN category, on all 100 posterior trees (lines 335–342; revised Figure S1). Fig. S1 in R1 used a 20-tree subsample.

| | SeqDef alone | Priority |
| --- | --- | --- |
| ρ with auto_max across λ ∈ [1, 25] (median; minimum) | 0.82; 0.78 | 0.94; 0.89 |
| Top-10 agreement with auto_max (Jaccard, ties included), median across λ | 0.25 | 0.72 |
| ρ, auto_max vs by_genus | 1.00 [0.97, 1.00] | 1.00 [0.99, 1.00] |
| ρ with the exponential: Gaussian / linear / Brownian | 0.91 / 0.92 / 0.83 | 0.88 / 0.94 / 0.78 |

Within IUCN categories, the median ρ across λ is 0.82 (LC), 0.73 (NT), 0.84 (VU), 0.87 (EN) and 0.75 (CR). Across species, ρ(Priority, GE) = 0.86 [0.82, 0.87], whereas ρ(Priority, SeqDef) = 0.28 [0.24, 0.33]. For 95% [93, 96] of species pairs in different IUCN categories, the Priority order simply follows the category. We now state plainly that, with base 2, the IUCN weight orders taxa between categories and SeqDef orders them within categories. We recommend reporting SeqDef alongside Priority, and we note that `calc_priority(base = …)` lets users reduce the weight of the trait (lines 274–277; lines 398–399).

**Top of the ranking at low λ.** Agreed; we replaced "highly stable". Whole-ranking agreement remains high, but the top of the SeqDef-only list depends on λ: the top-10 Jaccard index is 0.14 at λ = 1–2, 1.0 at the selected λ, and 0.25 (median) across [1, 25]. We now advise users to examine top targets across a band of λ. The values in Figure S1B are now medians across all 100 trees, with 2.5–97.5% bands as in panel A (in R1 they were means across 20 trees). A new panel C shows that *C. atromarginatus* is top-ranked in 69–95% of trees at every λ.

**Brownian kernel.** The reviewer is right, and we thank them for the derivation. With T = root age, C_fi = T − d_fi/2 and C_ff = C_ii = T, so the Brownian correlation is 1 − d_fi/(2T), which is exactly the linear kernel with λ = 0.5 on our depth-normalized scale. We confirmed this numerically (identical on all 100 trees up to rounding, maximum difference in SeqDef 1.6 × 10⁻⁸) and with a package unit test. We corrected both statements (lines 149–151; lines 496–501). We now attribute the Brownian kernel's lower agreement (Priority ρ = 0.78; top target in 23 of 100 trees) to its wide horizon (small λ), consistent with the sensitivity at small λ, rather than to a different model of information decay (lines 351–354).

### 4. Data provenance

We added a data-assembly paragraph (lines 230–245) covering each point:

- **Name matching.** Names were matched exactly, after removing extra whitespace and replacing spaces with underscores, without synonym resolution. Of the 1,192 tips, 980 matched a Red List assessment and 212 did not, mostly because of genus-level revisions since 2018; 167 of those match a Red List synonym. 237 assessed species are absent from the trees.
- **Data Deficient and Not Evaluated.** Species assessed as Data Deficient (102 tips) or Extinct (1) were excluded because no extinction-risk weight can be assigned. No tips were Not Evaluated.
- **Imputed tips.** From Stein et al. (2018, Supplementary Table 1): of the 877 analysed species, 498 were placed with molecular data and 379 by taxonomic constraints. *C. atromarginatus* was placed with molecular data. Among the Figure 3 species, *Rhinobatos albomaculatus*, *R. irvinei*, *R. lionotus*, *Squatina argentina* and *S. oculata* are imputed, and the others were placed with molecular data. Of the 55 sequenced species, 49 were placed with molecular data and 6 were imputed.
- **NCBI query date and definition of an assembly.** A species was coded S = 1 if the NCBI Assembly database held at least one genome-assembly record for it, of any assembly level and in GenBank or RefSeq, released on or before 30 December 2025, the date the case-study table was compiled. We re-ran this query on 23 November 2026 (search term "Chondrichthyes[Organism]", restricted by release date). It reproduced the original coding for 876 of the 877 species; the exception is *Harriotta haeckeli*, whose assembly NCBI files under its synonym *H. raleighana*, a separate tip in the tree.
- **Tips coded S = 1:** 55 of 877 (1–9 assemblies per species). 30 of 63 families, including Centrophoridae, have none.
- **IUCN version.** Latest global assessments (Red List version 2025-2), accessed 10 December 2025.

Because DD and name-mismatched tips were pruned before scoring, any genome they carry was ignored. We therefore repeated the analysis scoring SeqDef on all 1,192 tips, with S mapped through NCBI and Red List synonyms. This gave 60 tips with an assembly, six of them outside the main analysis (*Rhina ancylostoma*, *Stegostoma fasciatum*, *Manta birostris*, *Dasyatis sabina*, *D. akajei* and *Torpedo nobiliana*). Priority was then formed for all 1,005 assessed tips. The results were essentially unchanged: Priority ρ with the main analysis = 0.985 [0.945, 0.994], and *C. atromarginatus* was the unique top-ranked species in 97 of 100 trees (lines 389–392).

### 5. Statements to be corrected

- **d/T_depth spans [0, 2].** Corrected (lines 179–181). The package now defines T as the maximum root-to-tip distance, which equals the root age on ultrametric trees and is also defined for other trees.
- **"95% CI".** Replaced by "2.5–97.5% quantiles across the 100 posterior trees" in the Methods (lines 257–258), the Figure 2 caption (lines 327–328) and the Results.
- **"A significant methodological advance".** Removed. We now say that summarizing over posterior trees follows established practice (e.g. Gumbs et al., 2023) (lines 467–469).
- **Title.** Narrowed to "SeqDef: An R package for phylogenetically weighted prioritization of genome-sequencing targets". We also added a Discussion paragraph on tree requirements (lines 513–519). The package stops with an explicit message when branch lengths are missing, handles polytomies exactly (unit-tested), and applies the kernels to patristic distances on non-ultrametric trees, in which case the within-clade equality described in point 1 no longer holds exactly.

---

## Reviewer 2 — minor comments

**1. Figure 1.** We agree and thank the reviewer. Figure 1 now shows scores computed by the package on the toy tree at λ = 1, 3 and 10 (lines 215–218). It shows the horizon effect the editor asked about: the sister of the sequenced tip scores 0.21, 0.45 and 0.86 at the three values, while the most distant tips approach 1.

**2. Repository and test.** The package is at https://github.com/Andresdbp/SeqDef (release v0.2.0) and the analyses at https://github.com/Andresdbp/SeqDef_project. Both are cited in the Abstract (lines 22–23), Results, Conclusions and Data Availability Statement. The package's testthat suite now checks, on ultrametric and non-ultrametric random trees with binary and continuous S and with polytomies, that the O(n) exponential and Brownian traversals equal the full-matrix calculation to 10⁻¹⁰, and that the Brownian kernel equals the linear kernel at λ = 0.5 on ultrametric trees. `R CMD check` passes without errors, warnings or notes. This is described in the Methods (lines 190–192).

**3. Min–max normalization.** The reviewer is right; the R1 rationale was wrong. We now state that min–max scaling is a within-tree rescaling that aids ranking and display. It makes scores *less* comparable across trees or data snapshots, because each score is relative to that tree's best- and least-covered tips (lines 93–98). `SeqDef()` now also returns the unnormalized A_f (`$availability`) for such comparisons. When all A_f are equal (e.g. no tip has data), the denominator is zero: the package assigns SeqDef = 1 to every tip, which treats all tips as equally and maximally deficient, and issues a warning.

---

We thank both reviewers again. Reviewer 2's comments led to a clearer statement of what SeqDef does and does not do, and to corrections that make the paper more accurate.
