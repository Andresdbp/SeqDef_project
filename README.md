# SeqDef — analyses for the JEB paper

This repository holds the analysis pipeline, data, results and figures for

> Barboza, A. (2026). SeqDef: an R package for phylogenetically weighted prioritization of
> genome-sequencing targets. *Journal of Evolutionary Biology* (in revision, JEB-2026-00140).

The R package itself is at **https://github.com/Andresdbp/SeqDef** (release v0.2.0).
`function.R` here is an exact copy of the package code (`R/SeqDef.R` + `R/calc_priority.R`), so the
analyses run the published implementation.

## What SeqDef does, in one paragraph

For each tip of a time-calibrated phylogeny, SeqDef measures how much of that taxon's genomic
information is already represented by sequenced relatives: it is the complement of the
phylogenetically weighted availability `A_f = sum_i w_fi S_i`, with weights from a distance-decay
kernel whose parameter λ sets a phylogenetic horizon. All tips of a clade without sequenced members
get the same score, so SeqDef ranks lineages; `calc_priority()` combines it with a trait (here IUCN
extinction risk) to rank species. As λ grows, SeqDef converges on the ranking by distance to the
nearest sequenced relative. The default exponential kernel runs in exact O(n) time and memory.

## Repository layout

```
function.R        SeqDef() + calc_priority() (copy of the package code, v0.2.0)
analyses/         one script per figure or in-text result (see below); run_all.R runs them all
data/             input data and their sources (data/README.md); the IUCN file is not redistributed
results/          result tables (CSV), cached intermediate objects (RDS), logs/
figures/          manuscript figures (PDF) and figures/png/ (300 dpi)
revision/         responses to reviewers and change records for each revision round
analysis.R        original end-to-end data download + Figure 4 (needs API tokens from environment variables)
validation.R      early checks of the statistic
```

## Reproducing the analyses

From the repository root, with `data/iucn_assessment_data.rds` re-created as described in
`data/README.md`:

```sh
Rscript analyses/run_all.R
```

| Script | Produces |
| --- | --- |
| `figS1_lambda_stability.R` | posterior cache (`results/posterior_priority.rds`), λ-method comparison |
| `figS2_edge_benchmark.R` | EDGE/EDGE2 benchmark, like-for-like ranks, Fig. S2 |
| `figS3_runtime.R` | runtime/memory benchmark, Fig. S3 (large-n dense points: `grace_dense_kernels.R` on an HPC node) |
| `kernel_comparison.R` | four-kernel Priority comparison |
| `fig1_toy_tree.R` | Fig. 1 (toy tree at λ = 1, 3, 10) |
| `fig3_topological_robustness.R` | Fig. 3 (top-ranked species per tree, ties handled) |
| `simple_measures_comparison.R` | SeqDef vs nearest-relative distance and family rule, Fig. S4 |
| `stability_seqdef_vs_priority.R` | stability for SeqDef, Priority and within IUCN categories; kernels; Fig. S1 |
| `data_provenance.R` | name matching, exclusions, imputed tips, NCBI/IUCN details (`results/provenance.txt`) |
| `brownian_comparison.R` | Brownian kernel = linear kernel at λ = 0.5; Brownian vs exponential |
| `ncbi_map_all_tips.R` | NCBI snapshot mapped onto all 1,192 tips |
| `full_tree_sensitivity.R` | analysis repeated on all 1,192 tips |
| `priority_table.R` | full ranked table (`results/priority_ranking_posterior.csv`) |
| `ncbi_query_all_tips.R` | re-runs the NCBI query (needs internet; not run by `run_all.R`) |

All analyses use all 100 posterior trees; random steps use `set.seed(42)`.

## Licence

Code: MIT (see `LICENSE`). Results and figures: CC BY 4.0. Input data keep their providers' terms
(`data/README.md`); IUCN Red List data are not redistributed here.
