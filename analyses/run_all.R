# =============================================================================
# run_all.R  -- reproduce every analysis (run from project root).
#   Rscript analyses/run_all.R
# Heavy steps cache to results/*.rds; delete those to force a clean recompute.
# Never sources analysis.R (API keys); loads cached data directly.
# =============================================================================
scripts <- sprintf("analyses/%s", c(
  # ---- R1 analyses ----
  "figS1_lambda_stability.R",          # posterior cache + lambda/method CSVs (R1)
  "figS2_edge_benchmark.R",            # Fig. S2 + EDGE CSVs (needs picante)
  "figS3_runtime.R",                   # Fig. S3 (reads results/dense_kernel_grace.csv)
  "kernel_comparison.R",               # kernel-agreement values (R1, Priority only)
  # ---- R2 analyses (JEB-2026-00140.R1 decision, Nov 2026) ----
  "fig1_toy_tree.R",                   # Fig. 1 with computed values (Minor 1)
  "fig3_topological_robustness.R",     # Fig. 3, tie-aware (Major 1)
  "simple_measures_comparison.R",      # SeqDef vs nearest-relative distance / family gap; Fig. S4 (Major 1)
  "stability_seqdef_vs_priority.R",    # SeqDef-alone stability, within-IUCN, kernels; Fig. S1 (Major 3)
  "data_provenance.R",                 # counts for the data-assembly paragraph (Major 4)
  "brownian_comparison.R",             # (re-run) Brownian kernel on all 100 trees
  "ncbi_map_all_tips.R",               # NCBI snapshot mapped onto all 1,192 tips (Major 4)
  "full_tree_sensitivity.R",           # SeqDef on all 1,192 tips (Major 4 sensitivity)
  "priority_table.R"                   # full ranked table cited in the Results
))
for (s in scripts) {
  message("\n========== ", s, " ==========")
  source(s, echo = FALSE)
}
# Not run here:
#  * grace_dense_kernels.R  -- O(n^2) HPC points for Fig. S3 (large-memory node)
#  * ncbi_query_all_tips.R  -- needs internet + rentrez; regenerates
#    data/ncbi_assembly_species_summary.tsv from NCBI (committed copy: 2026-11-23)
message("\nAll analyses complete.")
