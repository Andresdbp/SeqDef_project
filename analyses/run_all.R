# =============================================================================
# run_all.R  -- reproduce every revision analysis (run from project root).
#   Rscript analyses/run_all.R
# Heavy steps cache to results/*.rds; delete those to force a clean recompute.
# Never sources analysis.R (API keys); loads cached data directly.
# =============================================================================
scripts <- sprintf("analyses/%s", c(
  "figS1_lambda_stability.R",  # Fig. S1 + lambda/method CSVs
  "figS2_edge_benchmark.R",    # Fig. S2 + EDGE CSVs
  "figS3_runtime.R",           # Fig. S3 (reads results/dense_kernel_grace.csv)
  "kernel_comparison.R",       # kernel-agreement values (reported in text)
  "brownian_comparison.R"      # Brownian-comparison values (reported in text)
))
for (s in scripts) {
  message("\n========== ", s, " ==========")
  source(s, echo = FALSE)
}
# The O(n^2) HPC points in Fig. S3 are produced separately on a large-memory
# node by grace_dense_kernels.R (see grace_dense_kernels.slurm) and read back
# from results/dense_kernel_grace.csv.
message("\nAll analyses complete. See results/SUMMARY.md for headline numbers.")
