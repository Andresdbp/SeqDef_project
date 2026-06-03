# =============================================================================
# run_all.R  -- reproduce every revision analysis (run from project root).
#   Rscript analyses/run_all.R
# Heavy steps cache to results/*.rds; delete those to force a clean recompute.
# Never sources analysis.R (API keys); loads cached data directly.
# =============================================================================
scripts <- sprintf("analyses/%s", c(
  "01_lambda_stability.R",
  "02_kernel_comparison.R",
  "03_edge_benchmark.R",
  "04_runtime_benchmark.R",
  "05_fig1_regen.R",
  "06_continuous_S.R"
))
for (s in scripts) {
  message("\n========== ", s, " ==========")
  source(s, echo = FALSE)
}
message("\nAll analyses complete. See results/SUMMARY.md for headline numbers.")
