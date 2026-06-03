# =============================================================================
# grace_runtime_compute.R  -- runtime benchmark for the 3 TB / 80-core Grace node
# COMPUTE ONLY (writes a CSV; the figure is made locally).
#
# Measures the cost of ONE SeqDef call (single lambda, exponential kernel) on
# simulated trees. Single-threaded BLAS for a clean O(n^2) scaling curve.
# Dense geometric n grid for a smooth curve + large n to probe the ceiling.
#
# Run via SLURM (see grace_runtime.slurm). No data needed (trees are simulated).
# =============================================================================

suppressWarnings(suppressMessages({
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
}))
suppressPackageStartupMessages(library(ape))
source("function.R")        # kernel-enabled SeqDef() (scp'd from local)

set.seed(42)
SINGLE_LAMBDA <- 10

# memory budget from SLURM (MB) if available, else assume ~2.4 TB usable
mem_mb_alloc <- as.numeric(Sys.getenv("SLURM_MEM_PER_NODE", "2457600"))
MEM_CAP_BYTES <- 0.80 * mem_mb_alloc * 1024^2
TIME_CAP_S    <- 900       # stop growing n once a single median call exceeds this

peak_mb <- function() { g <- gc(); (g["Vcells","max used"]*8 + g["Ncells","max used"]*56)/2^20 }

reps_for <- function(n) if (n <= 12000) 100 else if (n <= 30000) 25 else if (n <= 60000) 8 else 3

# dense geometric grid 100..12000 + larger probes
small <- unique(round(exp(seq(log(100), log(12000), length.out = 30))))
large <- c(15000, 20000, 30000, 50000, 75000, 100000, 150000, 200000)
ns <- sort(unique(c(small, large)))

cat(sprintf("mem cap = %.0f GB ; n grid = %d points (max %d)\n", MEM_CAP_BYTES/1024^3, length(ns), max(ns)))

rows <- list()
for (n in ns) {
  est_bytes <- 3 * (n^2) * 8                     # cophenetic + weight + transient
  if (est_bytes > MEM_CAP_BYTES) { cat(sprintf("n=%7d SKIP (est %.0f GB > cap)\n", n, est_bytes/1024^3)); next }
  reps <- reps_for(n)
  tt <- numeric(reps); mm <- numeric(reps)
  ok <- TRUE
  for (r in seq_len(reps)) {
    tr <- rtree(n)
    df <- data.frame(sp = tr$tip.label, s = sample(c(0, 1), n, replace = TRUE))
    invisible(gc(reset = TRUE))
    tm <- tryCatch(system.time(SeqDef(tr, df, lambda = SINGLE_LAMBDA, kernel = "exponential"))[["elapsed"]],
                   error = function(e) { cat("  ERR n=", n, ": ", conditionMessage(e), "\n"); NA })
    if (is.na(tm)) { ok <- FALSE; break }
    tt[r] <- tm; mm[r] <- peak_mb()
    rm(tr, df); invisible(gc())
  }
  if (!ok) { rows[[as.character(n)]] <- data.frame(n=n, reps=reps, time_s=NA, time_lo=NA, time_hi=NA, peak_mb=NA, status="error"); next }
  med <- median(tt)
  rows[[as.character(n)]] <- data.frame(
    n = n, reps = reps, time_s = med,
    time_lo = as.numeric(quantile(tt, .05)), time_hi = as.numeric(quantile(tt, .95)),
    peak_mb = median(mm), status = "ok")
  cat(sprintf("n=%7d reps=%3d time=%9.4fs [%.4f, %.4f] peak=%7.0f MB\n", n, reps, med, quantile(tt,.05), quantile(tt,.95), median(mm)))
  if (med > TIME_CAP_S) { cat("  reached time cap; stopping larger n.\n"); break }
}
res <- do.call(rbind, rows)
res$theory_gb <- 8 * (res$n^2) / 1e9
write.csv(res, "results/runtime_benchmark_grace.csv", row.names = FALSE)
cat("\nDONE -> results/runtime_benchmark_grace.csv\n")
sink("results/SESSION_grace.txt"); cat("Grace runtime benchmark\n"); print(Sys.info()); cat("SLURM_MEM_PER_NODE(MB)=", Sys.getenv("SLURM_MEM_PER_NODE"), "\nSLURM_CPUS=", Sys.getenv("SLURM_CPUS_ON_NODE"), "\n"); print(sessionInfo()); sink()
