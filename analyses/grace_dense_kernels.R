# =============================================================================
# grace_dense_kernels.R -- benchmark the O(n^2) kernels (GAUSSIAN + LINEAR) on a
# large-memory node, pushing n past the laptop limit. COMPUTE ONLY (writes CSV).
#
# These two kernels build the full n x n distance matrix, so peak memory ~6*n^2*8
# bytes (=447 GB at n=1e5). With ~3 TB we can reach n ~ 200-250k.
#
# Run via SLURM (grace_dense_kernels.slurm). Needs function.R scp'd alongside.
# =============================================================================
suppressWarnings(suppressMessages({
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
}))
suppressPackageStartupMessages(library(ape))
source("function.R")
set.seed(42)
SINGLE_LAMBDA <- 10

# Edit this grid to trade data points against node time. Each large-n point is
# expensive: ~13 min at n=1e5, ~30 min at 1.5e5, ~55 min at 2e5 (single-threaded).
ns <- c(25000, 50000, 100000, 150000)        # add 200000 / 250000 if mem + time allow
reps_for <- function(n) if (n <= 50000) 2 else 1

mem_mb_alloc  <- as.numeric(Sys.getenv("SLURM_MEM_PER_NODE", "2900000"))
MEM_CAP_BYTES <- 0.85 * mem_mb_alloc * 1024^2
peak_mb <- function() { g <- gc(); (g["Vcells","max used"]*8 + g["Ncells","max used"]*56)/2^20 }

cat(sprintf("mem cap = %.0f GB\n", MEM_CAP_BYTES/1024^3))
rows <- list()
for (n in ns) {
  est <- 6 * (n^2) * 8                        # ~6 n-sized double arrays at peak
  if (est > MEM_CAP_BYTES) { cat(sprintf("n=%7d SKIP (est %.0f GB > cap)\n", n, est/1024^3)); next }
  for (k in c("gaussian", "linear")) {
    reps <- reps_for(n); tt <- numeric(reps); mm <- NA; ok <- TRUE
    for (r in seq_len(reps)) {
      tr <- rtree(n); df <- data.frame(sp = tr$tip.label, s = sample(c(0,1), n, replace = TRUE))
      invisible(gc(reset = TRUE))
      tm <- tryCatch(system.time(SeqDef(tr, df, lambda = SINGLE_LAMBDA, kernel = k))[["elapsed"]],
                     error = function(e) { cat("  ERR", k, n, ":", conditionMessage(e), "\n"); NA })
      if (is.na(tm)) { ok <- FALSE; break }
      tt[r] <- tm; mm <- peak_mb(); rm(tr, df); invisible(gc())
    }
    if (!ok) next
    rows[[paste(k, n)]] <- data.frame(n = n, time_s = median(tt), peak_mb = mm, kernel = k)
    cat(sprintf("%-9s n=%7d time=%9.1fs peak=%8.0f MB (%.0f GB)\n", k, n, median(tt), mm, mm/1024))
    flush.console()
  }
}
res <- do.call(rbind, rows)
write.csv(res, "results/dense_kernel_grace.csv", row.names = FALSE)
cat("DONE -> results/dense_kernel_grace.csv\n")
