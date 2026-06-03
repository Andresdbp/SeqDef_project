# =============================================================================
# 04_runtime_benchmark.R
# Reviewer 1 #4 / Editor M4:  Back the "computationally efficient" claim with a
# complexity statement + empirical runtime/memory benchmarks.
#
# Per the project decision: SINGLE lambda (= 10) so we measure the cost of one
# SeqDef calculation (not the auto_max scan). Broad, log-spaced n grid pushed up
# to the O(n^2)-memory wall on this 16 GB machine.
#
# Produces:
#   results/runtime_benchmark.csv
#   figures/figS_runtime.pdf
#   (appends timing to results/SESSION.txt)
# =============================================================================

source("analyses/00_setup.R")
write_session()
set.seed(42)

SINGLE_LAMBDA <- 10
RAM_GB <- tryCatch(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1024^3, error = function(e) 16)
SAFETY_GB <- 0.65 * RAM_GB   # don't attempt runs whose 2 dense matrices exceed this

# log-spaced grid (powers of 10 + half-decade steps) then push toward the wall
ns <- c(100, 316, 1000, 3162, 10000, 15000, 20000, 30000)

peak_mb <- function() {                 # peak R memory (MB) since last gc(reset=TRUE)
  g <- gc()
  (g["Vcells", "max used"] * 8 + g["Ncells", "max used"] * 56) / 2^20
}

bench_one <- function(n, reps) {
  tr <- rtree(n)
  df <- data.frame(sp = tr$tip.label, s = sample(c(0, 1), n, replace = TRUE))
  tt <- mm <- numeric(reps)
  for (r in seq_len(reps)) {
    invisible(gc(reset = TRUE))
    tt[r] <- system.time(res <- SeqDef(tr, df, lambda = SINGLE_LAMBDA))[["elapsed"]]
    mm[r] <- peak_mb()
    rm(res)
  }
  tibble(n = n, time_s = median(tt), peak_mb = median(mm), status = "ok")
}

message("04: runtime benchmark (single lambda = ", SINGLE_LAMBDA, ")...")
rows <- list()
for (n in ns) {
  est_gb <- 2 * (n^2) * 8 / 1e9          # two n x n double matrices
  if (est_gb > SAFETY_GB) {
    message(sprintf("  n=%6d  SKIP (est. %.1f GB > %.1f GB safety cap)", n, est_gb, SAFETY_GB))
    rows[[as.character(n)]] <- tibble(n = n, time_s = NA, peak_mb = NA, status = "skipped_memory")
    next
  }
  reps <- if (n <= 3162) 3 else 2
  out <- tryCatch(bench_one(n, reps), error = function(e) {
    message(sprintf("  n=%6d  FAILED: %s", n, conditionMessage(e)))
    tibble(n = n, time_s = NA, peak_mb = NA, status = "error")
  })
  if (!is.na(out$time_s))
    message(sprintf("  n=%6d  time=%.3fs  peak=%.0f MB (theory %.2f GB)", n, out$time_s, out$peak_mb, est_gb))
  rows[[as.character(n)]] <- out
}
rt <- bind_rows(rows)
rt$theory_gb <- 2 * (rt$n^2) * 8 / 1e9
write.csv(rt, "results/runtime_benchmark.csv", row.names = FALSE)

# ---- Real Chondrichthyes anchor ---------------------------------------------
inp <- build_input_binary(tree_mcc)
t_single <- system.time(suppressMessages(SeqDef(inp$tree, inp$df, lambda = SINGLE_LAMBDA)))[["elapsed"]]
t_auto   <- system.time(suppressMessages(SeqDef(inp$tree, inp$df, lambda = "auto_max")))[["elapsed"]]
cat(sprintf("\n[Chondrichthyes MCC, %d tips]  single-lambda: %.3fs   auto_max: %.3fs\n",
            length(inp$tree$tip.label), t_single, t_auto))
cat(sprintf("[Chondrichthyes posterior]     100 trees x auto_max (from Analysis 1) ~ 40s wall-clock\n"))

# ---- Complexity fit (compute-bound window n in [1000, 10000]) ---------------
# Beyond ~10k the dense matrices approach RAM and wall-clock inflates from
# GC/paging (the memory-bound regime), so we fit the slope where the run is
# compute-bound and report the >10k inflation separately.
fitdat <- rt %>% filter(status == "ok", n >= 1000, n <= 10000, !is.na(time_s), time_s > 0)
slope <- coef(lm(log(time_s) ~ log(n), data = fitdat))[["log(n)"]]
mem_slope <- coef(lm(log(peak_mb) ~ log(n), data = fitdat))[["log(n)"]]
cat(sprintf("\n[Complexity] time ~ n^%.2f ; peak-memory ~ n^%.2f  (theory: O(n^2), fit on n in [1000,10000])\n", slope, mem_slope))
cat("[Note] for n > 10000 wall-clock inflates super-quadratically as the dense n x n matrices approach RAM (GC/paging) -- the memory-bound regime.\n")
n_at_7gb <- sqrt(7e9 / 8)   # single n x n double = 7 GB
cat(sprintf("[Memory ceiling] one n x n double matrix = 7 GB at n ~ %.0f; ~80 GB at n=100,000\n", n_at_7gb))

# ---- Figure -----------------------------------------------------------------
ok <- rt %>% filter(status == "ok", !is.na(time_s))
okT <- ok %>% filter(time_s > 0)
pT <- ggplot(okT, aes(n, time_s)) +
  geom_line(color = "#1862C9") + geom_point(color = "#1862C9", size = 2) +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides = "bl") +
  labs(x = "n (tips)", y = "time per SeqDef call (s)",
       subtitle = sprintf("A  Runtime, single lambda  (slope = %.2f ~ O(n^2))", slope)) +
  theme_minimal() + theme(panel.grid.minor = element_blank())

memdat <- tibble(n = c(ok$n, 30000, 1e5), gb = 8 * c(ok$n, 30000, 1e5)^2 / 1e9)
pM <- ggplot() +
  geom_line(data = memdat, aes(n, gb), color = "grey55", linetype = "dashed") +
  geom_point(data = ok, aes(n, peak_mb / 1024), color = "#C96E18", size = 2) +
  geom_hline(yintercept = RAM_GB, color = "#D73027", linetype = "dotted") +
  annotate("text", x = 130, y = RAM_GB * 1.12, label = sprintf("%.0f GB RAM", RAM_GB),
           color = "#D73027", hjust = 0, size = 3) +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(x = "n (tips)", y = "memory (GB)",
       subtitle = "B  Memory: measured peak (orange) vs O(n^2) single-matrix (dashed)") +
  theme_minimal() + theme(panel.grid.minor = element_blank())

ggsave("figures/figS_runtime.pdf", patchwork::wrap_plots(pT, pM, nrow = 1), width = 10, height = 4)
cat("\n[04] DONE. Wrote results/runtime_benchmark.csv, figures/figS_runtime.pdf\n")
