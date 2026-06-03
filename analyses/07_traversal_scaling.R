# =============================================================================
# 07_traversal_scaling.R
# Benchmark the exact O(n) exponential traversal against the dense O(n^2) path,
# and draw the comparison figure for the revised M4 / "why exponential" argument.
#
# Produces:
#   results/traversal_benchmark.csv
#   figures/figS_runtime.pdf   (overwrites: now dense O(n^2) vs traversal O(n))
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(ggplot2); library(dplyr); library(patchwork) })
source("function.R")
set.seed(42)

peak_mb <- function() { g <- gc(); (g["Vcells","max used"]*8 + g["Ncells","max used"]*56) / 2^20 }

# --- Benchmark the exponential traversal (single lambda), O(n) -> up to 1e6 ----
ns <- unique(round(exp(seq(log(100), log(1e6), length.out = 16))))
message("Benchmarking exponential traversal to n = 1e6 ...")
trav <- do.call(rbind, lapply(ns, function(n) {
  reps <- if (n <= 1e4) 5 else if (n <= 1e5) 3 else 1
  tt <- mm <- numeric(reps)
  for (r in seq_len(reps)) {
    tr <- rtree(n); s <- sample(c(0,1), n, replace = TRUE); td <- max(branching.times(tr))
    invisible(gc(reset = TRUE))
    tt[r] <- system.time(.seqdef_exp_avail(tr, s, 10, td))[["elapsed"]]
    mm[r] <- peak_mb(); rm(tr, s)
  }
  data.frame(n = n, time_s = median(tt), peak_mb = median(mm), method = "traversal (exponential)")
}))
write.csv(trav, "results/traversal_benchmark.csv", row.names = FALSE)
cat(sprintf("traversal: n=1e6 -> %.2fs, %.0f MB\n",
            trav$time_s[trav$n == max(trav$n)], trav$peak_mb[trav$n == max(trav$n)]))

# --- Dense path (from the Grace benchmark) ------------------------------------
dense <- read.csv("results/runtime_benchmark_grace.csv")
dense <- subset(dense, status == "ok")[, c("n","time_s","peak_mb")]
dense$method <- "dense (gaussian/linear)"
both <- rbind(trav[, c("n","time_s","peak_mb","method")], dense)
cols <- c("dense (gaussian/linear)" = "#C96E18", "traversal (exponential)" = "#2CA25F")

# --- Reference slopes: O(n^2) anchored to dense, O(n) anchored to traversal ----
aN2 <- exp(mean(log(dense$time_s[dense$n>=1000]) - 2*log(dense$n[dense$n>=1000])))
aN1 <- exp(mean(log(trav$time_s[trav$n>=1e3])    - 1*log(trav$n[trav$n>=1e3])))
mN2 <- exp(mean(log(dense$peak_mb[dense$n>=1000]/1024) - 2*log(dense$n[dense$n>=1000])))
mN1 <- exp(mean(log(trav$peak_mb[trav$n>=1e4]/1024)    - 1*log(trav$n[trav$n>=1e4])))
xr  <- 10^seq(2, 6, length.out = 100)
ref <- rbind(data.frame(n=xr, y=aN2*xr^2, lab="O(n^2)"), data.frame(n=xr, y=aN1*xr, lab="O(n)"))
refm<- rbind(data.frame(n=xr, y=mN2*xr^2, lab="O(n^2)"), data.frame(n=xr, y=mN1*xr, lab="O(n)"))

pT <- ggplot(subset(both, time_s > 0), aes(n, time_s, color = method)) +
  geom_line(data = ref, aes(n, y, group = lab), inherit.aes = FALSE, linetype = "dotted", color = "grey55") +
  geom_point(size = 1.7) + scale_color_manual(values = cols, name = NULL) +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(x = "n (tips)", y = "time per call (s)", subtitle = "A  Runtime: dense O(n^2) vs exponential traversal O(n)") +
  theme_minimal() + theme(panel.grid.minor = element_blank())

pM <- ggplot(both, aes(n, peak_mb/1024, color = method)) +
  geom_line(data = refm, aes(n, y, group = lab), inherit.aes = FALSE, linetype = "dotted", color = "grey55") +
  geom_point(size = 1.7) + scale_color_manual(values = cols, name = NULL) +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(x = "n (tips)", y = "peak memory (GB)", subtitle = "B  Peak memory: dense O(n^2) vs traversal O(n)") +
  theme_minimal() + theme(panel.grid.minor = element_blank())

ggsave("figures/figS_runtime.pdf", (pT + pM) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 11, height = 4.4)
cat("[07] DONE -> figures/figS_runtime.pdf, results/traversal_benchmark.csv\n")
