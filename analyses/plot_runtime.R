# =============================================================================
# plot_runtime.R  -- make figures/figS_runtime.pdf from the Grace benchmark CSV.
# Points only (no connecting lines); dotted O(n^2) reference; NO hardware RAM line.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

rt <- read.csv("results/runtime_benchmark_grace.csv")
ok <- rt %>% filter(status == "ok", !is.na(time_s), time_s > 0)

# Best-fit O(n^2) reference: time = a * n^2 (least squares on the compute-bound points)
cb <- ok %>% filter(n >= 1000)
a_t <- exp(mean(log(cb$time_s)  - 2 * log(cb$n)))
a_m <- exp(mean(log(cb$peak_mb / 1024) - 2 * log(cb$n)))
slope_t <- coef(lm(log(time_s) ~ log(n), data = cb))[["log(n)"]]
ref <- data.frame(n = 10^seq(log10(min(ok$n)), log10(max(ok$n)), length.out = 200))
ref$time <- a_t * ref$n^2
ref$gb   <- a_m * ref$n^2

pT <- ggplot(ok, aes(n, time_s)) +
  geom_line(data = ref, aes(n, time), linetype = "dotted", color = "grey35", linewidth = 0.6) +
  geom_point(color = "#1862C9", size = 1.5) +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(x = "n (tips)", y = "time per SeqDef call (s)",
       subtitle = bquote("A  Runtime, single " * lambda * "  (dotted = O(" * n^2 * "); fitted slope " * .(sprintf("%.2f", slope_t)) * ")")) +
  theme_minimal() + theme(panel.grid.minor = element_blank())

pM <- ggplot(ok, aes(n, peak_mb / 1024)) +
  geom_line(data = ref, aes(n, gb), linetype = "dotted", color = "grey35", linewidth = 0.6) +
  geom_point(color = "#C96E18", size = 1.5) +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(x = "n (tips)", y = "peak memory (GB)",
       subtitle = bquote("B  Peak memory  (dotted = O(" * n^2 * "))")) +
  theme_minimal() + theme(panel.grid.minor = element_blank())

ggsave("figures/figS_runtime.pdf", pT + pM, width = 10, height = 4)
cat(sprintf("figS_runtime.pdf written. n range %d-%d (%d points); fitted time slope %.2f\n",
            min(ok$n), max(ok$n), nrow(ok), slope_t))
