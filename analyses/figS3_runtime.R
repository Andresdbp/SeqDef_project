# =============================================================================
# figS3_runtime.R
# Computational cost of SeqDef across the four kernels (Supplementary Fig. S3).
#
#   - Exponential & Brownian: exact O(n) tree traversals (no n x n matrix);
#     timed on this machine up to n = 1e6.
#   - Gaussian & linear: O(n^2) dense computation (full cophenetic matrix);
#     timed on this machine up to the laptop memory limit, and on a large-memory
#     HPC node (results/dense_kernel_grace.csv, written by grace_dense_kernels.R)
#     up to n = 1.5e5.
#
# Produces:
#   results/traversal_benchmark.csv     (exponential + brownian, O(n))
#   results/dense_kernel_benchmark.csv  (gaussian + linear, O(n^2), this machine)
#   figures/figS3_runtime.pdf
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(ggplot2); library(dplyr); library(patchwork) })
source("function.R")
set.seed(42)
peak_mb <- function() { g <- gc(); (g["Vcells","max used"]*8 + g["Ncells","max used"]*56) / 2^20 }

# --- O(n) traversals: exponential + brownian, to n = 1e6 ----------------------
# Build each random tree once and time the traversal repeatedly (its cost is
# topology-independent); take the median.
ns <- sort(unique(c(round(exp(seq(log(100), log(1e6), length.out = 16))), 1e4, 1e5, 1e6)))
reps_for <- function(n) if (n <= 1e4) 30 else if (n <= 1e5) 15 else 7
message("Benchmarking O(n) traversals (exponential + brownian) to n = 1e6 ...")
trav <- do.call(rbind, lapply(ns, function(n) {
  reps <- reps_for(n)
  tr <- rtree(n); s <- sample(c(0, 1), n, replace = TRUE); td <- max(branching.times(tr))
  invisible(gc(reset = TRUE))
  te <- sapply(seq_len(reps), function(i) system.time(.seqdef_exp_avail(tr, s, 10, td))[["elapsed"]])
  me <- peak_mb()
  invisible(gc(reset = TRUE))
  tb <- sapply(seq_len(reps), function(i) system.time(.seqdef_bm_avail(tr, s))[["elapsed"]])
  mb <- peak_mb()
  rm(tr, s); invisible(gc())
  rbind(data.frame(n = n, time_s = median(te), peak_mb = me, kernel = "exponential"),
        data.frame(n = n, time_s = median(tb), peak_mb = mb, kernel = "brownian"))
}))
write.csv(trav, "results/traversal_benchmark.csv", row.names = FALSE)

# --- O(n^2) dense: gaussian + linear, to the laptop memory limit --------------
ns_d <- c(100, 300, 1000, 2000, 4000, 6000, 8000, 10000, 12000)
reps_d <- function(n) if (n <= 2000) 5 else if (n <= 8000) 3 else 2
message("Benchmarking O(n^2) dense (gaussian + linear) on this machine ...")
lapd <- do.call(rbind, lapply(ns_d, function(n) {
  out <- list()
  for (k in c("gaussian", "linear")) {
    tr <- rtree(n); df <- data.frame(taxa = tr$tip.label, score = sample(c(0, 1), n, replace = TRUE))
    tt <- numeric(reps_d(n))
    for (r in seq_len(reps_d(n))) { invisible(gc(reset = TRUE))
      tt[r] <- system.time(suppressMessages(SeqDef(tr, df, lambda = 10, kernel = k)))[["elapsed"]] }
    out[[k]] <- data.frame(n = n, time_s = median(tt), peak_mb = peak_mb(), kernel = k)
    rm(tr, df); invisible(gc())
  }
  do.call(rbind, out)
}))
write.csv(lapd, "results/dense_kernel_benchmark.csv", row.names = FALSE)

# --- O(n^2) at HPC scale (gaussian + linear; from grace_dense_kernels.R) -------
hpc <- read.csv("results/dense_kernel_grace.csv")

# --- Figure: four kernels; circles = laptop, squares = large-memory server ----
trav$source <- "laptop"; lapd$source <- "laptop"; hpc$source <- "hpc"
all <- rbind(trav[, c("n","time_s","peak_mb","kernel","source")],
             lapd[, c("n","time_s","peak_mb","kernel","source")],
             hpc[,  c("n","time_s","peak_mb","kernel","source")])
all$kernel <- factor(all$kernel, levels = c("exponential","gaussian","linear","brownian"))
all$source <- factor(all$source, levels = c("laptop","hpc"))
all <- all[order(match(all$kernel, c("brownian","linear","gaussian","exponential"))), ]
cols4 <- c(exponential = "#1B9E77", gaussian = "#D95F02", linear = "#7570B3", brownian = "#E7298A")
klab  <- c(exponential = "Exponential", gaussian = "Gaussian", linear = "Linear", brownian = "Brownian")
scl <- list(scale_color_manual(values = cols4, labels = klab, name = NULL),
            scale_shape_manual(values = c(laptop = 16, hpc = 15),
                               labels = c(laptop = "Laptop", hpc = "Large-memory server"), name = NULL))
tagTR <- function(lab) annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.4, vjust = 1.5, fontface = "bold", size = 5)

ex <- subset(all, kernel == "exponential"); ga <- subset(all, kernel == "gaussian")
xrA <- 10^seq(2, 6, length.out = 120); xrD <- 10^seq(2, log10(1.5e5), length.out = 90)
aT1 <- median(ex$time_s[ex$n>=1e3]/ex$n[ex$n>=1e3]);         aT2 <- median(ga$time_s[ga$n>=1e3]/ga$n[ga$n>=1e3]^2)
aM1 <- median((ex$peak_mb[ex$n>=1e4]/1024)/ex$n[ex$n>=1e4]); aM2 <- median((ga$peak_mb[ga$n>=1e3]/1024)/ga$n[ga$n>=1e3]^2)
refT <- rbind(data.frame(n = xrA, y = aT1*xrA, g = "On"), data.frame(n = xrD, y = aT2*xrD^2, g = "On2"))
refM <- rbind(data.frame(n = xrA, y = aM1*xrA, g = "On"), data.frame(n = xrD, y = aM2*xrD^2, g = "On2"))

pT <- ggplot(subset(all, time_s > 0), aes(n, time_s, color = kernel, shape = source)) +
  geom_line(data = refT, aes(n, y, group = g), inherit.aes = FALSE, linetype = "dotted", color = "grey55") +
  geom_point(size = 2.9, alpha = 0.6) + scl +
  annotate("text", x = 4e5,   y = aT1*4e5*6,     label = "O(n)",   parse = TRUE, color = "grey35", size = 3.5) +
  annotate("text", x = 4.5e4, y = aT2*4.5e4^2*4, label = "O(n^2)", parse = TRUE, color = "grey35", size = 3.5) +
  tagTR("A") + scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(x = "Number of tips (n)", y = "Runtime (s)") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
pM <- ggplot(all, aes(n, peak_mb/1024, color = kernel, shape = source)) +
  geom_line(data = refM, aes(n, y, group = g), inherit.aes = FALSE, linetype = "dotted", color = "grey55") +
  geom_point(size = 2.9, alpha = 0.6) + scl +
  annotate("text", x = 4e5,   y = aM1*4e5*7,     label = "O(n)",   parse = TRUE, color = "grey35", size = 3.5) +
  annotate("text", x = 4.5e4, y = aM2*4.5e4^2*4, label = "O(n^2)", parse = TRUE, color = "grey35", size = 3.5) +
  tagTR("B") + scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(x = "Number of tips (n)", y = "Memory use (GB)") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("figures/figS3_runtime.pdf", (pT + pM) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 11, height = 4.4)
cat("[figS3] DONE -> figures/figS3_runtime.pdf, results/{traversal_benchmark.csv, dense_kernel_benchmark.csv}\n")
