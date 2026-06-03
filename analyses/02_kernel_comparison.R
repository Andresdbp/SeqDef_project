# =============================================================================
# 02_kernel_comparison.R
# Reviewer 1 #2b / Editor M2b:  Justify the exponential kernel by showing that
# SeqDef rankings are robust to the kernel choice.
#
# Uses the new kernel= argument (exponential / gaussian / linear), each with its
# own auto_max lambda, on the MCC tree + a posterior sample.
#
# Produces:
#   results/kernel_comparison.csv
#   figures/figS_kernel.pdf
# =============================================================================

source("analyses/00_setup.R")
top_set <- function(pri) names(pri)[pri == max(pri, na.rm = TRUE)]
kernels <- c("exponential", "gaussian", "linear")

run_kernels <- function(phy_raw, label) {
  inp <- build_input_binary(phy_raw)
  if (is.null(inp)) return(NULL)
  ge <- ge_for_tree(inp$tree)
  res <- lapply(kernels, function(k) {
    r <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = "auto_max", kernel = k))
    list(seqdef = r$seqdef, pri = calc_priority(r, ge, "exponential", base = 2), lambda = r$lambda)
  })
  names(res) <- kernels
  common <- names(res$exponential$seqdef)
  sp <- function(a, b) suppressWarnings(cor(a[common], b[common], method = "spearman"))
  tibble(
    tree = label,
    lambda_exp = res$exponential$lambda, lambda_gauss = res$gaussian$lambda, lambda_lin = res$linear$lambda,
    rho_pri_exp_gauss = sp(res$exponential$pri, res$gaussian$pri),
    rho_pri_exp_lin   = sp(res$exponential$pri, res$linear$pri),
    rho_pri_gauss_lin = sp(res$gaussian$pri,   res$linear$pri),
    rho_sd_exp_gauss  = sp(res$exponential$seqdef, res$gaussian$seqdef),
    rho_sd_exp_lin    = sp(res$exponential$seqdef, res$linear$seqdef),
    target_exp   = TARGET %in% top_set(res$exponential$pri),
    target_gauss = TARGET %in% top_set(res$gaussian$pri),
    target_lin   = TARGET %in% top_set(res$linear$pri)
  )
}

# MCC tree (store priority vectors for the scatter)
inp_mcc <- build_input_binary(tree_mcc); ge_mcc <- ge_for_tree(inp_mcc$tree)
mcc_pri <- sapply(kernels, function(k) {
  r <- suppressMessages(SeqDef(inp_mcc$tree, inp_mcc$df, lambda = "auto_max", kernel = k))
  calc_priority(r, ge_mcc, "exponential", base = 2)[inp_mcc$tree$tip.label]
})

# MCC + posterior sample
set.seed(42)
samp_ids <- sort(sample(seq_along(chond), 20))
message("02: kernel comparison on MCC + 20 posterior trees...")
kc <- bind_rows(
  run_kernels(tree_mcc, "MCC"),
  map_dfr(samp_ids, ~ run_kernels(chond[[.x]], paste0("post_", .x)), .progress = TRUE)
)
write.csv(kc, "results/kernel_comparison.csv", row.names = FALSE)

post <- kc %>% filter(tree != "MCC")
cat("\n[Kernel comparison] (posterior sample, n =", nrow(post), ")\n")
cat(sprintf("  Spearman rho of PRIORITY rankings vs exponential:\n    gaussian: median=%.3f [%.3f, %.3f]\n    linear  : median=%.3f [%.3f, %.3f]\n",
            median(post$rho_pri_exp_gauss), quantile(post$rho_pri_exp_gauss,.025), quantile(post$rho_pri_exp_gauss,.975),
            median(post$rho_pri_exp_lin),   quantile(post$rho_pri_exp_lin,.025),   quantile(post$rho_pri_exp_lin,.975)))
cat(sprintf("  target in top set: exp %d/%d | gauss %d/%d | linear %d/%d\n",
            sum(post$target_exp), nrow(post), sum(post$target_gauss), nrow(post), sum(post$target_lin), nrow(post)))

# ---- Figure -----------------------------------------------------------------
long <- post %>%
  select(tree, `exp vs gaussian` = rho_pri_exp_gauss, `exp vs linear` = rho_pri_exp_lin,
         `gaussian vs linear` = rho_pri_gauss_lin) %>%
  pivot_longer(-tree, names_to = "pair", values_to = "rho")

pA <- ggplot(long, aes(pair, rho)) +
  geom_boxplot(fill = "#75A1FF", alpha = 0.5, outlier.size = 0.6, width = 0.55) +
  coord_cartesian(ylim = c(min(0.7, min(long$rho)), 1)) +
  labs(x = NULL, y = expression(Spearman~rho~"(Priority)"), subtitle = "A  Cross-kernel ranking agreement (posterior)") +
  theme_minimal() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 15, hjust = 1))

scat <- as.data.frame(mcc_pri); scat$tip <- rownames(scat); scat$is_target <- scat$tip == TARGET
pB <- ggplot(scat, aes(exponential, gaussian)) +
  geom_point(color = "grey60", size = 0.8, alpha = 0.6) +
  geom_point(data = subset(scat, is_target), color = "#D73027", size = 2.6) +
  geom_text(data = subset(scat, is_target), aes(label = "C. atromarginatus"),
            color = "#D73027", hjust = 1.1, vjust = 0.3, size = 3, fontface = "italic") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(x = "Priority (exponential)", y = "Priority (gaussian)", subtitle = "B  MCC priority: exponential vs gaussian") +
  theme_minimal() + theme(panel.grid.minor = element_blank())

ggsave("figures/figS_kernel.pdf", patchwork::wrap_plots(pA, pB, nrow = 1), width = 10, height = 4)
cat("\n[02] DONE. Wrote results/kernel_comparison.csv, figures/figS_kernel.pdf\n")
