# =============================================================================
# figS1_lambda_stability.R
# Reviewer 1 #2a / Editor M2a:  Does lambda change the PRIORITIZATION (not just
# the variance)?  And how do the two lambda-selection methods compare on outcomes?
#
# Produces:
#   results/posterior_priority.rds     (cached per-tree auto_max & by_genus runs)
#   results/lambda_rank_stability.csv
#   results/method_comparison.csv
#   figures/figS1_lambda_stability.pdf
# =============================================================================

source("analyses/00_setup.R")

# ---- helper: priority winner set (handles ties) -----------------------------
top_set <- function(pri) names(pri)[pri == max(pri, na.rm = TRUE)]

# =============================================================================
# PART A.  Per-tree auto_max & by_genus across all 100 posterior trees
#          (cached; also reproduces the Fig-3 robustness distribution)
# =============================================================================
post_path <- "results/posterior_priority.rds"
if (file.exists(post_path)) {
  posterior <- readRDS(post_path)
} else {
  message("PART A: running auto_max + by_genus on 100 posterior trees...")
  posterior <- map(seq_along(chond), function(i) {
    inp <- build_input_binary(chond[[i]])
    if (is.null(inp)) return(NULL)
    ge <- ge_for_tree(inp$tree)

    r_auto  <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = "auto_max"))
    r_genus <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = "by_genus"))

    p_auto  <- calc_priority(r_auto,  ge, model = "exponential", base = 2)
    p_genus <- calc_priority(r_genus, ge, model = "exponential", base = 2)

    list(tree_id = i,
         lambda_auto = r_auto$lambda, lambda_genus = r_genus$lambda,
         seqdef_auto = r_auto$seqdef, seqdef_genus = r_genus$seqdef,
         pri_auto = p_auto, pri_genus = p_genus)
  }, .progress = TRUE)
  posterior <- compact(posterior)
  saveRDS(posterior, post_path)
}

# ---- Robustness distribution (auto_max winners; reproduces Fig 3) -----------
winners_auto <- map(posterior, ~ top_set(.x$pri_auto))
win_tbl <- tibble(sp = unlist(winners_auto)) %>% count(sp, sort = TRUE)
n_trees <- length(posterior)
target_wins <- win_tbl$n[win_tbl$sp == TARGET]
target_wins <- ifelse(length(target_wins) == 0, 0, target_wins)
cat(sprintf("\n[Robustness] trees analysed: %d | %s is top in %d trees\n",
            n_trees, TARGET, target_wins))
print(head(win_tbl, 6))

# =============================================================================
# PART B.  Method head-to-head: auto_max vs by_genus (across 100 trees)
# =============================================================================
method_cmp <- map_dfr(posterior, function(p) {
  common <- intersect(names(p$pri_auto), names(p$pri_genus))
  rho <- suppressWarnings(cor(p$pri_auto[common], p$pri_genus[common], method = "spearman"))
  tibble(
    tree_id      = p$tree_id,
    lambda_auto  = p$lambda_auto,
    lambda_genus = p$lambda_genus,
    rho_priority = rho,
    auto_is_target  = TARGET %in% top_set(p$pri_auto),
    genus_is_target = TARGET %in% top_set(p$pri_genus),
    agree_top = length(intersect(top_set(p$pri_auto), top_set(p$pri_genus))) > 0
  )
})
write.csv(method_cmp, "results/method_comparison.csv", row.names = FALSE)

cat(sprintf("\n[Method head-to-head]\n  lambda_auto  median=%.2f [%.2f, %.2f]\n  lambda_genus median=%.2f [%.2f, %.2f]\n",
            median(method_cmp$lambda_auto),  quantile(method_cmp$lambda_auto, .025),  quantile(method_cmp$lambda_auto, .975),
            median(method_cmp$lambda_genus), quantile(method_cmp$lambda_genus, .025), quantile(method_cmp$lambda_genus, .975)))
cat(sprintf("  Spearman rho(auto, genus) priority: median=%.3f [%.3f, %.3f]\n",
            median(method_cmp$rho_priority), quantile(method_cmp$rho_priority, .025), quantile(method_cmp$rho_priority, .975)))
cat(sprintf("  auto_max top-1 == target : %d/%d\n", sum(method_cmp$auto_is_target),  nrow(method_cmp)))
cat(sprintf("  by_genus top-1 == target : %d/%d\n", sum(method_cmp$genus_is_target), nrow(method_cmp)))
cat(sprintf("  methods agree on top set : %d/%d\n", sum(method_cmp$agree_top), nrow(method_cmp)))

# =============================================================================
# PART C.  Lambda rank-stability:  fix lambda on a grid and compare the
#          resulting rankings against each tree's own auto_max solution.
# =============================================================================
set.seed(42)
samp_ids <- sort(sample(seq_along(posterior), min(20, length(posterior))))
lambda_grid <- seq(1, 25, 0.5)

overlap_k <- function(a, b, k = 10) {
  length(intersect(names(sort(a, decreasing = TRUE))[1:k],
                   names(sort(b, decreasing = TRUE))[1:k]))
}

message("PART C: lambda rank-stability over grid [1,25]...")
stab <- map_dfr(samp_ids, function(i) {
  p   <- posterior[[i]]
  inp <- build_input_binary(chond[[p$tree_id]])
  ge  <- ge_for_tree(inp$tree)
  sd_auto  <- p$seqdef_auto
  pri_auto <- p$pri_auto
  common   <- names(sd_auto)

  map_dfr(lambda_grid, function(lam) {
    r   <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = lam))
    pri <- calc_priority(r, ge, model = "exponential", base = 2)
    sd  <- r$seqdef[common]; pr <- pri[common]
    tibble(
      tree_id = p$tree_id, lambda = lam,
      rho_seqdef   = suppressWarnings(cor(sd, sd_auto[common],  method = "spearman")),
      rho_priority = suppressWarnings(cor(pr, pri_auto[common], method = "spearman")),
      top10_seqdef   = overlap_k(r$seqdef, sd_auto, 10),
      top10_priority = overlap_k(pri, pri_auto, 10),
      top1_is_target = TARGET %in% top_set(pri)
    )
  })
}, .progress = TRUE)
write.csv(stab, "results/lambda_rank_stability.csv", row.names = FALSE)

# Summaries
stab_summary <- stab %>%
  group_by(lambda) %>%
  summarise(
    rho_pri_mean = mean(rho_priority), rho_pri_lo = quantile(rho_priority, .025), rho_pri_hi = quantile(rho_priority, .975),
    rho_sd_mean  = mean(rho_seqdef),
    top10_pri_mean = mean(top10_priority),
    frac_target_top1 = mean(top1_is_target), .groups = "drop"
  )
# lambda band where target is the top-1 in a majority of sampled trees
maj_band <- stab_summary$lambda[stab_summary$frac_target_top1 >= 0.5]
cat(sprintf("\n[Lambda rank-stability over [1,25]]\n  rho(priority vs auto_max): median across grid = %.3f (min %.3f)\n",
            median(stab_summary$rho_pri_mean), min(stab_summary$rho_pri_mean)))
cat(sprintf("  rho(seqdef   vs auto_max): median across grid = %.3f (min %.3f)\n",
            median(stab_summary$rho_sd_mean), min(stab_summary$rho_sd_mean)))
if (length(maj_band))
  cat(sprintf("  target is top-1 in >=50%% of sampled trees for lambda in [%.1f, %.1f]\n",
              min(maj_band), max(maj_band)))

# =============================================================================
# PART D.  Figure (Supplementary Fig. S1)
# =============================================================================
lam_med <- median(method_cmp$lambda_auto)
lam_lo  <- quantile(method_cmp$lambda_auto, .025)
lam_hi  <- quantile(method_cmp$lambda_auto, .975)
band <- list(
  annotate("rect", xmin = lam_lo, xmax = lam_hi, ymin = -Inf, ymax = Inf, fill = "#D73027", alpha = 0.12),
  geom_vline(xintercept = lam_med, color = "#D73027", linetype = "longdash"))  # auto-selected lambda (explained in caption)
tagTR <- function(lab) annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.4, vjust = 1.5, fontface = "bold", size = 5)
p1 <- ggplot(stab_summary, aes(lambda, rho_pri_mean)) + band +
  geom_ribbon(aes(ymin = rho_pri_lo, ymax = rho_pri_hi), fill = "#1862C9", alpha = 0.2) +
  geom_line(color = "#1862C9", linewidth = 1) +
  coord_cartesian(ylim = c(0.5, 1)) + tagTR("A") +
  labs(x = expression(lambda), y = expression("Ranking agreement with auto_max (Spearman " * rho * ")")) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
p2 <- ggplot(stab_summary, aes(lambda, top10_pri_mean)) + band +
  geom_point(color = "#C96E18", size = 1.7) +
  scale_y_continuous(breaks = 0:10, minor_breaks = NULL, limits = c(0, 10)) + tagTR("B") +
  labs(x = expression(lambda), y = "Top-10 overlap with auto_max (of 10)") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("figures/figS1_lambda_stability.pdf", patchwork::wrap_plots(p1, p2, nrow = 1), width = 10, height = 4)

cat("\n[figS1] DONE. Wrote results/{posterior_priority.rds, method_comparison.csv, lambda_rank_stability.csv}, figures/figS1_lambda_stability.pdf\n")
