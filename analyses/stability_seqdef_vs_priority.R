# =============================================================================
# stability_seqdef_vs_priority.R
# R2 Reviewer 2, Major 3: the robustness statistics in the R1 manuscript were
# computed on Priority = SeqDef x 2^GE, where the IUCN term (1..16) dominates.
# Here we report every stability statistic for (i) SeqDef alone, (ii) Priority,
# and (iii) Priority within each IUCN category, plus how much of the Priority
# ranking is determined by GE vs SeqDef.  All 100 posterior trees are used
# (the R1 Fig. S1 used a 20-tree subsample).  Ties are handled explicitly.
#
# Produces:
#   results/stability_lambda_grid.csv   (tree x lambda)
#   results/stability_summary.csv       (one row per comparison x quantity)
#   results/kernel_comparison_r2.csv    (4 kernels, SeqDef + Priority, 100 trees)
#   figures/figS1_lambda_stability.pdf  (revised: SeqDef and Priority; ribbons on both panels)
# =============================================================================
source("analyses/00_setup.R")
posterior <- readRDS("results/posterior_priority.rds")

sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
TOL <- 1e-6
top_set  <- function(x) names(x)[x >= max(x) * (1 - TOL)]
# tie-inclusive top-k: everything at least as large as the k-th largest value
topk_set <- function(x, k) { v <- sort(x, decreasing = TRUE)[k]; names(x)[x >= v * (1 - TOL)] }
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
cats <- names(risk_index)

# Statistic bundle comparing a candidate (sd, pri) with a reference (sd0, pri0)
compare <- function(sd, sd0, pri, pri0, ge, s) {
  u <- s == 0
  within <- sapply(cats, function(cc) { k <- ge == risk_index[cc]; if (sum(k) > 2) sp(pri[k], pri0[k]) else NA })
  c(rho_seqdef = sp(sd, sd0), rho_seqdef_unseq = sp(sd[u], sd0[u]), rho_priority = sp(pri, pri0),
    setNames(within, paste0("rho_pri_within_", cats)),
    jac10_seqdef = jacc(topk_set(sd, 10), topk_set(sd0, 10)),
    jac10_priority = jacc(topk_set(pri, 10), topk_set(pri0, 10)),
    top1_same = as.numeric(length(intersect(top_set(pri), top_set(pri0))) > 0))
}

# ---- (1) lambda grid on all 100 trees ---------------------------------------
lambda_grid <- seq(1, 25, 0.5)
message("lambda grid x 100 trees ...")
grid <- map_dfr(posterior, function(p) {
  inp <- build_input_binary(chond[[p$tree_id]]); ge <- ge_for_tree(inp$tree); s <- inp$df$score
  nm <- inp$tree$tip.label; sd0 <- p$seqdef_auto[nm]; pri0 <- p$pri_auto[nm]
  map_dfr(lambda_grid, function(l) {
    r <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = l)); sd <- r$seqdef[nm]; pri <- sd * 2^ge[nm]
    as_tibble(as.list(c(tree_id = p$tree_id, lambda = l, compare(sd, sd0, pri, pri0, ge[nm], s),
                        target_top1 = as.numeric(TARGET %in% top_set(pri)))))
  })
}, .progress = TRUE)
write.csv(grid, "results/stability_lambda_grid.csv", row.names = FALSE)

# ---- (2) auto_max vs by_genus (cached runs) ----------------------------------
meth <- map_dfr(posterior, function(p) {
  inp <- build_input_binary(chond[[p$tree_id]]); ge <- ge_for_tree(inp$tree); nm <- inp$tree$tip.label
  as_tibble(as.list(compare(p$seqdef_genus[nm], p$seqdef_auto[nm], p$pri_genus[nm], p$pri_auto[nm], ge[nm], inp$df$score)))
})

# ---- (3) kernels: each with its own auto_max lambda; plus Brownian == linear(0.5)
message("kernels x 100 trees ...")
kern <- map_dfr(posterior, function(p) {
  inp <- build_input_binary(chond[[p$tree_id]]); ge <- ge_for_tree(inp$tree); nm <- inp$tree$tip.label
  sd0 <- p$seqdef_auto[nm]; pri0 <- p$pri_auto[nm]
  out <- map_dfr(c("gaussian", "linear", "brownian"), function(k) {
    r <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = "auto_max", kernel = k)); sd <- r$seqdef[nm]
    as_tibble(as.list(c(compare(sd, sd0, sd * 2^ge[nm], pri0, ge[nm], inp$df$score), lambda = r$lambda))) %>% mutate(kernel = k)
  })
  bm  <- suppressMessages(SeqDef(inp$tree, inp$df, kernel = "brownian"))$seqdef[nm]
  l05 <- suppressMessages(SeqDef(inp$tree, inp$df, kernel = "linear", lambda = 0.5))$seqdef[nm]
  out %>% mutate(tree_id = p$tree_id, max_abs_diff_bm_vs_linear05 = max(abs(bm - l05)))
}, .progress = TRUE)
write.csv(kern, "results/kernel_comparison_r2.csv", row.names = FALSE)

# ---- (4) How much of Priority is GE? -----------------------------------------
dom <- map_dfr(posterior, function(p) {
  inp <- build_input_binary(chond[[p$tree_id]]); ge <- ge_for_tree(inp$tree); nm <- inp$tree$tip.label
  pri <- p$pri_auto[nm]; sd <- p$seqdef_auto[nm]; g <- ge[nm]
  # share of species pairs with different GE whose Priority order follows GE
  set.seed(42); i <- sample(length(nm), 2e5, TRUE); j <- sample(length(nm), 2e5, TRUE)
  diffge <- g[i] != g[j]
  tibble(rho_pri_ge = sp(pri, g), rho_pri_seqdef = sp(pri, sd), rho_seqdef_ge = sp(sd, g),
         pairs_follow_ge = mean(sign(pri[i] - pri[j])[diffge] == sign(g[i] - g[j])[diffge]))
})

# ---- summaries -----------------------------------------------------------------
q <- function(v) sprintf("%.2f [%.2f, %.2f]", median(v, na.rm = TRUE), quantile(v, .025, na.rm = TRUE), quantile(v, .975, na.rm = TRUE))
gsum <- grid %>% group_by(lambda) %>% summarise(across(-tree_id, list(m = ~median(.x, na.rm = TRUE),
                                                                        lo = ~quantile(.x, .025, na.rm = TRUE),
                                                                        hi = ~quantile(.x, .975, na.rm = TRUE),
                                                                        mean = ~mean(.x, na.rm = TRUE))), .groups = "drop")
cat("\n================ STABILITY: SeqDef vs Priority ================\n")
cat("[lambda in 1..25, per-lambda median across 100 trees; then median (min) across the grid]\n")
for (v in c("rho_seqdef", "rho_seqdef_unseq", "rho_priority", paste0("rho_pri_within_", cats), "jac10_seqdef", "jac10_priority", "top1_same")) {
  m <- gsum[[paste0(v, "_m")]]; cat(sprintf("  %-22s median %.3f  min %.3f\n", v, median(m), min(m)))
}
cat("  low-lambda detail (lambda = 1, 2, 3, 5, 10):\n")
print(gsum %>% filter(lambda %in% c(1, 2, 3, 5, 10)) %>% select(lambda, rho_seqdef_m, rho_priority_m, jac10_seqdef_m, jac10_priority_m, top1_same_mean, target_top1_mean))
cat("\n[auto_max vs by_genus]\n"); for (v in names(meth)) cat(sprintf("  %-22s %s\n", v, q(meth[[v]])))
cat("\n[kernels vs exponential]\n")
for (k in unique(kern$kernel)) { kk <- filter(kern, kernel == k)
  cat(sprintf("  %-9s lambda %s | rho SeqDef %s | rho Priority %s | Jaccard top10 SeqDef %s | same top-1 %d/100\n",
              k, q(kk$lambda), q(kk$rho_seqdef), q(kk$rho_priority), q(kk$jac10_seqdef), sum(kk$top1_same))) }
cat(sprintf("  Brownian vs linear(lambda = 0.5): max |difference| over all trees = %.2e\n", max(kern$max_abs_diff_bm_vs_linear05)))
cat("\n[Priority composition]\n"); for (v in names(dom)) cat(sprintf("  %-18s %s\n", v, q(dom[[v]])))

summ <- bind_rows(
  tibble(comparison = "lambda grid (median across grid of per-lambda medians)", stat = names(grid)[-(1:2)],
         value = sapply(names(grid)[-(1:2)], function(v) median(gsum[[paste0(v, "_m")]]))),
  tibble(comparison = "auto_max vs by_genus", stat = names(meth), value = sapply(meth, median, na.rm = TRUE)),
  map_dfr(unique(kern$kernel), function(k) { kk <- filter(kern, kernel == k) %>% select(-kernel, -tree_id)
    tibble(comparison = paste("kernel", k, "vs exponential"), stat = names(kk), value = sapply(kk, median, na.rm = TRUE)) }),
  tibble(comparison = "Priority composition", stat = names(dom), value = sapply(dom, median)))
write.csv(summ, "results/stability_summary.csv", row.names = FALSE)

# ---- revised Fig S1 ------------------------------------------------------------
lam_auto <- map_dbl(posterior, "lambda_auto")
band <- list(annotate("rect", xmin = quantile(lam_auto, .025), xmax = quantile(lam_auto, .975), ymin = -Inf, ymax = Inf, fill = "#D73027", alpha = .12),
             geom_vline(xintercept = median(lam_auto), colour = "#D73027", linetype = "longdash"))
long <- function(stub) bind_rows(
  gsum %>% transmute(lambda, what = "SeqDef", m = .data[[paste0(stub[1], "_m")]], lo = .data[[paste0(stub[1], "_lo")]], hi = .data[[paste0(stub[1], "_hi")]]),
  gsum %>% transmute(lambda, what = "Priority", m = .data[[paste0(stub[2], "_m")]], lo = .data[[paste0(stub[2], "_lo")]], hi = .data[[paste0(stub[2], "_hi")]]))
cols <- c(SeqDef = "#1862C9", Priority = "#C96E18")
pA <- ggplot(long(c("rho_seqdef", "rho_priority")), aes(lambda, m, colour = what, fill = what)) + band +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) + geom_line(linewidth = .9) +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = expression(lambda), y = expression("Spearman " * rho * " with the auto_max solution"), tag = "A") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
pB <- ggplot(long(c("jac10_seqdef", "jac10_priority")), aes(lambda, m, colour = what, fill = what)) + band +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) + geom_line(linewidth = .9) +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = expression(lambda), y = "Top-10 agreement with auto_max (Jaccard; ties included)", tag = "B") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
pC <- ggplot(gsum, aes(lambda, target_top1_mean)) + band + geom_line(colour = "#D73027", linewidth = .9) + geom_point(colour = "#D73027", size = 1) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = expression(lambda), y = expression("Fraction of trees with " * italic("C. atromarginatus") * " top-ranked"), tag = "C") +
  theme_minimal(base_size = 10)
ggsave("figures/figS1_lambda_stability.pdf", patchwork::wrap_plots(pA, pB, pC, nrow = 1), width = 13, height = 4.3)
cat("\n[stability] DONE.\n")
