# =============================================================================
# brownian_comparison.R
# The Brownian-motion kernel versus the other kernels, on all 100 posterior
# trees (R2; R1 used the MCC tree only).
#   * identity check: on an ultrametric tree the Brownian correlation equals
#     the linear kernel at lambda = 0.5 (Reviewer 2, Major 3c);
#   * discrimination: SeqDef variance, Brownian vs exponential at auto_max;
#   * agreement with the exponential default (SeqDef and Priority) and the
#     top-ranked species.
# Produces: results/brownian_comparison.csv (one row per tree) and prints the
# values quoted in the manuscript.
# =============================================================================
source("analyses/00_setup.R")
posterior <- readRDS("results/posterior_priority.rds")
TOL <- 1e-6; top_set <- function(x) names(x)[x >= max(x) * (1 - TOL)]
sp <- function(a, b) cor(a, b, method = "spearman")
bc <- map_dfr(posterior, function(p) {
  inp <- build_input_binary(chond[[p$tree_id]]); nm <- inp$tree$tip.label; ge <- ge_for_tree(inp$tree)[nm]
  bm  <- suppressMessages(SeqDef(inp$tree, inp$df, kernel = "brownian"))
  l05 <- suppressMessages(SeqDef(inp$tree, inp$df, kernel = "linear", lambda = 0.5))
  e05 <- suppressMessages(SeqDef(inp$tree, inp$df, kernel = "exponential", lambda = 0.5))
  sd_exp <- p$seqdef_auto[nm]; pri_exp <- p$pri_auto[nm]; pri_bm <- bm$seqdef[nm] * 2^ge
  tibble(tree_id = p$tree_id,
         max_abs_diff_bm_linear05 = max(abs(bm$availability[nm] - l05$availability[nm])),
         rho_bm_exp05 = sp(bm$seqdef[nm], e05$seqdef[nm]),
         var_bm = var(bm$seqdef), var_exp_auto = var(sd_exp), lambda_auto = p$lambda_auto,
         rho_seqdef_bm_exp = sp(bm$seqdef[nm], sd_exp), rho_priority_bm_exp = sp(pri_bm, pri_exp),
         top1_bm = paste(sort(top_set(pri_bm)), collapse = ";"),
         target_top1_bm = identical(sort(top_set(pri_bm)), TARGET))
})
write.csv(bc, "results/brownian_comparison.csv", row.names = FALSE)
q <- function(v, d = 3) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), median(v), quantile(v, .025), quantile(v, .975))
cat(sprintf("\n[Brownian, 100 trees] max |A_BM - A_linear(0.5)| = %.2e\n", max(bc$max_abs_diff_bm_linear05)))
cat(sprintf("  rho(BM, exponential lambda=0.5) SeqDef: %s\n", q(bc$rho_bm_exp05)))
cat(sprintf("  SeqDef variance: BM %s vs exponential auto_max %s\n", q(bc$var_bm, 4), q(bc$var_exp_auto, 4)))
cat(sprintf("  rho vs exponential auto_max: SeqDef %s | Priority %s\n", q(bc$rho_seqdef_bm_exp), q(bc$rho_priority_bm_exp)))
cat(sprintf("  C. atromarginatus unique top-1 under BM: %d / 100\n", sum(bc$target_top1_bm)))
print(sort(table(sub("_.*", "", sapply(strsplit(bc$top1_bm, ";"), `[`, 1))), decreasing = TRUE))
