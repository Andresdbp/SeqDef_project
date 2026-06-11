# =============================================================================
# brownian_comparison.R
# Add the Brownian-motion kernel to the comparison and show empirically that it
# is the flat, low-discrimination limit of the exponential kernel.
#
# Produces: results/brownian_comparison.csv  (Brownian comparison; values reported in the manuscript text)
# =============================================================================
source("analyses/00_setup.R")
inp <- build_input_binary(tree_mcc); tr <- inp$tree; ge <- ge_for_tree(tr)

# Variance of SeqDef vs lambda (exponential) -- fast now via the O(n) traversal
lam_grid <- seq(0.1, 25, 0.1)
vv <- sapply(lam_grid, function(l) var(suppressMessages(SeqDef(tr, inp$df, lambda = l, kernel = "exponential"))$seqdef))

res_auto <- suppressMessages(SeqDef(tr, inp$df, lambda = "auto_max", kernel = "exponential"))
e05 <- suppressMessages(SeqDef(tr, inp$df, lambda = 0.5, kernel = "exponential"))
bm  <- suppressMessages(SeqDef(tr, inp$df, kernel = "brownian"))

common   <- names(bm$seqdef)
rho_flat <- cor(bm$seqdef[common], e05$seqdef[common], method = "spearman")
pri_bm <- calc_priority(bm, ge); pri_ex <- calc_priority(res_auto, ge)
rho_pri <- cor(pri_bm[common], pri_ex[common], method = "spearman")
top_bm <- names(sort(pri_bm, decreasing = TRUE))[1]; top_ex <- names(sort(pri_ex, decreasing = TRUE))[1]

cat(sprintf("\n[Brownian] rho(BM, exp lambda=0.5)=%.3f (flat-limit) | var(BM)=%.4f vs var(exp auto_max lambda=%.1f)=%.4f\n",
            rho_flat, var(bm$seqdef), res_auto$lambda, var(res_auto$seqdef)))
cat(sprintf("           rho(BM, exp auto_max) priority=%.3f | top: BM=%s, exp=%s\n", rho_pri, top_bm, top_ex))

write.csv(data.frame(
  metric = c("rho_BM_vs_exp0.5","var_BM","lambda_auto_max","var_exp_auto","rho_BM_vs_exp_auto_priority","top_BM","top_exp_auto"),
  value  = c(round(rho_flat,4), round(var(bm$seqdef),4), res_auto$lambda, round(var(res_auto$seqdef),4),
             round(rho_pri,4), top_bm, top_ex)), "results/brownian_comparison.csv", row.names = FALSE)

# --- (Figure removed) --------------------------------------------------------
# Brownian comparison values are reported inline in the manuscript text (data: results/brownian_comparison.csv).
cat("[brownian] DONE -> results/brownian_comparison.csv (values reported in the manuscript text).\n")
