# =============================================================================
# kernel_comparison.R
# Reviewer 1 #2b / Editor M2b:  Justify the exponential kernel by showing that
# SeqDef rankings are robust to the kernel choice.
#
# Compares all four kernels (exponential / gaussian / linear / brownian), each
# with its own auto_max lambda where applicable, on the MCC tree + all 100
# posterior trees.  Values are reported in the manuscript text.
#
# Produces:
#   results/kernel_comparison.csv
# =============================================================================

source("analyses/00_setup.R")
kernels <- c("exponential", "gaussian", "linear", "brownian")
top_set <- function(pri) names(pri)[pri == max(pri, na.rm = TRUE)]

run_kernels <- function(phy_raw, label) {
  inp <- build_input_binary(phy_raw)
  if (is.null(inp)) return(NULL)
  ge <- ge_for_tree(inp$tree)
  res <- lapply(kernels, function(k) {
    r <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = "auto_max", kernel = k))
    list(pri = calc_priority(r, ge, "exponential", base = 2), lambda = r$lambda)
  })
  names(res) <- kernels
  common <- names(res$exponential$pri)
  rho <- function(k) suppressWarnings(cor(res$exponential$pri[common], res[[k]]$pri[common], method = "spearman"))
  tibble(
    tree = label,
    lambda_exp = res$exponential$lambda, lambda_gauss = res$gaussian$lambda, lambda_lin = res$linear$lambda,
    rho_gauss = rho("gaussian"), rho_lin = rho("linear"), rho_bm = rho("brownian"),
    target_exp   = TARGET %in% top_set(res$exponential$pri),
    target_gauss = TARGET %in% top_set(res$gaussian$pri),
    target_lin   = TARGET %in% top_set(res$linear$pri),
    target_bm    = TARGET %in% top_set(res$brownian$pri)
  )
}

message("02: kernel comparison on MCC + all 100 posterior trees...")
kc <- bind_rows(
  run_kernels(tree_mcc, "MCC"),
  map_dfr(seq_along(chond), ~ run_kernels(chond[[.x]], paste0("post_", .x)), .progress = TRUE)
)
write.csv(kc, "results/kernel_comparison.csv", row.names = FALSE)

post <- kc %>% filter(tree != "MCC")
qci <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .025), quantile(x, .975))
cat(sprintf("\n[Kernel comparison] (n = %d posterior trees)\n", nrow(post)))
cat(sprintf("  auto_max lambda (median): exp %.2f | gaussian %.2f | linear %.2f  (brownian: parameter-free)\n",
            median(post$lambda_exp), median(post$lambda_gauss), median(post$lambda_lin)))
cat(sprintf("  rho(Priority vs exponential): gaussian %s | linear %s | brownian %s\n",
            qci(post$rho_gauss), qci(post$rho_lin), qci(post$rho_bm)))
cat(sprintf("  C. atromarginatus top target (of %d): exp %d | gaussian %d | linear %d | brownian %d\n",
            nrow(post), sum(post$target_exp), sum(post$target_gauss), sum(post$target_lin), sum(post$target_bm)))
cat("\n[kernel] DONE. Wrote results/kernel_comparison.csv (values reported in the manuscript text).\n")
