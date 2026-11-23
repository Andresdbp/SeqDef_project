# =============================================================================
# fig1_toy_tree.R
# R2 Reviewer 2, Minor 1 / Editor: Figure 1 now shows REAL SeqDef values
# computed by the package on the toy tree at three lambda values (plus the
# auto_max choice), instead of illustrative evenly spaced numbers.
# Produces: figures/fig1.pdf, figures/png/fig1.png, results/fig1_values.csv
# =============================================================================
suppressPackageStartupMessages(library(ape))
source("function.R")
tr <- read.tree(text = "(taxa1:100,(taxa2:90,(taxa3:80,(taxa4:70,(taxa5:60,(taxa6:50,(taxa7:40,(taxa8:30,(taxa9:10,taxa10:10):20):10):10):10):10):10):10):10);")
df <- data.frame(label = paste0("taxa", 1:10), S = c(rep(0, 9), 1))
lams <- c(1, 3, 10)
vals <- sapply(lams, function(l) suppressMessages(SeqDef(tr, df, lambda = l))$seqdef[tr$tip.label])
auto <- suppressMessages(SeqDef(tr, df, lambda = "auto_max"))
out <- data.frame(tip = tr$tip.label, S = df$S[match(tr$tip.label, df$label)], round(vals, 2))
names(out)[3:5] <- paste0("lambda_", lams)
out$auto_max <- round(auto$seqdef[tr$tip.label], 2)
write.csv(out, "results/fig1_values.csv", row.names = FALSE)
cat("auto_max lambda on toy tree:", auto$lambda, "\n"); print(out)

draw <- function() {
  par(mar = c(1, 1, 1, 1))
  plot(tr, x.lim = 230, y.lim = c(0.5, 12), align.tip.label = TRUE, font = 3, label.offset = 2, edge.width = 2)
  y <- seq_along(tr$tip.label); o <- match(tr$tip.label, out$tip)
  xs <- c(138, 168, 192, 216)
  text(xs[1], y, out$S[o]); text(xs[1], 10.9, "Sequence\navailability (S)", font = 2, cex = .75, adj = c(.5, 0))
  for (k in 1:3) { text(xs[k + 1], y, sprintf("%.2f", vals[o, k]))
    text(xs[k + 1], 10.9, bquote(atop(bold(SeqDef), lambda == .(lams[k]))), cex = .75, adj = c(.5, 0)) }
}
pdf("figures/fig1.pdf", width = 7.5, height = 6); draw(); dev.off()
dir.create("figures/png", showWarnings = FALSE)
png("figures/png/fig1.png", width = 7.5, height = 6, units = "in", res = 300); draw(); dev.off()
