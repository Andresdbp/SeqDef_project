# =============================================================================
# 05_fig1_regen.R
# Reviewer 2 + Editor:  Figure 1 must report the lambda used. Rebuilt from an
# actual SeqDef() run on the 10-taxon toy tree. Three SeqDef columns show the
# linear kernel, the exponential auto_max lambda, and exponential lambda = 5.
# Title/notes go in the caption (manuscript_changes.md), not on the figure.
#
# Produces:  figures/fig1.pdf  (overwrites);  results/fig1_values.csv
# =============================================================================

source("analyses/00_setup.R")

# Toy tree + availability (taxa10 = the only sequenced tip; taxa9 is its sister)
tree_text <- "(taxa1:100,(taxa2:90,(taxa3:80,(taxa4:70,(taxa5:60,(taxa6:50,(taxa7:40,(taxa8:30,(taxa9:10,taxa10:10):20):10):10):10):10):10):10):10);"
tr  <- read.tree(text = tree_text)
dat <- read.csv("seqdef.csv", header = FALSE, col.names = c("label", "Availability"))
dat <- dat[match(tr$tip.label, dat$label), ]

# Three SeqDef columns: (1) linear kernel, (2) exponential auto_max, (3) exponential lambda = 5
fit_lin  <- suppressMessages(SeqDef(tr, dat, data.col = "Availability", lambda = "auto_max", kernel = "linear"))
fit_auto <- suppressMessages(SeqDef(tr, dat, data.col = "Availability", lambda = "auto_max", kernel = "exponential"))
fit_l5   <- suppressMessages(SeqDef(tr, dat, data.col = "Availability", lambda = 5,          kernel = "exponential"))
lam_auto <- round(fit_auto$lambda, 1)
lam_lin  <- round(fit_lin$lambda, 1)

sc <- list(fit_lin$seqdef[tr$tip.label], fit_auto$seqdef[tr$tip.label], fit_l5$seqdef[tr$tip.label])
h2 <- list("linear", bquote(lambda == .(lam_auto)), bquote(lambda == 5))

out <- data.frame(taxon = tr$tip.label, sequencing_availability = dat$Availability,
                  seqdef_linear = sc[[1]], seqdef_exp_auto = sc[[2]], seqdef_exp_l5 = sc[[3]])
write.csv(out, "results/fig1_values.csv", row.names = FALSE)
cat("\n[Fig 1] exponential auto_max lambda =", lam_auto, "| linear auto_max lambda =", lam_lin, "\n")
print(out, row.names = FALSE)

# ---- Plot: tree + 'Sequencing Availability' column + 3 SeqDef columns -------
pdf("figures/fig1.pdf", width = 9, height = 5)
op <- par(mar = c(1, 1, 2.5, 1))
plot(tr, x.lim = 240, y.lim = c(0.8, 11.3), align.tip.label = TRUE,
     font = 3, label.offset = 2, edge.width = 2, cex = 0.95)
pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
yy <- pp$yy[seq_len(length(tr$tip.label))]      # tip y-coords in tip-index order

cols_x <- c(130, 168, 200, 230)
h1 <- max(yy) + 1.0    # header line 1
h2y <- max(yy) + 0.5   # header line 2

# data columns
text(cols_x[1], yy, labels = out$sequencing_availability, cex = 0.95, xpd = NA)
for (j in 1:3) text(cols_x[j + 1], yy, labels = sprintf("%.2f", sc[[j]]), cex = 0.95, xpd = NA)

# headers (two lines each, aligned)
text(cols_x[1], h1,  "Sequencing",   font = 2, cex = 0.78, xpd = NA)
text(cols_x[1], h2y, "Availability", font = 2, cex = 0.78, xpd = NA)
for (j in 1:3) {
  text(cols_x[j + 1], h1,  "SeqDef", font = 2, cex = 0.82, xpd = NA)
  text(cols_x[j + 1], h2y, h2[[j]], font = 2, cex = 0.82, xpd = NA)
}
par(op)
dev.off()

cat("\n[05] DONE. Overwrote figures/fig1.pdf; wrote results/fig1_values.csv\n")
