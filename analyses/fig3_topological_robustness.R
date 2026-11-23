# =============================================================================
# fig3_topological_robustness.R
# R2 Reviewer 2, Major 1 (Fig. 3 sums to 107): the R1 figure counted every
# species tied for the top Priority in a tree, and tied species were separated
# only by floating-point noise (~1e-9) because tips in a clade with no
# sequenced member receive identical SeqDef (and identical Priority if they
# share an IUCN category).  Here ties are detected with a relative tolerance
# (1e-6) and each tree contributes exactly one unit, assigned to its winning
# set.  Produces figures/fig3.pdf, figures/png/fig3.png, results/fig3_top1.csv
# =============================================================================
source("analyses/00_setup.R")
posterior <- readRDS("results/posterior_priority.rds")
TOL <- 1e-6
win <- map(posterior, ~ { x <- .x$pri_auto; sort(names(x)[x >= max(x) * (1 - TOL)]) })
fmt <- function(w) if (length(w) == 1) sub("_", " ", w) else {
  g <- unique(sub("_.*", "", w)); sprintf("%d-way tie: CR %s", length(w), paste(g, collapse = " + ")) }
tab <- tibble(tree_id = map_int(posterior, "tree_id"), n_tied = lengths(win), winners = map_chr(win, paste, collapse = ";"),
              label = map_chr(win, fmt))
write.csv(tab, "results/fig3_top1.csv", row.names = FALSE)
cnt <- tab %>% count(label, sort = TRUE)
cat(sprintf("trees: %d | trees with a tie for top-1: %d (sizes: %s)\n", nrow(tab), sum(tab$n_tied > 1),
            paste(names(table(tab$n_tied[tab$n_tied > 1])), "x", table(tab$n_tied[tab$n_tied > 1]), collapse = ", ")))
print(cnt)
p <- ggplot(cnt, aes(n, reorder(label, n))) + geom_col(fill = c("#D73027", rep("grey60", nrow(cnt) - 1))) +
  geom_text(aes(label = n), hjust = -0.3, size = 3.4) + scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .02))) +
  labs(x = "Number of posterior trees (of 100) in which the set ranks first", y = NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "italic"))
ggsave("figures/fig3.pdf", p, width = 7.5, height = 2.6)
ggsave("figures/png/fig3.png", p, width = 7.5, height = 2.6, dpi = 300)
