# =============================================================================
# 06_continuous_S.R
# Reviewer 1 #1 / Reviewer 2:  Demonstrate that S need not be binary. We build an
# illustrative CONTINUOUS availability score from the NCBI assembly COUNT via a
# saturating map S = count / (count + 1) (diminishing returns; no API call), and
# compare the resulting prioritization with the binary-S case study.
#
# Produces:
#   results/continuous_S_example.csv
#   figures/figS_continuousS.pdf
# =============================================================================

source("analyses/00_setup.R")
top_set <- function(pri) names(pri)[pri == max(pri, na.rm = TRUE)]

inp <- build_input_binary(tree_mcc)
tr  <- inp$tree
ge  <- ge_for_tree(tr)

count <- ncbi$assembly_count[match(tr$tip.label, ncbi$scientific_name)]
count[is.na(count)] <- 0
S_bin  <- as.integer(count >= 1)
S_cont <- count / (count + 1)            # 0 -> 0, 1 -> 0.5, 2 -> 0.67, 5 -> 0.83

cat("[Continuous S] assembly_count distribution among sequenced spp:\n")
print(table(count[count > 0]))

df_bin  <- tibble(taxa = tr$tip.label, score = S_bin)
df_cont <- tibble(taxa = tr$tip.label, score = S_cont)

r_bin  <- suppressMessages(SeqDef(tr, df_bin,  lambda = "auto_max"))
r_cont <- suppressMessages(SeqDef(tr, df_cont, lambda = "auto_max"))
p_bin  <- calc_priority(r_bin,  ge, "exponential", base = 2)[tr$tip.label]
p_cont <- calc_priority(r_cont, ge, "exponential", base = 2)[tr$tip.label]

rho_sd  <- suppressWarnings(cor(r_bin$seqdef[tr$tip.label], r_cont$seqdef[tr$tip.label], method = "spearman"))
rho_pri <- suppressWarnings(cor(p_bin, p_cont, method = "spearman"))

rnk <- function(x) rank(-x, ties.method = "min")
res <- tibble(
  species = tr$tip.label, assembly_count = count, S_binary = S_bin, S_cont = round(S_cont, 3),
  seqdef_bin = r_bin$seqdef[tr$tip.label], seqdef_cont = r_cont$seqdef[tr$tip.label],
  priority_bin = p_bin, priority_cont = p_cont,
  rank_bin = rnk(p_bin), rank_cont = rnk(p_cont)
) %>% mutate(rank_shift = rank_cont - rank_bin) %>% arrange(rank_bin)
write.csv(res, "results/continuous_S_example.csv", row.names = FALSE)

cat(sprintf("\n[Continuous vs binary S, MCC tree]\n  lambda: binary=%.2f  continuous=%.2f\n", r_bin$lambda, r_cont$lambda))
cat(sprintf("  Spearman rho(SeqDef)   binary vs continuous = %.3f\n", rho_sd))
cat(sprintf("  Spearman rho(Priority) binary vs continuous = %.3f\n", rho_pri))
cat(sprintf("  top-10 priority overlap = %d/10\n",
            length(intersect(res$species[order(-p_bin)][1:10], res$species[order(-p_cont)][1:10]))))
cat(sprintf("  %s priority rank: binary %d -> continuous %d\n",
            TARGET, res$rank_bin[res$species == TARGET], res$rank_cont[res$species == TARGET]))
cat("\n  Largest priority increases under continuous S (relatives only partially covered):\n")
print(head(res %>% filter(assembly_count == 0) %>% arrange(rank_shift) %>%
             select(species, rank_bin, rank_cont, rank_shift), 5), row.names = FALSE)

# ---- Figure -----------------------------------------------------------------
res$flag <- ifelse(res$species == TARGET, "target", "other")
p1 <- ggplot(res, aes(seqdef_bin, seqdef_cont)) +
  geom_abline(slope = 1, intercept = 0, color = "grey80") +
  geom_point(aes(color = flag, size = flag), alpha = 0.6) +
  scale_color_manual(values = c(other = "grey55", target = "#D73027"), guide = "none") +
  scale_size_manual(values = c(other = 0.8, target = 2.8), guide = "none") +
  labs(x = "SeqDef (binary S)", y = "SeqDef (continuous S)",
       subtitle = sprintf("A  SeqDef: binary vs continuous  (rho = %.3f)", rho_sd)) +
  theme_minimal() + theme(panel.grid.minor = element_blank())

p2 <- ggplot(res, aes(priority_bin, priority_cont)) +
  geom_abline(slope = 1, intercept = 0, color = "grey80") +
  geom_point(aes(color = flag, size = flag), alpha = 0.6) +
  scale_color_manual(values = c(other = "grey55", target = "#D73027"), guide = "none") +
  scale_size_manual(values = c(other = 0.8, target = 2.8), guide = "none") +
  labs(x = "Priority (binary S)", y = "Priority (continuous S)",
       subtitle = sprintf("B  Priority: binary vs continuous  (rho = %.3f)", rho_pri)) +
  theme_minimal() + theme(panel.grid.minor = element_blank())

ggsave("figures/figS_continuousS.pdf", patchwork::wrap_plots(p1, p2, nrow = 1), width = 10, height = 4)
cat("\n[06] DONE. Wrote results/continuous_S_example.csv, figures/figS_continuousS.pdf\n")
