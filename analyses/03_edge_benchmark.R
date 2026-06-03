# =============================================================================
# 03_edge_benchmark.R
# Reviewer 1 #3 / Editor M3:  Benchmark SeqDef-Priority against EDGE (and EDGE2).
#
# Frame: COMPLEMENTARITY. EDGE ignores existing genomic data; SeqDef conditions
# on it. The divergence (high-EDGE / low-SeqDef because a relative is already
# sequenced) is the value proposition.
#
# Everything is computed on the MCC tree so ED, EDGE, EDGE2 and SeqDef share an
# identical topology and branch lengths.
#
# Produces:
#   results/edge_benchmark.csv
#   results/edge_divergence.csv
#   figures/figS_edge_vs_seqdef.pdf
# =============================================================================

source("analyses/00_setup.R")
suppressPackageStartupMessages({ library(picante) })

inp <- build_input_binary(tree_mcc)
tr  <- inp$tree
ge  <- ge_for_tree(tr)
cats <- iucn_clean$iucn_category[match(tr$tip.label, iucn_clean$scientific_name)]
has_asm <- ncbi$assembly_availability[match(tr$tip.label, ncbi$scientific_name)]
has_asm[is.na(has_asm)] <- 0

# ---- SeqDef + Priority on the MCC tree --------------------------------------
res  <- suppressMessages(SeqDef(tr, inp$df, lambda = "auto_max"))
seqdef <- res$seqdef[tr$tip.label]
priority <- calc_priority(res, ge, "exponential", base = 2)[tr$tip.label]

# ---- Evolutionary Distinctiveness (fair proportion + equal splits) ----------
ed_fp <- evol.distinct(tr, type = "fair.proportion")
ed_es <- evol.distinct(tr, type = "equal.splits")
ED  <- setNames(ed_fp$w, ed_fp$Species)[tr$tip.label]
EDe <- setNames(ed_es$w, ed_es$Species)[tr$tip.label]

# ---- Classic EDGE (Isaac et al. 2007):  log(1+ED) + GE*log(2) ----------------
EDGE <- log(1 + ED) + ge * log(2)

# ---- EDGE2 (Gumbs et al. 2023), expected-PD-loss formulation -----------------
# IUCN -> 50-yr extinction probability (GE2). Values per Mooers et al. / Gumbs
# et al. 2023; CITE in manuscript.
pext_map <- c(LC = 0.0009, NT = 0.0071, VU = 0.0513, EN = 0.4276, CR = 0.9688)
pext <- pext_map[as.character(cats)]; pext[is.na(pext)] <- pext_map["LC"]
names(pext) <- tr$tip.label

# Expected PD loss attributed to each tip: for each edge, prob(edge lost) =
# product of pext over descendant tips; expected loss = L * that, split equally
# among descendants. ED2 sums a tip's attributed expected loss along its path.
desc <- phangorn::Descendants(tr, seq_len(max(tr$edge)), type = "tips")
ED2 <- setNames(numeric(length(tr$tip.label)), tr$tip.label)
for (e in seq_len(nrow(tr$edge))) {
  child <- tr$edge[e, 2]; L <- tr$edge.length[e]
  dt <- desc[[child]]
  w  <- prod(pext[dt])
  ED2[dt] <- ED2[dt] + L * w / length(dt)
}
EDGE2 <- ED2   # already incorporates extinction probability (expected loss)

# ---- Assemble & rank --------------------------------------------------------
rnk <- function(x) rank(-x, ties.method = "min")
bench <- tibble(
  species = tr$tip.label, iucn = as.character(cats), GE = ge, has_assembly = has_asm,
  ED = ED, ED_equal = EDe, EDGE = EDGE, EDGE2 = EDGE2,
  SeqDef = seqdef, Priority = priority,
  rank_EDGE = rnk(EDGE), rank_EDGE2 = rnk(EDGE2),
  rank_SeqDef = rnk(seqdef), rank_Priority = rnk(priority)
) %>% arrange(rank_Priority)
write.csv(bench, "results/edge_benchmark.csv", row.names = FALSE)

sp_cor <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
top_k <- function(x, k = 25) bench$species[order(-x)][1:k]
ov <- function(a, b, k = 25) length(intersect(top_k(a, k), top_k(b, k)))

cat("\n[EDGE benchmark on MCC tree] n =", nrow(bench), "\n")
cat(sprintf("  Spearman rho(EDGE,  Priority) = %.3f\n", sp_cor(EDGE, priority)))
cat(sprintf("  Spearman rho(EDGE2, Priority) = %.3f\n", sp_cor(EDGE2, priority)))
cat(sprintf("  Spearman rho(EDGE,  SeqDef)   = %.3f\n", sp_cor(EDGE, seqdef)))
cat(sprintf("  Spearman rho(ED,    SeqDef)   = %.3f\n", sp_cor(ED, seqdef)))
cat(sprintf("  top-25 overlap EDGE vs Priority  = %d/25\n", ov(EDGE, priority)))
cat(sprintf("  top-25 overlap EDGE2 vs Priority = %d/25\n", ov(EDGE2, priority)))
tw <- bench[bench$species == TARGET, ]
cat(sprintf("  %s: Priority rank %d, EDGE rank %d, EDGE2 rank %d, SeqDef rank %d (IUCN %s)\n",
            TARGET, tw$rank_Priority, tw$rank_EDGE, tw$rank_EDGE2, tw$rank_SeqDef, tw$iucn))

# ---- Divergence: high-EDGE but low-SeqDef BECAUSE a relative is sequenced ----
# NOTE: derive genus/congener flags from bench's OWN (already-arranged) columns
# so they stay row-aligned with species.
bench$genus <- sub("_.*", "", bench$species)
genus_has_asm <- tapply(bench$has_assembly, bench$genus, max)  # 1 if any congener has an assembly
bench$congener_sequenced <- genus_has_asm[bench$genus] == 1 & bench$has_assembly == 0

# Rank the divergence directly: large positive gap = high EDGE rank but low SeqDef
# rank, conditioned on a congener already being sequenced (the mechanism).
bench$gap <- bench$rank_SeqDef - bench$rank_EDGE
divergence <- bench %>%
  filter(congener_sequenced, rank_EDGE <= 100) %>%
  arrange(desc(gap)) %>%
  select(species, iucn, ED, rank_EDGE, SeqDef, rank_SeqDef, rank_Priority, gap, genus)
write.csv(divergence, "results/edge_divergence.csv", row.names = FALSE)
exemplars <- head(divergence$species, 5)   # clearest cases (largest rank gap)
# label two well-separated exemplars (one Carcharhinus, one Sphyrna) to avoid overlap
label_sp  <- c(divergence$species[1], divergence$species[grep("^Sphyrna", divergence$species)[1]])
cat(sprintf("\n[Divergence] top high-EDGE / low-SeqDef species with a sequenced congener (by rank gap):\n"))
print(head(divergence, 6))

# ---- Figure: EDGE vs SeqDef (the orthogonality / complementarity is the point)
abbr <- function(s) sapply(strsplit(s, "_"), function(p) paste0(substr(p[1], 1, 1), ". ", p[2]))
bench$flag <- ifelse(bench$species == TARGET, "target",
                     ifelse(bench$species %in% exemplars, "divergent", "other"))
pdat <- bench %>% mutate(EDGE_s = scales::rescale(EDGE))
lab_tgt <- subset(pdat, species == TARGET)
lab_div <- subset(pdat, species %in% label_sp); lab_div$lab <- abbr(lab_div$species)
p <- ggplot(pdat, aes(EDGE_s, SeqDef)) +
  geom_point(data = subset(pdat, flag == "other"), color = "grey75", size = 0.8, alpha = 0.5) +
  geom_point(data = subset(pdat, flag == "divergent"), color = "#762a83", size = 2) +
  geom_point(data = lab_tgt, color = "#D73027", size = 3) +
  geom_text(data = lab_tgt, label = "C. atromarginatus", color = "#D73027", vjust = -1, hjust = 0.9, size = 3, fontface = "italic") +
  geom_text(data = lab_div, aes(label = lab), color = "#762a83", vjust = -0.9, size = 2.8, fontface = "italic") +
  labs(x = "EDGE (scaled)", y = "SeqDef (genomic deficiency)",
       subtitle = sprintf("EDGE and SeqDef are near-orthogonal (rho = %.2f). Purple: the clearest high-EDGE / low-SeqDef\ncases (a congener is already sequenced); C. atromarginatus (red) is high on both.",
                          sp_cor(EDGE, seqdef))) +
  theme_minimal() + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9))
ggsave("figures/figS_edge_vs_seqdef.pdf", p, width = 6.5, height = 5)

cat("\n[03] DONE. Wrote results/edge_benchmark.csv, results/edge_divergence.csv, figures/figS_edge_vs_seqdef.pdf\n")
