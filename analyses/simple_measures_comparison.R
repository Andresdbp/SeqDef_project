# =============================================================================
# simple_measures_comparison.R
# R2 (round 2) Reviewer 2, Major 1:  Does SeqDef add information over simpler
# measures?  Compares SeqDef (exponential kernel, auto_max lambda) with
#   (a) NSD  = phylogenetic distance from each species to its nearest sequenced
#              relative (depth-normalised, same scale as SeqDef's kernel), and
#   (b) FAM  = whether the species' family already has >= 1 genome assembly.
# Also verifies the clade-constancy property raised by the reviewer and
# characterises NSD as the lambda -> infinity limit of SeqDef.
#
# Produces:
#   results/simple_measures_by_tree.csv     (per posterior tree statistics)
#   results/simple_measures_mcc.csv         (per species, MCC tree)
#   results/simple_measures_lambda.csv      (rho(SeqDef, NSD) vs lambda)
#   results/simple_measures_discordant.csv  (clade pairs ranked differently)
#   figures/figS4_simple_measures.pdf
# =============================================================================
source("analyses/00_setup.R")

# ---- family lookup from the IUCN taxon records ------------------------------
fam_tbl <- map_dfr(a_data, function(el) tibble(
  scientific_name = gsub(" ", "_", str_squish(el$taxon$scientific_name[1])),
  family = tryCatch(str_to_title(el$taxon$family_name[1]), error = function(e) NA_character_)))
fam_of <- function(sp) fam_tbl$family[match(sp, fam_tbl$scientific_name)]

posterior <- readRDS("results/posterior_priority.rds")
top_set   <- function(x) names(x)[x >= max(x, na.rm = TRUE) * (1 - 1e-6)]   # relative tolerance: tied tips differ only by rounding
sp  <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))

# Maximal clades with no sequenced member: tips inside get identical SeqDef on
# an ultrametric tree (reviewer's point).  Returns a tip -> clade-id vector.
unseq_clades <- function(tr, s) {
  n <- Ntip(tr); desc <- phangorn::Descendants(tr, seq_len(n + tr$Nnode), "tips")
  allzero <- vapply(desc, function(d) all(s[d] == 0), logical(1))
  par <- setNames(tr$edge[, 1], tr$edge[, 2])
  id <- rep(NA_integer_, n)
  for (t in which(s == 0)) {           # climb while parent clade is still all-unsequenced
    v <- t
    while (!is.na(par[as.character(v)]) && allzero[par[as.character(v)]]) v <- par[as.character(v)]
    id[t] <- v
  }
  id
}

per_tree <- function(phy, sd_given = NULL) {
  inp <- build_input_binary(phy); tr <- inp$tree; s <- inp$df$score
  td  <- max(ape::node.depth.edgelength(tr))
  D   <- cophenetic.phylo(tr)[tr$tip.label, tr$tip.label] / td
  Dseq <- D[, s == 1, drop = FALSE]; Dseq[cbind(which(s == 1), seq_len(sum(s)))] <- Inf
  nsd <- apply(Dseq, 1, min); nsd[s == 1] <- 0           # sequenced tips: 0
  sd  <- if (is.null(sd_given)) suppressMessages(SeqDef(tr, inp$df, lambda = "auto_max"))$seqdef else sd_given
  sd  <- sd[tr$tip.label]
  fam <- fam_of(tr$tip.label)
  fam_has <- tapply(s, fam, max)[fam] == 1
  ge  <- ge_for_tree(tr)
  list(tr = tr, s = s, sd = sd, nsd = nsd, fam = fam, fam_has = fam_has, ge = ge,
       clade = unseq_clades(tr, s))
}

# =============================================================================
# A. Across the 100 posterior trees
# =============================================================================
message("A: SeqDef vs NSD / family-gap on 100 posterior trees ...")
by_tree <- map_dfr(posterior, function(p) {
  x <- per_tree(chond[[p$tree_id]], p$seqdef_auto)
  u <- x$s == 0
  pri <- x$sd * 2^x$ge
  pri_nsd <- (x$nsd / max(x$nsd)) * 2^x$ge
  # clade-level (one row per maximal unsequenced clade) correlation
  cl <- tibble(clade = x$clade[u], sd = x$sd[u], nsd = x$nsd[u]) %>%
    group_by(clade) %>% summarise(sd_range = diff(range(sd)), nsd_range = diff(range(nsd)),
                                  sd = mean(sd), nsd = mean(nsd), .groups = "drop")
  topk <- function(v, k) names(sort(v, decreasing = TRUE))[seq_len(k)]
  first_fam <- function(k) mean(!x$fam_has[match(topk(pri, k), names(x$sd))])
  t1 <- top_set(pri)
  tibble(
    tree_id = p$tree_id, lambda = p$lambda_auto,
    n_unseq = sum(u), n_unseq_clades = nrow(cl),
    n_distinct_sd_unseq = n_distinct(round(x$sd[u], 8)),
    max_within_clade_sd_range = max(cl$sd_range),
    rho_sd_nsd_tips   = sp(x$sd[u], x$nsd[u]),
    rho_sd_nsd_clades = sp(cl$sd, cl$nsd),
    rho_sd_fam  = sp(x$sd[u], as.numeric(!x$fam_has[u])),
    rho_pri_prinsd = sp(pri, pri_nsd),
    top1_n_tied = length(t1), top1_is_target = TARGET %in% t1,
    top1_first_in_family = mean(!x$fam_has[match(t1, names(x$sd))]),
    top10_first_fam = first_fam(10), top25_first_fam = first_fam(25), top50_first_fam = first_fam(50),
    top25_overlap_nsd = length(intersect(topk(pri, 25), topk(pri_nsd, 25))),
    top1_same_as_nsd = length(intersect(t1, top_set(pri_nsd))) > 0,
    frac_unseq_in_seq_family = mean(x$fam_has[u])
  )
}, .progress = TRUE)
write.csv(by_tree, "results/simple_measures_by_tree.csv", row.names = FALSE)

q <- function(v, d = 2) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), median(v), quantile(v, .025), quantile(v, .975))
cat("\n================ SIMPLE-MEASURES COMPARISON (100 posterior trees) ================\n")
cat(sprintf("unsequenced tips: %d in %s maximal unsequenced clades\n", by_tree$n_unseq[1], q(by_tree$n_unseq_clades, 0)))
cat(sprintf("distinct SeqDef values among unsequenced tips: %s\n", q(by_tree$n_distinct_sd_unseq, 0)))
cat(sprintf("max within-clade SeqDef range (should be ~0): %.2e\n", max(by_tree$max_within_clade_sd_range)))
cat(sprintf("rho(SeqDef, NSD) unsequenced tips  : %s\n", q(by_tree$rho_sd_nsd_tips, 3)))
cat(sprintf("rho(SeqDef, NSD) clade level       : %s\n", q(by_tree$rho_sd_nsd_clades, 3)))
cat(sprintf("rho(SeqDef, family lacks assembly) : %s\n", q(by_tree$rho_sd_fam, 3)))
cat(sprintf("rho(Priority, NSD x 2^GE)          : %s\n", q(by_tree$rho_pri_prinsd, 3)))
cat(sprintf("top-1 Priority tied (>1 species)   : %d/100 trees; total winners counted = %d\n",
            sum(by_tree$top1_n_tied > 1), sum(by_tree$top1_n_tied)))
cat(sprintf("top-1 would be first family assembly: %.0f%% of trees\n", 100 * mean(by_tree$top1_first_in_family)))
cat(sprintf("top-10/25/50 first-in-family share : %s / %s / %s\n", q(by_tree$top10_first_fam), q(by_tree$top25_first_fam), q(by_tree$top50_first_fam)))
cat(sprintf("top-25 overlap with NSD-priority   : %s\n", q(by_tree$top25_overlap_nsd, 0)))
cat(sprintf("same top-1 as NSD-priority         : %d/100\n", sum(by_tree$top1_same_as_nsd)))
cat(sprintf("unsequenced spp. whose family has an assembly: %s\n", q(by_tree$frac_unseq_in_seq_family)))

# =============================================================================
# B. MCC tree: per-species table and discordant cases
# =============================================================================
x <- per_tree(tree_mcc)
mcc <- tibble(species = names(x$sd), family = x$fam, iucn = names(risk_index)[match(x$ge, risk_index)],
              S = x$s, SeqDef = x$sd, NSD = x$nsd, family_has_assembly = x$fam_has, clade = x$clade,
              Priority = x$sd * 2^x$ge) %>%
  mutate(rank_SeqDef = rank(-SeqDef, ties.method = "min"), rank_NSD = rank(-NSD, ties.method = "min"),
         rank_Priority = rank(-Priority, ties.method = "min"))
write.csv(mcc, "results/simple_measures_mcc.csv", row.names = FALSE)

# Clade-level discordance: pairs of unsequenced clades ordered differently by
# SeqDef and NSD.  Explain via the number of sequenced tips within a horizon.
tr <- x$tr; td <- max(node.depth.edgelength(tr))
D <- cophenetic.phylo(tr)[tr$tip.label, tr$tip.label] / td
cl <- mcc %>% filter(S == 0) %>% group_by(clade) %>%
  summarise(n_tips = n(), example = first(species), families = paste(unique(family), collapse = "/"),
            SeqDef = mean(SeqDef), NSD = mean(NSD), .groups = "drop") %>%
  mutate(n_seq_within_2NSD = map_int(example, ~ sum(D[.x, x$s == 1] <= 2 * NSD[example == .x] + 1e-12)),
         n_seq_within_0.5 = map_int(example, ~ sum(D[.x, x$s == 1] <= 0.5)))
pairs <- expand_grid(a = seq_len(nrow(cl)), b = seq_len(nrow(cl))) %>% filter(a < b) %>%
  mutate(dNSD = cl$NSD[a] - cl$NSD[b], dSD = cl$SeqDef[a] - cl$SeqDef[b]) %>%
  filter(sign(dNSD) * sign(dSD) < 0)
cat(sprintf("\n[MCC] %d unsequenced clades; %d of %d clade pairs (%.1f%%) are ordered differently by SeqDef and NSD\n",
            nrow(cl), nrow(pairs), choose(nrow(cl), 2), 100 * nrow(pairs) / choose(nrow(cl), 2)))
# most informative discordances: large NSD advantage but lower SeqDef
disc <- pairs %>% mutate(score = abs(dNSD) * abs(dSD)) %>% arrange(desc(score)) %>% head(15) %>%
  mutate(clade_A = cl$example[a], fam_A = cl$families[a], NSD_A = cl$NSD[a], SD_A = cl$SeqDef[a], nseq_A = cl$n_seq_within_0.5[a],
         clade_B = cl$example[b], fam_B = cl$families[b], NSD_B = cl$NSD[b], SD_B = cl$SeqDef[b], nseq_B = cl$n_seq_within_0.5[b]) %>%
  select(clade_A, fam_A, NSD_A, SD_A, nseq_A, clade_B, fam_B, NSD_B, SD_B, nseq_B)
write.csv(disc, "results/simple_measures_discordant.csv", row.names = FALSE)
print(disc, width = 200)

# family-measure discordances on MCC
cat("\n[MCC] high-SeqDef species in families that DO have an assembly (top 8):\n")
print(mcc %>% filter(S == 0, family_has_assembly) %>% arrange(desc(SeqDef)) %>% select(species, family, SeqDef, NSD, iucn) %>% head(8))
cat("\n[MCC] low-SeqDef species in families with NO assembly (bottom 8):\n")
print(mcc %>% filter(S == 0, !family_has_assembly) %>% arrange(SeqDef) %>% select(species, family, SeqDef, NSD, iucn) %>% head(8))
cat("\nfamilies with no assembly:", n_distinct(mcc$family[!mcc$family_has_assembly]), "of", n_distinct(mcc$family), "\n")

# =============================================================================
# C. NSD is the lambda -> infinity limit: rho(SeqDef(lambda), NSD) vs lambda
# =============================================================================
lam_grid <- c(0.25, 0.5, 1, 2, 3.4, 5, 10, 25, 50)   # beyond ~50 exp(-lambda x) underflows for distant tips and ties appear
lam_tab <- map_dfr(posterior, function(p) {   # all 100 posterior trees
  inp <- build_input_binary(chond[[p$tree_id]]); xx <- per_tree(chond[[p$tree_id]], p$seqdef_auto); u <- xx$s == 0
  map_dfr(lam_grid, function(l) {
    r <- suppressMessages(SeqDef(inp$tree, inp$df, lambda = l))$seqdef[names(xx$sd)]
    tibble(tree_id = p$tree_id, lambda = l, rho = sp(r[u], xx$nsd[u]))
  })
})
write.csv(lam_tab, "results/simple_measures_lambda.csv", row.names = FALSE)
cat("\n[lambda limit] rho(SeqDef(lambda), NSD) among unsequenced tips (100 trees):\n")
print(lam_tab %>% group_by(lambda) %>% summarise(rho_med = median(rho), lo = quantile(rho, .025), hi = quantile(rho, .975)))

# =============================================================================
# D. Figure S4
# =============================================================================
pA <- ggplot(mcc %>% filter(S == 0), aes(NSD, SeqDef, colour = family_has_assembly)) +
  geom_point(alpha = .6, size = 1.3) +
  geom_point(data = subset(mcc, species == TARGET), colour = "#D73027", size = 3) +
  scale_colour_manual(values = c(`TRUE` = "grey55", `FALSE` = "#1862C9"),
                      labels = c(`TRUE` = "family has an assembly", `FALSE` = "family has no assembly"), name = NULL) +
  labs(x = "Distance to nearest sequenced relative (d / tree depth)", y = expression("SeqDef (auto_max " * lambda * ")"), tag = "A") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
lam_sum <- lam_tab %>% group_by(lambda) %>% summarise(m = median(rho), lo = quantile(rho, .025), hi = quantile(rho, .975))
pB <- ggplot(lam_sum, aes(lambda, m)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1862C9", alpha = .2) +
  geom_line(colour = "#1862C9") + geom_point(colour = "#1862C9") + scale_x_log10() +
  geom_vline(xintercept = median(by_tree$lambda), linetype = "dashed", colour = "#D73027") +
  labs(x = expression(lambda * " (log scale)"), y = expression("Spearman " * rho * "(SeqDef, nearest-sequenced distance)"), tag = "B") + theme_minimal(base_size = 10)
ggsave("figures/figS4_simple_measures.pdf", patchwork::wrap_plots(pA, pB, nrow = 1), width = 10, height = 4.3)
cat("\n[simple measures] DONE.\n")
