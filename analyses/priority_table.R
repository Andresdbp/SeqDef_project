# =============================================================================
# priority_table.R
# The full ranked table of the case study (cited in the Results): one row per
# species with posterior-median SeqDef, its 2.5-97.5% quantiles across the 100
# posterior trees, Priority = median SeqDef x 2^GE (Eq. in Methods), the
# quantiles of per-tree Priority, the number of trees in which the species is
# (or ties for) the top Priority, family, assembly status and placement type
# (molecular / imputed, Stein et al. 2018).
# Output: results/priority_ranking_posterior.csv
# =============================================================================
source("analyses/00_setup.R")
posterior <- readRDS("results/posterior_priority.rds")
TOL <- 1e-6
sp <- names(posterior[[1]]$seqdef_auto)
SD <- sapply(posterior, function(p) p$seqdef_auto[sp])
PR <- sapply(posterior, function(p) p$pri_auto[sp])
top <- rowSums(sapply(posterior, function(p) { x <- p$pri_auto[sp]; x >= max(x) * (1 - TOL) }))
fam <- map_chr(sp, function(s) { i <- which(map_chr(a_data, ~ gsub(" ", "_", str_squish(.x$taxon$scientific_name[1]))) == s)
  if (length(i)) str_to_title(a_data[[i[1]]]$taxon$family_name[1]) else NA_character_ })
st <- if (file.exists("data/stein2018_data_type.csv")) read.csv("data/stein2018_data_type.csv") else NULL
placement <- if (!is.null(st)) ifelse(st$data_type[match(sp, gsub(" ", "_", st$species))] == "molecular", "molecular", "imputed") else NA
cat_ <- iucn_clean$iucn_category[match(sp, iucn_clean$scientific_name)]
ge <- risk_index[cat_]
tab <- tibble(species = sp, family = fam, iucn = cat_, GE = unname(ge),
              has_assembly = ncbi$assembly_availability[match(sp, ncbi$scientific_name)],
              placement = placement,
              seqdef_median = apply(SD, 1, median), seqdef_q025 = apply(SD, 1, quantile, .025), seqdef_q975 = apply(SD, 1, quantile, .975),
              priority = apply(SD, 1, median) * 2^unname(ge),
              priority_q025 = apply(PR, 1, quantile, .025), priority_q975 = apply(PR, 1, quantile, .975),
              trees_top1 = top) %>%
  arrange(desc(priority)) %>% mutate(rank = rank(-signif(priority, 6), ties.method = "min"), .before = 1)  # medians equal to 6 significant digits (clade members) share a rank
write.csv(tab, "results/priority_ranking_posterior.csv", row.names = FALSE)
print(head(tab, 12), width = 200)
