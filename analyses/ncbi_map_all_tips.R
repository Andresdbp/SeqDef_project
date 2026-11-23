# =============================================================================
# ncbi_map_all_tips.R
# Maps the NCBI Assembly snapshot (data/ncbi_assembly_species_summary.tsv,
# produced by ncbi_query_all_tips.R; one row per NCBI species) onto ALL 1,192
# tips of the Stein et al. (2018) trees.
# Matching order for each NCBI species with >= 1 assembly released on or
# before the snapshot date (2025-12-30):
#   1. tip name == NCBI scientific name
#   2. tip name == an NCBI taxonomy synonym of that species (if several tips
#      are NCBI synonyms, the one whose IUCN accepted name is the NCBI name)
#   3. tip whose IUCN accepted name (via IUCN synonym list) == NCBI name
# Also checks that the original 877-species table (data/ncbi_assembly_data.rds)
# is reproduced.
# Output: data/ncbi_assembly_all_tips.csv
# =============================================================================
source("analyses/00_setup.R")
s <- read.delim("data/ncbi_assembly_species_summary.tsv", stringsAsFactors = FALSE)
snap <- s[s$n_assemblies_le_snapshot > 0, ]
tips <- clean_tree(chond[[1]])$tip.label
unm <- read.csv("results/provenance_unmatched_tips.csv", stringsAsFactors = FALSE)   # from data_provenance.R

hits <- map_dfr(seq_len(nrow(snap)), function(i) {
  nm  <- gsub(" ", "_", snap$ncbi_species[i])
  syn <- gsub(" ", "_", strsplit(snap$ncbi_synonyms[i], ";")[[1]])
  t1 <- intersect(nm, tips); t2 <- intersect(syn, tips)
  t3 <- unm$tip[!is.na(unm$iucn_accepted_name) & unm$iucn_accepted_name == nm]
  t23 <- intersect(t2, t3)
  tip <- if (length(t1)) t1 else if (length(t23)) t23[1] else if (length(t2)) t2[1] else if (length(t3)) t3[1] else NA
  how <- if (length(t1)) "exact" else if (length(t2)) "ncbi_synonym" else if (length(t3)) "iucn_synonym" else "not_on_tree"
  tibble(ncbi_species = snap$ncbi_species[i], tip = tip, match = how,
         assembly_count = snap$n_assemblies_le_snapshot[i], best_level = snap$best_level_le_snapshot[i])
})
cat("NCBI species with an assembly by 2025-12-30:", nrow(snap), "\n"); print(table(hits$match))
print(hits %>% filter(match != "exact"))

all_tips <- tibble(tip = tips) %>%
  left_join(hits %>% filter(!is.na(tip)) %>% select(tip, ncbi_species, match, assembly_count, best_level), by = "tip") %>%
  mutate(assembly_count = coalesce(assembly_count, 0L), assembly_availability = as.integer(assembly_count > 0))
write.csv(all_tips, "data/ncbi_assembly_all_tips.csv", row.names = FALSE)

# reproduce the original 877-species table
chk <- inner_join(all_tips, ncbi, by = c("tip" = "scientific_name"))
cat(sprintf("\nOriginal table reproduced (presence/absence) for %d / %d species; assembly counts identical for %d\n",
            sum(chk$assembly_availability.x == chk$assembly_availability.y), nrow(chk),
            sum(chk$assembly_count.x == chk$assembly_count.y)))
print(chk %>% filter(assembly_availability.x != assembly_availability.y | assembly_count.x != assembly_count.y) %>%
        select(tip, ncbi_species, count_now = assembly_count.x, count_R1 = assembly_count.y))
cat("\nTips with an assembly that were OUTSIDE the 877 analysed species in R1:\n")
print(all_tips %>% filter(assembly_availability == 1, !tip %in% ncbi$scientific_name))
