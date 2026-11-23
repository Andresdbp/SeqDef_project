# =============================================================================
# data_provenance.R
# R2 Reviewer 2, Major 4: name matching and losses at each join, DD/NE handling,
# number of tips coded S = 1, assembly-count distribution, IUCN assessment
# dates, and (if the Stein et al. 2018 genetic-data list is supplied)
# the number of imputed tips.  Offline; uses only cached data.
#
# Optional input: data/stein2018_genetic_species.csv  -- one column 'species'
#   listing the 610 species placed with molecular data in Stein et al. (2018)
#   (from their Supplementary Data / VertLife).  If absent, that part is skipped.
#
# Produces: results/provenance.txt, results/provenance_unmatched_tips.csv
# =============================================================================
source("analyses/00_setup.R")
out <- character(); say <- function(...) { line <- sprintf(...); cat(line, "\n"); out <<- c(out, line) }

tips <- clean_tree(chond[[1]])$tip.label
trees <- lapply(seq_along(chond), function(i) chond[[i]])   # [[ ]] restores tip labels of a compressed multiPhylo
same_tips <- all(vapply(trees, function(p) setequal(clean_tree(p)$tip.label, tips), logical(1)))
spread <- vapply(trees, function(p) { d <- node.depth.edgelength(p)[seq_len(Ntip(p))]; diff(range(d)) / max(d) }, numeric(1))
raw <- map_dfr(a_data, function(el) tibble(
  name = gsub(" ", "_", str_squish(el$taxon$scientific_name[1])),
  cat = el$red_list_category$code[1],
  family = str_to_title(el$taxon$family_name[1]),
  year_published = as.character(el$year_published[1]),
  assessment_date = substr(as.character(el$assessment_date[1]), 1, 10),
  syn = list(tryCatch(gsub(" ", "_", str_squish(paste(el$taxon_synonyms$genus_name, el$taxon_synonyms$species_name))),
                      error = function(e) character()))))

say("== PHYLOGENY ==")
say("VertLife/Stein et al. 2018 posterior: %d trees, %d tips each (identical tip sets across trees: %s)", length(chond), length(tips), same_tips)
say("All trees binary: %s; ultrametric up to rounding (max relative spread of root-to-tip depths: %.1e)", all(vapply(trees, is.binary, logical(1))), max(spread))
say("== IUCN (class Chondrichthyes, latest global assessments) ==")
say("assessment records: %d (unique names: %d)", nrow(raw), n_distinct(raw$name))
tab <- table(factor(raw$cat, levels = c("CR", "EN", "VU", "NT", "LC", "DD", "NE", "EX")))
say("categories: %s", paste(names(tab), tab, sep = "=", collapse = ", "))
acc <- unique(regmatches(raw_cit <- sapply(a_data, function(e) as.character(e$citation[1])), regexpr("Accessed on [0-9]+ [A-Za-z]+ [0-9]{4}", raw_cit)))
say("IUCN access date(s) recorded in the assessment citations: %s (Red List version current on that date: 2025-2, released 10 Oct 2025)", paste(acc, collapse = "; "))
say("assessment years published: %s to %s; most recent assessment date: %s",
    min(raw$year_published, na.rm = TRUE), max(raw$year_published, na.rm = TRUE), max(raw$assessment_date, na.rm = TRUE))
say("== NAME MATCHING (exact string match after whitespace squish, space -> underscore; no synonym resolution) ==")
m_any <- tips %in% raw$name
say("tree tips matching an IUCN accepted name (any category): %d / %d", sum(m_any), length(tips))
say("tree tips with no exact IUCN match (lost at the tree-IUCN join): %d", sum(!m_any))
cat_of <- raw$cat[match(tips, raw$name)]
say("matched tips excluded as DD: %d; NE: %d; EX: %d", sum(cat_of == "DD", na.rm = TRUE), sum(cat_of == "NE", na.rm = TRUE), sum(cat_of == "EX", na.rm = TRUE))
keep <- !is.na(cat_of) & cat_of %in% c("CR", "EN", "VU", "NT", "LC")
say("analysed species (tree n IUCN, categories LC-CR): %d", sum(keep))
iu_ok <- raw$name[raw$cat %in% c("CR", "EN", "VU", "NT", "LC")]
say("IUCN-assessed species (LC-CR) absent from the tree by exact name: %d", sum(!iu_ok %in% tips))
# synonym rescue diagnostic
syn_long <- raw %>% select(name, cat, syn) %>% unnest(syn) %>% filter(syn %in% tips[!m_any])
say("of the %d unmatched tips, %d match an IUCN synonym (DD among them: %d); %d of these synonyms point to an accepted name that is itself a separate tip (taxonomic lumping)",
    sum(!m_any), n_distinct(syn_long$syn), n_distinct(syn_long$syn[syn_long$cat == "DD"]), n_distinct(syn_long$syn[syn_long$name %in% tips]))
write.csv(tibble(tip = tips[!m_any]) %>% left_join(syn_long %>% distinct(syn, .keep_all = TRUE), by = c("tip" = "syn")) %>%
            rename(iucn_accepted_name = name, iucn_category = cat), "results/provenance_unmatched_tips.csv", row.names = FALSE)
say("== DD / NE handling ==")
say("DD and NE species are removed before the intersection (no Priority weight can be assigned); they are therefore also absent from the tree on which SeqDef is computed.")
say("== NCBI ASSEMBLIES ==")
say("NCBI table rows: %d (queried for the analysed species only); tips coded S = 1: %d; S = 0: %d", nrow(ncbi), sum(ncbi$assembly_availability == 1), sum(ncbi$assembly_availability == 0))
ac <- table(ncbi$assembly_count[ncbi$assembly_count > 0]); say("assembly count per sequenced species: %s", paste(names(ac), ac, sep = ":", collapse = ", "))
say("NCBI query: Assembly database, all records for the species (any level, GenBank/RefSeq, any status) released on or before 2025-12-30 (R1 table compiled 30 Dec 2025); reproduced on 2026-11-23 with ncbi_query_all_tips.R (term 'Chondrichthyes[Organism]') -- see ncbi_map_all_tips.R for the reconciliation")
fam <- raw$family[match(tips[keep], raw$name)]; s_k <- ncbi$assembly_availability[match(tips[keep], ncbi$scientific_name)]
fh <- tapply(s_k, fam, max)
say("families in analysis: %d; with >= 1 assembly: %d; without: %d", length(fh), sum(fh == 1), sum(fh == 0))
say("Centrophoridae species in analysis: %d, with assembly: %d", sum(fam == "Centrophoridae"), sum(s_k[fam == "Centrophoridae"]))
say("== IMPUTED TIPS (Stein et al. 2018, Supplementary Table 1, column 'Data type') ==")
if (file.exists("data/stein2018_data_type.csv")) {
  st <- read.csv("data/stein2018_data_type.csv", stringsAsFactors = FALSE)
  st$tip <- gsub(" ", "_", str_squish(st$species))
  dt <- setNames(ifelse(st$data_type == "molecular", "molecular", "imputed"), st$tip)
  say("all 1,192 tips: molecular %d, imputed (taxonomic constraints) %d, not found in table %d",
      sum(dt[tips] == "molecular", na.rm = TRUE), sum(dt[tips] == "imputed", na.rm = TRUE), sum(is.na(dt[tips])))
  an <- tips[keep]
  say("877 analysed tips: molecular %d, imputed %d", sum(dt[an] == "molecular", na.rm = TRUE), sum(dt[an] == "imputed", na.rm = TRUE))
  seqd <- ncbi$scientific_name[ncbi$assembly_availability == 1]
  say("55 sequenced (S = 1) tips: molecular %d, imputed %d", sum(dt[seqd] == "molecular", na.rm = TRUE), sum(dt[seqd] == "imputed", na.rm = TRUE))
  fig3 <- c(TARGET, "Rhinobatos_albomaculatus", "Rhinobatos_jimbaranensis", "Rhinobatos_annandalei", "Acroteriobatus_variegatus",
            "Rhinobatos_rhinobatos", "Rhinobatos_irvinei", "Rhinobatos_lionotus", "Rhinobatos_schlegelii",
            "Squatina_argentina", "Squatina_armata", "Squatina_japonica", "Squatina_legnota", "Squatina_occulta", "Squatina_oculata")
  for (sp in fig3) say("  %-30s %s", sp, dt[sp])
  cen <- an[fam_tbl_family <- raw$family[match(an, raw$name)] == "Centrophoridae"]
  say("Centrophoridae in analysis: %s", paste(sprintf("%s (%s)", cen, dt[cen]), collapse = ", "))
} else say("data/stein2018_data_type.csv not found -- extract Supplementary Table 1 of Stein et al. (2018).")
writeLines(out, "results/provenance.txt")
