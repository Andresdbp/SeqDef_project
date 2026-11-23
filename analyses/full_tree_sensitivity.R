# =============================================================================
# full_tree_sensitivity.R
# R2 Reviewer 2, Major 4 (+ "run all sensitivity tests"): in the main analysis
# SeqDef is computed on the tree pruned to the 877 species with an exact IUCN
# name match in categories LC-CR. Pruning also removed any genome carried by a
# Data Deficient or name-mismatched tip. Here SeqDef is computed on ALL 1,192
# tips, with S taken from data/ncbi_assembly_all_tips.csv (same NCBI snapshot,
# 2025-12-30, mapped through NCBI and IUCN synonyms), and Priority is formed for
# every assessed (LC-CR) tip, matched exactly or through an IUCN synonym.
# All 100 posterior trees, auto_max lambda.
# Output: results/full_tree_sensitivity.csv, results/full_tree_top1.csv
# =============================================================================
source("analyses/00_setup.R")
allS <- read.csv("data/ncbi_assembly_all_tips.csv", stringsAsFactors = FALSE)
unm  <- read.csv("results/provenance_unmatched_tips.csv", stringsAsFactors = FALSE)
posterior <- readRDS("results/posterior_priority.rds")
TOL <- 1e-6; top_set <- function(x) names(x)[x >= max(x) * (1 - TOL)]
tips <- clean_tree(chond[[1]])$tip.label

# IUCN category per tip: exact accepted-name match, else an IUCN synonym whose
# accepted name is not itself a tip (avoids double-counting lumped taxa)
cat_exact <- iucn_clean$iucn_category[match(tips, iucn_clean$scientific_name)]
acc <- unm$iucn_accepted_name[match(tips, unm$tip)]
cat_syn <- ifelse(!is.na(acc) & !acc %in% tips, iucn_clean$iucn_category[match(acc, iucn_clean$scientific_name)], NA)
cat_all <- setNames(ifelse(is.na(cat_exact), cat_syn, cat_exact), tips)
assessed <- names(cat_all)[!is.na(cat_all)]
cat(sprintf("tips: %d | with assembly: %d | assessed LC-CR: %d (exact %d + via IUCN synonym %d)\n",
            length(tips), sum(allS$assembly_availability), length(assessed), sum(!is.na(cat_exact)), sum(is.na(cat_exact) & !is.na(cat_syn))))

res <- map_dfr(seq_along(chond), function(i) {
  tr <- clean_tree(chond[[i]])
  df <- tibble(taxa = tr$tip.label, S = allS$assembly_availability[match(tr$tip.label, allS$tip)])
  r  <- suppressMessages(SeqDef(tr, df, lambda = "auto_max"))
  ge <- risk_index[cat_all[assessed]]
  pri <- r$seqdef[assessed] * 2^ge; names(pri) <- assessed
  old <- posterior[[which(map_int(posterior, "tree_id") == i)]]
  common <- intersect(names(pri), names(old$pri_auto))
  t1 <- sort(top_set(pri))
  tibble(tree_id = i, lambda = r$lambda, lambda_R1 = old$lambda_auto,
         rho_priority_vs_R1 = cor(pri[common], old$pri_auto[common], method = "spearman"),
         rho_seqdef_vs_R1 = cor(r$seqdef[common], old$seqdef_auto[common], method = "spearman"),
         top1 = paste(t1, collapse = ";"), n_top1 = length(t1), target_top1 = identical(t1, TARGET),
         target_seqdef = r$seqdef[TARGET], target_rank = rank(-pri, ties.method = "min")[TARGET])
}, .progress = TRUE)
write.csv(res, "results/full_tree_sensitivity.csv", row.names = FALSE)
q <- function(v, d = 3) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), median(v), quantile(v, .025), quantile(v, .975))
cat(sprintf("lambda (full tree) %s vs R1 %s\n", q(res$lambda, 2), q(res$lambda_R1, 2)))
cat(sprintf("rho Priority vs R1 (common species): %s\n", q(res$rho_priority_vs_R1)))
cat(sprintf("rho SeqDef   vs R1 (common species): %s\n", q(res$rho_seqdef_vs_R1)))
cat(sprintf("C. atromarginatus unique top-1: %d / 100 trees (R1: 85); median rank %s\n", sum(res$target_top1), median(res$target_rank)))
lab <- function(w) { w <- strsplit(w, ";")[[1]]; if (length(w) == 1) w else sprintf("%d-way tie: %s", length(w), paste(unique(sub("_.*", "", w)), collapse = "+")) }
top <- tibble(label = map_chr(res$top1, lab)) %>% count(label, sort = TRUE)
print(top); write.csv(top, "results/full_tree_top1.csv", row.names = FALSE)
