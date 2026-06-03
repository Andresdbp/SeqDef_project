# =============================================================================
# 00_setup.R  -- Shared setup for the JEB revision analyses
# Run everything from the project root: SeqDef_project/
#
# IMPORTANT: This never sources analysis.R (which contains API keys). Data are
# loaded directly from the cached .rds / .nex files. No network calls.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

set.seed(42)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- Kernel-enabled SeqDef() + calc_priority() (mirror of the package source) -
source("function.R")

# --- Data --------------------------------------------------------------------
chond  <- read.nexus("data/chondrichthyes.nex")        # 100 posterior trees
ncbi   <- readRDS("data/ncbi_assembly_data.rds")        # NCBI assembly table
a_data <- readRDS("data/iucn_assessment_data.rds")      # IUCN assessments (list)

# IUCN clean table (same extraction logic as analysis.R)
iucn_clean <- map_dfr(a_data, function(el) {
  tibble(
    scientific_name = tryCatch(el$taxon$scientific_name[1], error = function(e) NA_character_),
    iucn_category   = tryCatch(el$red_list_category$code[1], error = function(e) NA_character_),
    order           = tryCatch(el$taxon$order_name[1],       error = function(e) NA_character_)
  )
}) %>%
  filter(!is.na(scientific_name), !is.na(iucn_category),
         !iucn_category %in% c("EX", "DD", "NE")) %>%
  mutate(scientific_name = str_squish(scientific_name)) %>%
  distinct(scientific_name, .keep_all = TRUE)

# Normalize names to underscores everywhere
iucn_clean$scientific_name <- gsub(" ", "_", iucn_clean$scientific_name)
ncbi$scientific_name       <- gsub(" ", "_", str_squish(ncbi$scientific_name))

# IUCN -> Global Endangerment weight
risk_index <- c(LC = 0, NT = 1, VU = 2, EN = 3, CR = 4)

TARGET <- "Centrophorus_atromarginatus"   # the case-study top target

# --- Helpers -----------------------------------------------------------------
clean_tree <- function(phy) {
  phy$tip.label <- gsub(" ", "_", str_squish(phy$tip.label))
  phy
}

# Prune a tree to the IUCN-intersection and attach binary assembly availability
build_input_binary <- function(phy) {
  phy <- clean_tree(phy)
  common <- intersect(phy$tip.label, iucn_clean$scientific_name)
  if (length(common) < 10) return(NULL)
  phy <- keep.tip(phy, common)
  s <- ncbi$assembly_availability[match(phy$tip.label, ncbi$scientific_name)]
  s[is.na(s)] <- 0
  list(tree = phy, df = tibble(taxa = phy$tip.label, score = s))
}

# GE vector (named) aligned to a tree's tips
ge_for_tree <- function(phy) {
  cats <- iucn_clean$iucn_category[match(phy$tip.label, iucn_clean$scientific_name)]
  ge   <- risk_index[as.character(cats)]
  ge[is.na(ge)] <- 0
  setNames(ge, phy$tip.label)
}

# --- MCC tree (cached: maxCladeCred is the slow step) ------------------------
mcc_path <- "results/mcc_tree.rds"
if (file.exists(mcc_path)) {
  tree_mcc <- readRDS(mcc_path)
} else {
  message("Computing MCC tree (cached after first run)...")
  tree_mcc <- clean_tree(phangorn::maxCladeCred(chond))
  saveRDS(tree_mcc, mcc_path)
}

# --- Session / hardware record ----------------------------------------------
write_session <- function() {
  ram <- tryCatch(round(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1024^3, 1),
                  error = function(e) NA)
  con <- file("results/SESSION.txt", "w")
  on.exit(close(con))
  writeLines(c(
    paste("Generated:", format(Sys.time())),
    paste("Machine  :", paste(Sys.info()[c("sysname", "release", "machine", "nodename")], collapse = " ")),
    paste("Cores    :", parallel::detectCores()),
    paste("RAM (GB) :", ram),
    "",
    "--- sessionInfo() ---",
    capture.output(sessionInfo())
  ), con)
}

message(sprintf("Setup loaded: %d posterior trees, %d IUCN species, %d NCBI rows, MCC tips = %d",
                length(chond), nrow(iucn_clean), nrow(ncbi), length(tree_mcc$tip.label)))
