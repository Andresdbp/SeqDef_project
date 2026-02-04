# ============================================
# ANALYSIS 1: GLOBAL CHONDRICHTHYES PRIORITIZATION
# ============================================
# Workflow:
# 1. Load IUCN & NCBI Data (or download via API)
# 2. Load Posterior Tree Distribution (100 trees)
# 3. Run SeqDef (Auto-Lambda) on ALL trees
# 4. Aggregate Scores (Median & CI)
# 5. Visualize on MCC Tree
# ============================================

# Libraries
library(iucnredlist)
library(ape)
library(phangorn)
library(taxize)
library(rentrez)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(viridis)
library(scales)
library(stringr)

# ============================================
# STEP 1: IUCN RED LIST DATA
# ============================================
# To run from scratch, uncomment the API calls below:

# Sys.setenv(
#   ENTREZ_KEY = "9ddac7f8b6b68674a120d97191ad836b3908",
#   NCBI_API_KEY = "9ddac7f8b6b68674a120d97191ad836b3908"
# )
# rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))
# api <- init_api("JHxnLuf9N7182D578jLL2QPWYjrfz9ZjMF8Y")
# CHANGE BEFORE PUBLISHING *********************************************
# Sys.setenv(IUCN_REDLIST_KEY = "YOUR_KEY_HERE")
# api <- init_api(Sys.getenv("IUCN_REDLIST_KEY"))
# 
# ids <- assessments_by_taxonomy(
#   api, level = "class", name = "chondrichthyes",
#   latest = TRUE, scope_code = 1, wait_time = 2
# )
# a_data <- assessment_data_many(api, ids$assessment_id, wait_time = 1)
# saveRDS(a_data, "data/raw/iucn_assessment_data.rds")

# --- LOAD DATA ---
a_data <- readRDS("data/raw/iucn_assessment_data.rds")

iucn_clean <- map_dfr(a_data, function(el) {
  tibble(
    scientific_name = tryCatch(el$taxon$scientific_name[1], error = function(e) NA_character_),
    iucn_category   = tryCatch(el$red_list_category$code[1],  error = function(e) NA_character_),
    order           = tryCatch(el$taxon$order_name[1],        error = function(e) NA_character_)
  )
}) %>%
  filter(!is.na(scientific_name), !is.na(iucn_category), !iucn_category %in% c("EX", "DD", "NE")) %>%
  mutate(scientific_name = str_squish(scientific_name)) %>%
  distinct(scientific_name, .keep_all = TRUE)

# ============================================
# STEP 2: LOAD TREES (Posterior & MCC)
# ============================================
chond <- read.nexus("data/chondrichthyes.nex") # 100 Posterior Trees
tree_mcc <- phangorn::maxCladeCred(chond) # For Plotting Topology

# Normalize Names (Underscores for consistency)
tree_mcc$tip.label <- gsub(" ", "_", str_squish(tree_mcc$tip.label))
iucn_clean$scientific_name <- gsub(" ", "_", iucn_clean$scientific_name)

# ============================================
# STEP 3: NCBI ASSEMBLY DATA
# ============================================
# To run from scratch, uncomment the API calls below:

# Sys.setenv(ENTREZ_KEY = "YOUR_KEY_HERE")
# rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))
# 
# species_to_check <- unique(tree_mcc$tip.label)
# ncbi_results <- tibble(scientific_name = species_to_check, assembly_availability = 0L)
# 
# chunk_size <- 50
# n_chunks <- ceiling(length(species_to_check) / chunk_size)
# 
# for (i in seq_len(n_chunks)) {
#   idx <- ((i - 1) * chunk_size + 1):min(i * chunk_size, length(species_to_check))
#   chunk_species <- species_to_check[idx]
#   # ... (Insert original NCBI loop logic here) ...
#   # (See previous code for full loop details if re-implementing)
# }
# saveRDS(ncbi_results, "data/raw/ncbi_assembly_data.rds")

# --- LOAD DATA ---
ncbi_results <- readRDS("data/raw/ncbi_assembly_data.rds")
ncbi_results$scientific_name <- gsub(" ", "_", ncbi_results$scientific_name)

# ============================================
# STEP 4: RUN SEQDEF ON POSTERIOR DISTRIBUTION
# ============================================
cat("Calculating SeqDef across", length(chond), "posterior trees...\n")

posterior_scores <- purrr::map_dfr(1:length(chond), function(i) {
  phy <- chond[[i]]
  phy$tip.label <- gsub(" ", "_", str_squish(phy$tip.label))
  
  # Intersect Data
  common_taxa <- intersect(phy$tip.label, iucn_clean$scientific_name)
  if(length(common_taxa) < 10) return(tibble())
  
  phy <- ape::keep.tip(phy, common_taxa)
  
  # Build Input (Availability 0/1)
  input_df <- tibble(
    taxa = phy$tip.label,
    score = ncbi_results$assembly_availability[match(phy$tip.label, ncbi_results$scientific_name)]
  ) %>% mutate(score = replace_na(score, 0))
  
  # Run SeqDef (Auto Mode)
  res <- suppressMessages(SeqDef(phy, input_df, lambda = "auto_max"))
  
  tibble(tree_id = i, taxa = names(res$seqdef), seqdef_score = as.numeric(res$seqdef))
}, .progress = TRUE)

# Aggregate Scores (Median)
agg_scores <- posterior_scores %>%
  group_by(taxa) %>%
  summarise(median_seqdef = median(seqdef_score), .groups="drop")

# ============================================
# STEP 5: PREPARE PLOTTING DATA & CALCULATE PRIORITY
# ============================================
plot_data <- tibble(tip_label = tree_mcc$tip.label) %>%
  filter(tip_label %in% agg_scores$taxa) %>%
  left_join(agg_scores, by = c("tip_label" = "taxa")) %>%
  left_join(iucn_clean, by = c("tip_label" = "scientific_name")) %>%
  left_join(ncbi_results, by = c("tip_label" = "scientific_name")) %>%
  mutate(
    iucn_category = factor(iucn_category, levels = c("CR", "EN", "VU", "NT", "LC")),
    assembly_availability = replace_na(assembly_availability, 0L)
  )

# Calculate Priority using the package function 'calc_priority' inline
risk_index <- c(LC = 0, NT = 1, VU = 2, EN = 3, CR = 4)

plot_data <- plot_data %>%
  mutate(
    GE_score = replace_na(risk_index[as.character(iucn_category)], 0),
    
    # --- NEW ONE-LINER IMPLEMENTATION ---
    priority = calc_priority(
      seqdef_res   = list(seqdef = setNames(median_seqdef, tip_label)), 
      trait_values = setNames(GE_score, tip_label), 
      model        = "exponential"
    )[tip_label]
  )

# Save Final Results
# write_csv(plot_data, "output/chondrichthyes_final_priority.csv")

# ============================================
# STEP 6: VISUALIZATION
# ============================================

# Prepare tree
tree_plot <- ape::keep.tip(tree_mcc, plot_data$tip_label)
tree_plot <- ladderize(tree_plot, right = TRUE)

# Colors & Subsets
iucn_cols <- c(CR="#C9184A", EN="#FF4D6D", VU="#FF758F", NT="#FFB3C1", LC="#FFCCD5")
tips_with_asm <- plot_data$tip_label[plot_data$assembly_availability == 1]
top_priority <- plot_data %>% arrange(desc(priority)) %>% slice(1) %>% pull(tip_label)

# 1. Base Tree
p <- ggtree(tree_plot, layout="circular", size=0.15, color="#48bf91", open.angle=0) +
  theme_tree() + 
  theme(legend.position="right", legend.box="vertical") +
  
  # RING 1: Median SeqDef
  geom_fruit(data=plot_data, geom=geom_col, 
             mapping=aes(y=tip_label, x=median_seqdef, fill=median_seqdef),
             pwidth=0.3, offset=0.05, 
             axis.params=list(axis="x", text=NA, line=NA, ticks=NA)) +
  scale_fill_gradient2(name="SeqDef (Median)", low="#CCDCFF", mid="#75A1FF", high="#1862C9", midpoint=0.5,
                       guide=guide_colorbar(order=1)) +
  
  # RING 2: IUCN
  new_scale_fill() +
  geom_fruit(data=plot_data, geom=geom_col, 
             mapping=aes(y=tip_label, x=1, fill=iucn_category),
             pwidth=0.15, offset=0.02) +
  scale_fill_manual(name="IUCN", values=iucn_cols, na.value="grey90", guide=guide_legend(order=2)) +
  
  # RING 3: Priority Score
  new_scale_fill() +
  geom_fruit(data=plot_data, geom=geom_col, 
             mapping=aes(y=tip_label, x=priority, fill=priority),
             pwidth=0.3, offset=0.02) +
  scale_fill_gradient(name="Priority Score", low="#FFEACC", high="#C96E18", guide=guide_colorbar(order=3)) +
  
  # --- MARKERS (Legend Setup) ---
  new_scale_color() +
  
  # A. REAL GREY DOTS (Visible on tips)
  geom_tippoint(aes(subset=(label %in% tips_with_asm), color="Has assembly"), size=1.5) +
  
  # B. DUMMY RED DOTS (Invisible on tips, used ONLY to generate the legend)
  geom_tippoint(aes(subset=(label %in% top_priority), color="Top Priority"), 
                alpha=0, size=0) + 
  
  # C. Legend Definitions
  scale_color_manual(name="Status", values=c("Has assembly"="grey40", "Top Priority"="#D73027")) +
  
  # D. Force Legend Visibility (Make the invisible dot visible in key)
  guides(color = guide_legend(order=4, override.aes=list(
    alpha=1, size=3, shape=19
  )))

# --- OUTER RED MARKER (The Real One) ---
# We add this LAST so it appears outside the Priority ring
p <- p +
  new_scale_color() + 
  geom_fruit(
    data = plot_data %>% filter(tip_label %in% top_priority),
    geom = geom_point,
    mapping = aes(y = tip_label, x = 1), # x=1 centers it in this track
    orientation = "y",
    pwidth = 0,      # Width of the track (0 forces it to hug the previous ring)
    offset = 0.05,   # Spacing from the Priority ring
    size = 3,        # Size of the red dot
    color = "#D73027",
    shape = 19
  )

print(p)


# ## SeqDef Example
# # 1. Tree
# tree_text <- "(taxa1:100,(taxa2:90,(taxa3:80,(taxa4:70,(taxa5:60,(taxa6:50,(taxa7:40,(taxa8:30,(taxa9:10,taxa10:10):20):10):10):10):10):10):10):10);"
# tr <- read.tree(text = tree_text)
# 
# # 2. Data
# data <- data.frame(
#   label = paste0("taxa", 1:10),
#   Availability = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
#   SeqDef = c(1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0)
# )
# data_ordered <- data[match(tr$tip.label, data$label), ]
# 
# # 3. Plot
# plot(tr, 
#      x.lim = 170, 
#      y.lim = c(1, 12), 
#      align.tip.label = TRUE, 
#      font = 3, 
#      label.offset = 2, 
#      edge.width = 2)
# 
# # 4. Add Data Columns
# y_coords <- 1:length(tr$tip.label)
# 
# # Column 1
# text(x = 135, y = y_coords, labels = data_ordered$Availability)
# text(x = 135, y = 10.9, labels = "Sequencing\nAvailability", font = 2, cex = 0.8, adj = c(0.5, 0))
# 
# # Column 2
# text(x = 160, y = y_coords, labels = sprintf("%.1f", data_ordered$SeqDef))
# text(x = 160, y = 11, labels = "SeqDef. Statistic", font = 2, cex = 0.8, adj = c(0.5, 0))