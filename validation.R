# ============================================
# ANALYSIS 2: SENSITIVITY & VALIDATION
# ============================================

# ============================================
# TEST 1: DUAL SENSITIVITY ANALYSIS (Figure 2A & 2B)
# ============================================
library(patchwork) # For combining plots

# -------------------------------------------------------
# A. Generate Distributions for BOTH methods
# -------------------------------------------------------
cat("Generating Lambda distributions for Auto-Max and Genus methods...\n")

lambda_data <- map_dfr(1:length(chond), function(i) {
  phy <- chond[[i]]
  phy$tip.label <- gsub(" ", "_", stringr::str_squish(phy$tip.label))
  common <- intersect(phy$tip.label, ncbi_results$scientific_name)
  phy <- keep.tip(phy, common)
  
  input_df <- tibble(
    t = phy$tip.label, 
    s = ncbi_results$assembly_availability[match(phy$tip.label, ncbi_results$scientific_name)]
  )
  
  # Run BOTH methods on this tree
  res_max   <- suppressMessages(SeqDef(phy, input_df, lambda = "auto_max"))
  res_genus <- suppressMessages(SeqDef(phy, input_df, lambda = "by_genus"))
  
  tibble(tree_id = i, 
         val_max = res_max$lambda, 
         val_genus = res_genus$lambda)
}, .progress = TRUE)

# -------------------------------------------------------
# B. Calculate 95% Intervals
# -------------------------------------------------------
get_stats <- function(vec) {
  qs <- quantile(vec, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  list(min = qs["2.5%"], med = qs["50%"], max = qs["97.5%"])
}

stats_max   <- get_stats(lambda_data$val_max)
stats_genus <- get_stats(lambda_data$val_genus)

# -------------------------------------------------------
# C. Calculate Variance Curve (Blue Line)
# Note: We compute this ONCE as the background curve is identical
# -------------------------------------------------------
cat("Calculating variance curve...\n")
sample_trees <- sample(1:length(chond), 20) # 20 trees for speed (increase for final)
lambda_range <- seq(1, 25, 0.2)             # 0.2 step for smoother curve

variance_res <- map_dfr(sample_trees, function(i) {
  phy <- chond[[i]]
  phy$tip.label <- gsub(" ", "_", stringr::str_squish(phy$tip.label))
  common <- intersect(phy$tip.label, ncbi_results$scientific_name)
  phy <- keep.tip(phy, common)
  
  df_sub <- tibble(t=phy$tip.label, 
                   s=ncbi_results$assembly_availability[match(phy$tip.label, ncbi_results$scientific_name)])
  
  map_dfr(lambda_range, function(lam) {
    res <- suppressMessages(SeqDef(phy, df_sub, lambda=lam))
    tibble(tree_id=i, lambda=lam, var_seqdef=var(res$seqdef, na.rm=TRUE))
  })
}, .progress = TRUE)

# ============================================
# D. Create Plots
# ============================================

# Base Plot Function
create_panel <- function(stats, label_text, shade_color) {
  ggplot(variance_res, aes(x=lambda, y=var_seqdef)) +
    # Shaded Region (95% Interval)
    annotate("rect", xmin = stats$min, xmax = stats$max, ymin = -Inf, ymax = Inf,
             fill = shade_color, alpha = 0.2) +
    
    # Blue Curve
    stat_summary(fun.data = "mean_cl_boot", geom = "ribbon", alpha = 0.2, fill = "#1862C9") +
    stat_summary(fun = "mean", geom = "line", color = "#1862C9", linewidth = 1) +
    
    # Median Line
    geom_vline(xintercept = stats$med, color = shade_color, linetype = "longdash", linewidth = 1) +
    
    # A/B Label (Top-Right, Bold)
    # x=Inf, y=Inf places it in the far corner; hjust/vjust pulls it back slightly
    annotate("text", x = Inf, y = Inf, label = label_text, 
             hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold") +
    
    labs(x = "Lambda", y = "Variance in SeqDef") +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.line = element_line(color="black"))
}

# Plot A: Auto-Max (Red)
p1 <- create_panel(stats_max, "A", "#D73027")

# Plot B: Genus-Scale (Purple)
p2 <- create_panel(stats_genus, "B", "#762a83")

# Combine Side-by-Side
final_figure <- p1 + p2
print(final_figure)


# ============================================
# TEST 2: TOPOLOGICAL ROBUSTNESS TEST
# ============================================

get_top_priority_taxa <- function(tree, iucn_df, ncbi_df) {
  
  # 1. Clean Names
  tree$tip.label <- gsub(" ", "_", stringr::str_squish(tree$tip.label))
  iucn_names <- gsub(" ", "_", stringr::str_squish(as.character(iucn_df$scientific_name)))
  ncbi_names <- gsub(" ", "_", stringr::str_squish(as.character(ncbi_df$scientific_name)))
  
  # 2. Intersect Data
  common_taxa <- intersect(tree$tip.label, iucn_names)
  if(length(common_taxa) < 10) return(NA_character_)
  tree <- ape::keep.tip(tree, common_taxa)
  
  # 3. Build Input Dataframe
  s_scores <- ncbi_df$assembly_availability[match(tree$tip.label, ncbi_names)]
  s_scores[is.na(s_scores)] <- 0
  input_df <- tibble(taxa = tree$tip.label, score = s_scores)
  
  # 4. Run SeqDef (UPDATE: Use auto_max)
  res <- suppressMessages(SeqDef(tree, input_df, lambda = "auto_max"))
  
  # 5. Calculate Priority (CHANGE: Use package function)
  risk_index <- c(LC = 0, NT = 1, VU = 2, EN = 3, CR = 4)
  
  # Prepare GE vector
  cats <- iucn_df$iucn_category[match(tree$tip.label, iucn_names)]
  ge_scores <- risk_index[as.character(cats)]
  ge_scores[is.na(ge_scores)] <- 0
  names(ge_scores) <- tree$tip.label
  
  # Use function
  priority_scores <- calc_priority(res, ge_scores, model = "exponential", base = 2)
  
  # 6. Identify Winner
  max_val <- max(priority_scores, na.rm = TRUE)
  return(names(priority_scores)[priority_scores == max_val])
}

# Run the loop
top_taxa_list <- purrr::map(1:length(chond), function(i) {
  get_top_priority_taxa(chond[[i]], iucn_clean, ncbi_results)
}, .progress = TRUE)

# Flatten results
top_taxa_vec <- unlist(top_taxa_list)
top_taxa_vec <- top_taxa_vec[!is.na(top_taxa_vec)]

# Aggregate counts
robustness_counts <- tibble(Species = top_taxa_vec) %>%
  count(Species, name = "Count") %>%
  arrange(desc(Count)) %>%
  mutate(Species = gsub("_", " ", Species))

# Plot
p_rob <- ggplot(robustness_counts, aes(x = reorder(Species, Count), y = Count)) +
  geom_col(fill = "#C96E18", width = 0.7) + 
  coord_flip() +
  labs(x = "", y = "Number of Trees") +
  geom_text(aes(label = Count), hjust = -0.2, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.line.y = element_line(color = "black")
  )

print(p_rob)

