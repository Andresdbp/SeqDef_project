# ============================================
# DEFINITION: SeqDef Function
# ============================================
SeqDef <- function(tree, df, data.col = 2, invert = TRUE, scale = TRUE, lambda = "auto_max"){
  
  # 1. Convert & Align Data
  df <- as.data.frame(df)
  
  # Prune tree if it has tips not in data
  common_taxa <- intersect(tree$tip.label, df[,1])
  if(length(common_taxa) < length(tree$tip.label)){
    # message("Pruning tree to match dataframe...")
    tree <- ape::keep.tip(tree, common_taxa)
  }
  
  # Align dataframe to match tree order exactly
  df <- df[match(tree$tip.label, df[,1]), ]
  
  # 2. Calculate Tree Depth & Matrix
  # We need tree depth (td) to normalize distances to 0-1
  td <- max(ape::branching.times(tree))
  dist.matrix <- ape::cophenetic.phylo(tree)
  dist.matrix <- dist.matrix[tree$tip.label, tree$tip.label]
  
  # 3. Handle Lambda Selection (Manual vs Auto)
  final_lambda <- 0
  
  if (length(lambda) == 1 && lambda == "auto_max") {
    # --- MAXIMIZE VARIANCE LOGIC (With 10% Drop Tolerance) ---
    message("Optimizing Lambda: Searching for peak variance (Stopping if variance drops >10% below peak)...")
    
    # 1. Prepare Data
    s_vec <- df[[data.col]]
    norm_dists <- dist.matrix / td
    
    # 2. Iterative Scan
    lambda_seq <- seq(1, 50, 0.1)
    
    # Initialize with the first lambda value
    best_lambda <- lambda_seq[1]
    max_var <- -1
    
    # Helper to calc variance for a single lambda
    calc_var <- function(x) {
      w_mat <- exp(-x * norm_dists)
      raw_scores <- as.numeric(w_mat %*% s_vec)
      rng <- range(raw_scores)
      if(rng[2] - rng[1] < 1e-9) return(0)
      norm_scores <- (raw_scores - rng[1]) / (rng[2] - rng[1])
      return(var(1 - norm_scores))
    }
    
    # Calculate baseline for the first step
    max_var <- calc_var(best_lambda)
    
    # Start loop from the second value
    for(lam in lambda_seq[-1]) {
      curr_var <- calc_var(lam)
      
      if(curr_var > max_var) {
        # NEW PEAK FOUND: Update record and continue
        max_var <- curr_var
        best_lambda <- lam
        
      } else {
        # DROP DETECTED: Check if it is significant (>10% drop from the PEAK)
        drop_ratio <- (max_var - curr_var) / max_var
        
        if(drop_ratio > 0.10) {
          # STOP: We have fallen more than 10% from the peak. 
          # We assume the peak is real and we are now on the downward slope.
          # message(sprintf("  Stopping at Lambda=%.1f (Variance dropped %.1f%% below peak)", lam, drop_ratio*100))
          break
        }
        # If drop is < 10%, we treat it as noise/fluctuation and keep searching.
      }
    }
    
    final_lambda <- best_lambda
    message(sprintf("Selected Lambda: %.1f (Peak Variance: %.5f)", final_lambda, max_var))
    
  } else if (length(lambda) == 1 && lambda == "by_genus") {
    # --- GENUS-SCALE LOGIC ---
    
    # A. Extract Genus names
    genera <- sapply(strsplit(tree$tip.label, "[_ ]"), `[`, 1)
    unique_genera <- unique(genera)
    
    # B. Collect pairwise distances WITHIN Genera
    intra_genus_dists <- c()
    
    for(g in unique_genera) {
      tips_in_genus <- tree$tip.label[genera == g]
      
      # We need at least 2 species to measure a distance
      if(length(tips_in_genus) > 1) {
        # Subset distance matrix directly
        sub_mat <- dist.matrix[tips_in_genus, tips_in_genus]
        
        # Get unique distances > 0
        d_vals <- sub_mat[lower.tri(sub_mat)]
        d_vals <- d_vals[d_vals > 0]
        
        # Normalize by Tree Depth
        intra_genus_dists <- c(intra_genus_dists, d_vals / td)
      }
    }
    
    # C. Calculate Half-Life Lambda
    if(length(intra_genus_dists) > 0){
      median_dist <- median(intra_genus_dists)
      final_lambda <- log(2) / median_dist
      message(sprintf("Genus-Scale Lambda: %.2f (Half-life at median genus distance: %.4f)", 
                      final_lambda, median_dist))
    } else {
      # Edge case: Tree is genus-level (1 species per genus)
      warning("Biological calibration failed: No genera with >1 species found (input tree appears to be genus-level). Defaulting to Lambda=3. Consider using lambda='auto_max' for this dataset.")
      final_lambda <- 3
    }
    
  } else {
    # --- MANUAL LOGIC ---
    if(!is.numeric(lambda)) stop("Lambda must be 'auto_max', 'by_genus', or numeric.")
    final_lambda <- lambda
  }
  
  # 4. Kernel Calculation (Exponential Decay)
  # Note: dist.matrix / td scales input distances to 0-1
  dist.prop <- exp(-final_lambda * dist.matrix / td)
  
  # 5. Vectorized Calculation
  raw_scores <- dist.prop %*% as.numeric(df[, data.col])
  
  # Convert to named vector
  synscores <- as.vector(raw_scores)
  names(synscores) <- tree$tip.label
  
  # 6. Scale (Normalize to 0-1)
  if(scale){
    rng <- max(synscores, na.rm = TRUE) - min(synscores, na.rm = TRUE)
    if(rng == 0) {
      synscores[] <- 0 
    } else {
      synscores <- (synscores - min(synscores, na.rm = TRUE)) / rng
    }
  }
  
  # 7. Invert (Availability -> Deficiency)
  if(invert){
    synscores <- 1 - synscores
  }
  
  # 8. Return Results
  # We include the used 'lambda' in the output for reference
  results <- list(
    tree = tree, 
    seqdef = synscores, 
    empdata = df[, data.col], 
    lambda = final_lambda
  )
  class(results) <- "seqdef"
  return(results)
}

calc_priority <- function(seqdef_res, trait_values, model = c("exponential", "linear", "additive"), 
                          base = 2, na.fill = 0) {
  
  # 1. Input Validation
  model <- match.arg(model)
  
  # Validate that input looks like a SeqDef result
  if (!is.list(seqdef_res) || is.null(seqdef_res$seqdef)) {
    stop("seqdef_res must be a list containing a '$seqdef' element (output from SeqDef).")
  }
  
  # Extract scores and align names
  s_scores <- seqdef_res$seqdef
  target_taxa <- names(s_scores)
  
  if (is.null(target_taxa)) {
    stop("The 'seqdef' vector in seqdef_res must have names matching the taxa.")
  }
  
  # 2. Process Trait Values
  # Convert data frame to named vector if necessary
  if (is.data.frame(trait_values)) {
    # Assume col 1 is names, col 2 is values
    t_vec <- setNames(trait_values[[2]], trait_values[[1]])
  } else {
    t_vec <- trait_values
  }
  
  # Align trait vector to SeqDef vector
  # (Fill missing species with na.fill)
  aligned_traits <- t_vec[target_taxa]
  aligned_traits[is.na(aligned_traits)] <- na.fill
  
  # 3. Calculate Priority
  if (model == "exponential") {
    # Formula: Score * (Base ^ Trait)
    priority <- s_scores * (base ^ aligned_traits)
    
  } else if (model == "linear") {
    # Formula: Score * Trait
    priority <- s_scores * aligned_traits
    
  } else if (model == "additive") {
    # Formula: Score + Trait
    priority <- s_scores + aligned_traits
  }
  
  return(priority)
}