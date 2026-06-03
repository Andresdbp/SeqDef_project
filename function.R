# ============================================
# DEFINITION: SeqDef Function
# ============================================

# ----------------------------------------------------------------------------
# Linear-time, linear-memory availability for the EXPONENTIAL kernel.
#
# A_f = sum_i exp(-lambda * d_fi / td) * s_i, where d_fi is the cophenetic
# (path) distance. Because d_fi is additive along the tree path and the kernel
# is exponential, exp(-lambda d_fi/td) = product of per-edge factors
# phi_e = exp(-lambda L_e / td). A_f is then computed by a two-pass sum-product
# traversal (Felsenstein-style message passing) in O(n) time and O(n) memory --
# no n x n distance matrix is ever formed. This is EXACT (identical to the dense
# result), and only works for the exponential kernel (the only one that
# factorizes into independent per-branch factors).
# ----------------------------------------------------------------------------
.seqdef_exp_avail <- function(tree, s, lambda, td) {
  n   <- length(tree$tip.label)
  phy <- ape::reorder.phylo(tree, "postorder")   # children before parents
  E   <- phy$edge
  phi <- exp(-lambda * phy$edge.length / td)      # per-edge transmission factor
  nN  <- n + phy$Nnode

  # Up-sweep (post-order): Down[v] = sum over leaves below v of s_i * (product
  # of phi from v down to i).  Leaves: Down = s.  Internal: sum_children phi*Down.
  Down <- numeric(nN)
  Down[seq_len(n)] <- s
  for (k in seq_len(nrow(E))) {
    Down[E[k, 1]] <- Down[E[k, 1]] + phi[k] * Down[E[k, 2]]
  }

  # Down-sweep (pre-order = reverse post-order): Up[v] = signal from everything
  # NOT below v, measured at v.  Up[child] = phi * (Up[parent] + Down[parent] -
  # phi * Down[child]); the bracket is "parent's up-signal + parent's other
  # children" (siblings = Down[parent] - phi*Down[child]).
  Up <- numeric(nN)                                # Up[root] = 0
  for (k in rev(seq_len(nrow(E)))) {
    p <- E[k, 1]; c <- E[k, 2]
    Up[c] <- phi[k] * (Up[p] + Down[p] - phi[k] * Down[c])
  }

  # For a tip f: A_f = w_ff*s_f + sum_{i!=f} w_fi*s_i = s_f + Up[f]  (= Down[f]+Up[f])
  A <- Down[seq_len(n)] + Up[seq_len(n)]
  names(A) <- phy$tip.label
  A[tree$tip.label]                                # return in input tip order
}

SeqDef <- function(tree, df, data.col = 2, invert = TRUE, scale = TRUE, lambda = "auto_max",
                   kernel = c("exponential", "gaussian", "linear")){

  kernel <- match.arg(kernel)
  # Distance-decay kernel on normalized distance x = d / tree_depth (dense path).
  kern <- function(x, lam, type) {
    switch(type,
           exponential = exp(-lam * x),
           gaussian    = exp(-lam * x^2),
           linear      = { z <- 1 - lam * x; z * (z > 0) })  # clamp at 0, preserve matrix dims
  }

  # 1. Convert & Align Data
  df <- as.data.frame(df)
  common_taxa <- intersect(tree$tip.label, df[, 1])
  if (length(common_taxa) < length(tree$tip.label)) {
    tree <- ape::keep.tip(tree, common_taxa)
  }
  df <- df[match(tree$tip.label, df[, 1]), ]
  s_vec <- as.numeric(df[[data.col]])

  # 2. Tree depth (normalizes distances to a relative scale)
  td <- max(ape::branching.times(tree))

  # 3. Distance matrix is only needed for non-exponential kernels OR by_genus
  #    calibration. The exponential path is matrix-free (O(n)).
  is_bygenus  <- (length(lambda) == 1 && is.character(lambda) && lambda == "by_genus")
  need_matrix <- (kernel != "exponential") || is_bygenus
  norm_dists  <- NULL
  if (need_matrix) {
    dist.matrix <- ape::cophenetic.phylo(tree)
    dist.matrix <- dist.matrix[tree$tip.label, tree$tip.label]
    norm_dists  <- dist.matrix / td
  }

  # Raw phylogenetically-weighted availability A for a given lambda.
  avail <- function(lam) {
    if (kernel == "exponential")
      .seqdef_exp_avail(tree, s_vec, lam, td)                       # O(n) traversal
    else
      as.numeric(kern(norm_dists, lam, kernel) %*% s_vec)           # O(n^2) dense
  }

  # 4. Lambda selection
  final_lambda <- 0
  if (length(lambda) == 1 && is.character(lambda) && lambda == "auto_max") {
    message("Optimizing Lambda: Searching for peak variance (Stopping if variance drops >10% below peak)...")
    lambda_seq <- seq(1, 50, 0.1)

    calc_var <- function(x) {
      raw <- avail(x)
      rng <- range(raw)
      if (rng[2] - rng[1] < 1e-9) return(0)
      var(1 - (raw - rng[1]) / (rng[2] - rng[1]))
    }

    best_lambda <- lambda_seq[1]
    max_var <- calc_var(best_lambda)
    for (lam in lambda_seq[-1]) {
      curr_var <- calc_var(lam)
      if (curr_var > max_var) {
        max_var <- curr_var
        best_lambda <- lam
      } else if ((max_var - curr_var) / max_var > 0.10) {
        break
      }
    }
    final_lambda <- best_lambda
    message(sprintf("Selected Lambda: %.1f (Peak Variance: %.5f)", final_lambda, max_var))

  } else if (is_bygenus) {
    genera <- sapply(strsplit(tree$tip.label, "[_ ]"), `[`, 1)
    intra_genus_dists <- c()
    for (g in unique(genera)) {
      tips_in_genus <- tree$tip.label[genera == g]
      if (length(tips_in_genus) > 1) {
        sub_mat <- dist.matrix[tips_in_genus, tips_in_genus]
        d_vals <- sub_mat[lower.tri(sub_mat)]
        d_vals <- d_vals[d_vals > 0]
        intra_genus_dists <- c(intra_genus_dists, d_vals / td)
      }
    }
    if (length(intra_genus_dists) > 0) {
      median_dist <- median(intra_genus_dists)
      final_lambda <- log(2) / median_dist
      message(sprintf("Genus-Scale Lambda: %.2f (Half-life at median genus distance: %.4f)",
                      final_lambda, median_dist))
    } else {
      warning("Biological calibration failed: No genera with >1 species found. Defaulting to Lambda=3. Consider using lambda='auto_max'.")
      final_lambda <- 3
    }

  } else {
    if (!is.numeric(lambda)) stop("Lambda must be 'auto_max', 'by_genus', or numeric.")
    final_lambda <- lambda
  }

  # 5. Final availability + scale + invert
  synscores <- avail(final_lambda)
  names(synscores) <- tree$tip.label

  if (scale) {
    rng <- max(synscores, na.rm = TRUE) - min(synscores, na.rm = TRUE)
    if (rng == 0) synscores[] <- 0 else synscores <- (synscores - min(synscores, na.rm = TRUE)) / rng
  }
  if (invert) synscores <- 1 - synscores

  results <- list(tree = tree, seqdef = synscores, empdata = df[, data.col], lambda = final_lambda)
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
