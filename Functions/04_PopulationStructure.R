# ---- Population Structure ----

# run_pca_bigsnpr(obj, n_pcs, output_file, covariates_file, ncores)
# Runs PCA on genotype data using bigsnpr
# obj = bigSNP object
# n_pcs = Number of principal components to calculate
# output_file = Optional path to save full PCA results as .rds
# covariates_file = Optional path to save PC covariates as .csv
# ncores = Number of CPU cores to use
# Output: PCA scores, variance explained, and PC covariates
# Returns: List containing PCA object, scores, covariates, variance table,
# map, fam, number of markers, and number of samples
run_pca_bigsnpr <- function(
    obj,
    n_pcs = 10,
    output_file = NULL,
    covariates_file = NULL,
    ncores = 1
) {
  
  G <- obj$genotypes
  fam <- obj$fam
  map <- obj$map
  
  # Impute missing genotypes in bigsnpr format
  G_imp <- snp_fastImputeSimple(G, method = "mean2", ncores = ncores)
  
  # PCA using bigsnpr
  pca <- big_randomSVD(G_imp,
                       k = n_pcs,
                       fun.scaling = snp_scaleBinom(),
                       ncores = ncores)
  
  # Saving scores
  scores <- as.data.frame(pca$u)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores$sample.ID <- fam$sample.ID
  
  # Saving explained variances per each PC
  variance_df <- data.frame(
    PC = paste0("PC", seq_along(pca$d)),
    Eigenvalue = pca$d^2,
    Variance = (pca$d^2 / sum(pca$d^2)) * 100,
    CumulativeVariance = cumsum((pca$d^2 / sum(pca$d^2)) * 100)
  )
  
  # Extracting covariates
  covariates <- scores[, c("sample.ID", paste0("PC", seq_len(n_pcs))), drop = FALSE]
  
  # Saving results
  result <- list(
    pca = pca,
    scores = scores,
    covariates = covariates,
    variance = variance_df,
    map = map,
    fam = fam,
    n_markers_used = ncol(G),
    n_samples = nrow(G)
  )
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, output_file)
    message("PCA result saved to: ", output_file)
  }
  
  if (!is.null(covariates_file)) {
    dir.create(dirname(covariates_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(covariates, covariates_file, row.names = FALSE)
    message("PCA covariates saved to: ", covariates_file)
  }
  
  return(result)
}

# plot_pca(pca_result, metadata, sample_col, color_col, pc_x, pc_y)
# Plots selected principal components from PCA result
# pca_result = Output from run_pca_bigsnpr()
# metadata = Optional metadata data frame to merge with PCA scores for plotting
# sample_col = Metadata column containing sample IDs
# color_col = Optional metadata column used to color points
# pc_x = Principal component for x-axis
# pc_y = Principal component for y-axis
# Output: ggplot PCA scatterplot
# Returns: ggplot object
plot_pca <- function(
    pca_result,
    metadata = NULL,
    sample_col = "SeedID",
    color_col = NULL,
    pc_x = 1,
    pc_y = 2
) {
  
  # Extracting scores and variance
  scores <- pca_result$scores
  variance <- pca_result$variance
  
  # Selecting PCs for each axis
  x_col <- paste0("PC", pc_x)
  y_col <- paste0("PC", pc_y)
  
  if (!is.null(metadata)) {
    # Extracts individuals
    metadata[[sample_col]] <- as.character(metadata[[sample_col]])
    scores$sample.ID <- as.character(scores$sample.ID)
    
    # Joins individuals to scores
    scores <- left_join(scores, metadata,
                        by = setNames(sample_col, "sample.ID"))
  }
  
  # Adding variances to axis labels
  x_lab <- paste0(x_col, " (", round(variance$Variance[pc_x], 2), "%)")
  y_lab <- paste0(y_col, " (", round(variance$Variance[pc_y], 2), "%)")
  
  # if color_col is not NULL, data points colored by color_col
  if (!is.null(color_col)) {
    p <- ggplot(scores, aes(x = .data[[x_col]],
                            y = .data[[y_col]],
                            color = .data[[color_col]])) +
      geom_point(alpha = 0.8) +
      labs(title = "PCA of LD-pruned genotype data",
           x = x_lab,
           y = y_lab,
           color = color_col) +
      theme_classic()
  } else {
    p <- ggplot(scores, aes(x = .data[[x_col]],
                            y = .data[[y_col]])) +
      geom_point(alpha = 0.8) +
      labs(title = "PCA of LD-pruned genotype data",
           x = x_lab,
           y = y_lab,
           color = color_col) +
      theme_classic()
  }
  
  # Return plot
  return(p)
}

# calculate_kinship(geno, method, save, output_file, overwrite)
# Calculates a genetic kinship matrix from genotype data
# geno = Numeric genotype matrix, samples as rows and markers as columns
# method = Kinship method used by statgenGWAS::kinship()
# save = Whether to save the kinship matrix as .csv
# output_file = Path to save or load kinship matrix
# overwrite = If FALSE, reuse existing kinship file.
# Output: Kinship matrix
# Returns: Matrix of pairwise genetic relatedness among samples
calculate_kinship <- function(
    geno,
    method = "vanRaden",
    save = FALSE,
    output_file = NULL,
    overwrite = FALSE
) {
  
  # Overwrite existing kinship according to overwrite argument
  if (!is.null(output_file) && file.exists(output_file) && !overwrite) {
    message("Loading existing kinship matrix")
    return(as.matrix(read.csv(
      output_file,
      row.names = 1,
      check.names = FALSE
    )))
  }
  
  # Check input
  if (!is.matrix(geno)) {
    geno <- as.matrix(geno)
  }
  
  # Mean imputation per SNP
  SNPmeans <- colMeans(geno, na.rm = TRUE)
  
  geno_imputed <- geno
  
  idx <- which(is.na(geno_imputed), arr.ind = TRUE)
  geno_imputed[idx] <- SNPmeans[idx[, 2]]
  
  # Calculate kinship
  kinshipMat <- kinship(geno_imputed, method = method)
  
  # Add sample names
  if (!is.null(rownames(geno))) {
    rownames(kinshipMat) <- rownames(geno)
    colnames(kinshipMat) <- rownames(geno)
  }
  
  # if save, write to output file
  if (save) {
    write.csv(kinshipMat, file = output_file)
    message("Kinship matrix saved to: ", output_file)
  }
  
  return(kinshipMat)
}
