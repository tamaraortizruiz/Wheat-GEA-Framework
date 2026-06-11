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
  
  # Separate bigSNP object
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
  
  # Save PC scores
  scores <- as.data.frame(pca$u)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores$sample.ID <- fam$sample.ID
  
  # Save explained variances per each PC
  variance_df <- data.frame(
    PC = paste0("PC", seq_along(pca$d)),
    Eigenvalue = pca$d^2,
    Variance = (pca$d^2 / sum(pca$d^2)) * 100,
    CumulativeVariance = cumsum((pca$d^2 / sum(pca$d^2)) * 100)
  )
  
  # Extract covariates
  covariates <- scores[, c("sample.ID", paste0("PC", seq_len(n_pcs))), drop = FALSE]
  
  # Save results
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
  
  # Extract scores and variance
  scores <- pca_result$scores
  variance <- pca_result$variance
  
  # Select PCs for each axis
  x_col <- paste0("PC", pc_x)
  y_col <- paste0("PC", pc_y)
  
  if (!is.null(metadata)) {
    # Extract individuals
    metadata[[sample_col]] <- as.character(metadata[[sample_col]])
    scores$sample.ID <- as.character(scores$sample.ID)
    
    # Join individuals to scores
    scores <- left_join(scores, metadata,
                        by = setNames(sample_col, "sample.ID"))
  }
  
  # Add variances to axis labels
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

# run_gemma_kinship()
# Calculates a kinship matrix using GEMMA directly from PLINK genotype files
# Resulting kinship matrix is saved in GEMMA's native .cXX.txt format
# gemma = Path to the GEMMA executable
# Default = "gemma" assumes GEMMA is available in the system PATH
# bfile = PLINK binary file prefix (.bed/.bim/.fam)
# output_prefix = Prefix assigned to GEMMA output files
# output_dir = Directory to write GEMMA output files
# overwrite = If FALSE, reuse existing kinship file
# Output: GEMMA kinship matrix in centered relatedness format .cXX.txt
# Returns: Output prefix of generated kinship matrix file
run_gemma_kinship <- function(
    gemma = "gemma",
    bfile,
    output_prefix = "gemma_kinship",
    output_dir = "Output/GEA/GEMMA",
    overwrite = FALSE
) {
  
  # Create output file name
  kinship_file <- file.path(output_dir, paste0(output_prefix, ".cXX.txt"))
  
  # if there already is a GEMMA kinship file and overwrite = FALSE
  if (file.exists(kinship_file) && !overwrite) {
    message("Reusing existing GEMMA kinship file: ", kinship_file)
    return(kinship_file)
  }
  
  # Create output
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create dummy pheno file to include samples
  fam <- read.table(paste0(bfile, ".fam"), stringsAsFactors = FALSE)
  dummy_pheno <- rep(1, nrow(fam))
  dummy_pheno_file <- file.path(output_dir, "dummy_pheno.txt")
  write.table(dummy_pheno, dummy_pheno_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Define GEMMA arguments
  args <- c(
    "-bfile", bfile,
    "-p", dummy_pheno_file,
    "-gk", "1",
    "-o", output_prefix,
    "-outdir", output_dir
  )
  
  message("\nRunning GEMMA kinship:")
  message(gemma, " ", paste(args, collapse = " "))
  
  # Run GEMMA command from R using system2()
  status <- system2(gemma, args = args)
  
  if (status != 0) {
    stop("GEMMA kinship calculation failed. Check GEMMA log files.")
  }
  
  if (!file.exists(kinship_file)) {
    stop("Expected GEMMA kinship file not found: ", kinship_file)
  }
  
  message("GEMMA kinship saved to: ", kinship_file)
  
  return(kinship_file)
}

# summarize_gemma_kinship(kinship_file)
# Summarizes GEMMA kinship matrix results
# kinship_file = File path for GEMMA kinship matrix
# Output: A data frame with diagonal and off-diagonal kinship statistics
# Returns: Summary statistics data frame
summarize_gemma_kinship <- function(kinship_file) {
  
  # Read kinship matrix
  K <- as.matrix(read.table(kinship_file, header = FALSE))
  
  # Extract diagonal and off-diagonal values
  diag_values <- diag(K)
  off_diag_values <- K[upper.tri(K)]
  
  # Extract statistics
  stats <- data.frame(
    metric = c(
      "n_samples",
      "diag_min",
      "diag_mean",
      "diag_median",
      "diag_max",
      "offdiag_min",
      "offdiag_mean",
      "offdiag_median",
      "offdiag_max"
    ),
    value = c(
      nrow(K),
      min(diag_values, na.rm = TRUE),
      mean(diag_values, na.rm = TRUE),
      median(diag_values, na.rm = TRUE),
      max(diag_values, na.rm = TRUE),
      min(off_diag_values, na.rm = TRUE),
      mean(off_diag_values, na.rm = TRUE),
      median(off_diag_values, na.rm = TRUE),
      max(off_diag_values, na.rm = TRUE)
    )
  )
  
  return(stats)
}