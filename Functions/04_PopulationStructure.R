# ---- Population Structure ----

# run_pca_bigsnpr()
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
    covariates_file = NULL
) {
  
  G <- obj$genotypes
  fam <- obj$fam
  map <- obj$map
  
  # Impute missing genotypes
  G_imp <- snp_fastImputeSimple(G, method = "mean2", ncores = max(1, nb_cores() - 1))
  
  # PCA using bigsnpr
  pca <- big_randomSVD(G_imp,
                       k = n_pcs,
                       fun.scaling = snp_scaleBinom(),
                       ncores = max(1, nb_cores() - 1))
  
  # PC scores
  scores <- as.data.frame(pca$u)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores$sample.ID <- fam$sample.ID
  
  # Variances
  eigenvalues <- pca$d^2
  
  variance_percent <- (eigenvalues / sum(eigenvalues)) * 100
  
  variance_df <- data.frame(
    PC = paste0("PC", seq_along(pca$d)),
    PC_number = seq_along(pca$d),
    Eigenvalue = eigenvalues,
    Variance = variance_percent,
    CumulativeVariance = cumsum(variance_percent)
  )
  
  # Covariates
  covariates <- scores[, c("sample.ID", paste0("PC", seq_len(n_pcs))), drop = FALSE]
  
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

# plot_pca()
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
  
  scores <- pca_result$scores
  variance <- pca_result$variance
  
  x_col <- paste0("PC", pc_x)
  y_col <- paste0("PC", pc_y)
  
  if (!is.null(metadata)) {
    metadata[[sample_col]] <- as.character(metadata[[sample_col]])
    scores$sample.ID <- as.character(scores$sample.ID)
    scores <- left_join(scores, metadata, by = setNames(sample_col, "sample.ID"))
  }
  
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
  
  return(p)
}

# plot_pca_scree()
# Plots percentage of total genotype variance explained by each computed PC
# pca_result = Output from run_pca_bigsnpr()
# n_pcs = Optional number of PCs to display
# Returns: ggplot scree plot
plot_pca_scree <- function(
    pca_result,
    n_pcs = NULL
) {
  variance <- pca_result$variance
  
  if (!is.null(n_pcs)) {
    variance <- variance %>%
      filter(PC_number <= n_pcs)
  }
  
  ggplot(variance, aes(x = PC_number, y = Variance)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = variance$PC_number) +
    labs(
      title = "PCA scree plot",
      x = "Principal component",
      y = "Variance explained (%)"
    ) +
    theme_classic()
}

# plot_pca_cumulative_variance()
# Plots cumulative percentage of total genotype variance explained by computed PCs
# pca_result = Output from run_pca_bigsnpr()
# n_pcs = Optional number of PCs to display
# Returns: ggplot cumulative variance plot
plot_pca_cumulative_variance <- function(
    pca_result,
    n_pcs = NULL
) {
  variance <- pca_result$variance
  
  if (!is.null(n_pcs)) {
    variance <- variance %>%
      filter(PC_number <= n_pcs)
  }
  
  ggplot(variance, aes(x = PC_number, y = CumulativeVariance)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = variance$PC_number) +
    labs(
      title = "PCA cumulative variance",
      x = "Number of principal components",
      y = "Cumulative variance explained (%)"
    ) +
    theme_classic()
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
  
  kinship_file <- file.path(output_dir, paste0(output_prefix, ".cXX.txt"))
  
  # if there already is a GEMMA kinship file and overwrite = FALSE
  if (file.exists(kinship_file) && !overwrite) {
    message("Reusing existing GEMMA kinship file: ", kinship_file)
    return(kinship_file)
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  fam <- read.table(paste0(bfile, ".fam"), stringsAsFactors = FALSE)
  dummy_pheno <- rep(1, nrow(fam))
  dummy_pheno_file <- file.path(output_dir, "dummy_pheno.txt")
  write.table(dummy_pheno, dummy_pheno_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # GEMMA arguments
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

# summarize_gemma_kinship()
# Summarizes GEMMA centered kinship matrix results
# kinship_file = File path for GEMMA kinship matrix
# bfile = Optional PLINK binary prefix used to assign sample IDs
# Output: Data frame with diagonal and off-diagonal kinship statistics
# Returns: Summary statistics data frame
summarize_gemma_kinship <- function(
    kinship_file,
    bfile = NULL
) {
  
  K <- as.matrix(read.table(kinship_file, header = FALSE))
  
  # Add sample names
  if (!is.null(bfile)) {
    fam <- read.table(paste0(bfile, ".fam"), stringsAsFactors = FALSE)
    
    if (nrow(fam) != nrow(K)) {
      stop("Number of samples in .fam (", nrow(fam),
           ") does not match kinship matrix (", nrow(K), ").")
    }
    
    sample_ids <- fam[, 2]
    rownames(K) <- sample_ids
    colnames(K) <- sample_ids
  }
  
  # Separate diagonal and off diagonal values
  diag_values <- diag(K)
  off_diag_values <- K[upper.tri(K)]
  
  # Summary statistics
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

# diagnose_gemma_kinship()
# Performs quality-control diagnostics on a GEMMA centered kinship matrix
# kinship_file = GEMMA .cXX.txt kinship matrix
# bfile = PLINK binary file prefix used to recover sample IDs
# duplicate_threshold = Normalized similarity threshold used to flag near-duplicates
# diag_mad_threshold = Robust MAD threshold used to flag unusual diagonal values
# Output: Diagnostic plots and tables of potentially problematic samples
# Returns: List containing matrix diagnostics, duplicate pairs and flagged samples
diagnose_gemma_kinship <- function(
    kinship_file,
    bfile,
    duplicate_threshold = 0.99,
    diag_mad_threshold = 5,
    plot_diagonal = TRUE,
    plot_off_diagonal = TRUE,
    plot_similarity = TRUE,
    offdiag_plot_max = 500000
) {
  
  K <- as.matrix(read.table(kinship_file, header = FALSE))
  fam <- read.table(paste0(bfile, ".fam"), stringsAsFactors = FALSE)
  
  if (nrow(fam) != nrow(K)) {
    stop("Number of samples in .fam (", nrow(fam),
         ") does not match kinship matrix (", nrow(K), ").")
  }
  
  # Label matrix
  sample_ids <- fam[, 2]
  rownames(K) <- sample_ids
  colnames(K) <- sample_ids
  
  if (any(!is.finite(K))) {
    warning("Kinship matrix contains non-finite values.")
  }
  
  symmetry_error <- max(abs(K - t(K)), na.rm = TRUE)

  # Diagonal outlier detection
  diag_values <- diag(K)
  diag_median <- median(diag_values, na.rm = TRUE)
  diag_mad <- mad(diag_values, na.rm = TRUE)
  
  diag_robust_z <- (diag_values - diag_median) / diag_mad
  
  diagonal <- data.frame(
    sample = sample_ids,
    diagonal = diag_values,
    robust_z = diag_robust_z,
    flagged = abs(diag_robust_z) > diag_mad_threshold
  )

  # Duplicate detectio
  K_normalized <- K / sqrt(outer(diag(K), diag(K)))
  diag(K_normalized) <- NA_real_
  
  # Duplicate pairs
  duplicate_idx <- which(
    K_normalized >= duplicate_threshold & lower.tri(K_normalized),
    arr.ind = TRUE
    )
  
  if (nrow(duplicate_idx) > 0) {
    duplicate_pairs <- data.frame(
      sample1 = rownames(K)[duplicate_idx[, 1]],
      sample2 = colnames(K)[duplicate_idx[, 2]],
      similarity = K_normalized[duplicate_idx]
    )
    duplicate_pairs <- duplicate_pairs[order(-duplicate_pairs$similarity), ]
  } else {
    duplicate_pairs <- data.frame(
      sample1 = character(),
      sample2 = character(),
      similarity = numeric()
    )
  }
  
  # Maximum similarity per each accession
  max_similarity <- apply(K_normalized, 1, max, na.rm = TRUE)
  max_similarity_df <- data.frame(
    sample = names(max_similarity),
    max_similarity = max_similarity,
    duplicate_flag =
      max_similarity >= duplicate_threshold
  )

  off_diag <- K[upper.tri(K)]
  
  # Sample for plotting large matrices
  if (length(off_diag) > offdiag_plot_max) {
    set.seed(123)
    off_diag_plot <- sample(off_diag, offdiag_plot_max)
  } else {
    off_diag_plot <- off_diag
  }

  p_diag <- NULL
  p_offdiag <- NULL
  p_similarity <- NULL
  
  if (plot_diagonal) {
    p_diag <- ggplot(
      diagonal,
      aes(x = diagonal)
    ) +
      geom_histogram(bins = 50) +
      theme_classic() +
      labs(
        title = "GEMMA kinship diagonal distribution",
        x = "Kinship diagonal",
        y = "Number of samples"
      )
  }
  
  if (plot_off_diagonal) {
    p_offdiag <- ggplot(
      data.frame(kinship = off_diag_plot),
      aes(x = kinship)
    ) +
      geom_histogram(bins = 60) +
      geom_vline(
        xintercept = 0,
        linetype = 2
      ) +
      theme_classic() +
      labs(
        title = "GEMMA off-diagonal kinship distribution",
        x = "Pairwise kinship",
        y = "Number of sample pairs"
      )
  }
  
  if (plot_similarity) {
    p_similarity <- ggplot(
      max_similarity_df,
      aes(x = max_similarity)
    ) +
      geom_histogram(bins = 60) +
      geom_vline(
        xintercept = duplicate_threshold,
        linetype = 2
      ) +
      theme_classic() +
      labs(
        title = "Maximum genomic similarity per sample",
        x = "Maximum normalized similarity",
        y = "Number of samples"
      )
  }

  summary <- data.frame(
    metric = c(
      "n_samples",
      "symmetry_error",
      "n_diagonal_flags",
      "n_duplicate_pairs",
      "n_samples_with_duplicate"
    ),
    value = c(
      nrow(K),
      symmetry_error,
      sum(diagonal$flagged),
      nrow(duplicate_pairs),
      sum(max_similarity_df$duplicate_flag)
    )
  )
  
  return(list(
    summary = summary,
    diagonal = diagonal,
    duplicate_pairs = duplicate_pairs,
    max_similarity = max_similarity_df,
    plots = list(
      diagonal = p_diag,
      off_diagonal = p_offdiag,
      max_similarity = p_similarity
    )
  ))
}

# filter_kinship_samples()
# Filters samples from the QC PLINK dataset using kinship diagnostics
# The original QC .bed/.bim/.fam files are overwritten after filtering
# kin_diag = Output from diagnose_gemma_kinship()
# qc_prefix = PLINK QC prefix to filter
# plink = Path to PLINK executable
# remove_duplicates = If TRUE, removes repeated samples from duplicate groups
# remove_diagonal_outliers = If TRUE, removes diagonal kinship outliers
# Output: Prints samples selected for removal and overwrites QC PLINK files
# Returns: Removal table, sample counts, and whether filtering occurred
filter_kinship_samples <- function(
    kin_diag,
    qc_prefix,
    plink = "plink",
    remove_duplicates = TRUE,
    remove_diagonal_outliers = TRUE
) {
  
  fam_file <- paste0(qc_prefix, ".fam")
  
  if (!file.exists(fam_file)) {
    stop("PLINK .fam file not found: ", fam_file)
  }
  
  # Read pre filtered samples
  fam <- read.table(fam_file, stringsAsFactors = FALSE)
  sample_ids <- fam[, 2]
  
  removal_table <- data.frame(
    sample = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  # Diagonal outliers
  if (remove_diagonal_outliers) {
    diagonal_remove <- kin_diag$diagonal[kin_diag$diagonal$flagged, "sample"]
    if (length(diagonal_remove) > 0) {
      removal_table <- rbind(removal_table,
                             data.frame(
                               sample = diagonal_remove,
                               reason = "diagonal_outlier",
                               stringsAsFactors = FALSE
                               )
      )
    }
  }
  
  # Identify duplicates to remove 
  if (remove_duplicates && nrow(kin_diag$duplicate_pairs) > 0) {
    pairs <- kin_diag$duplicate_pairs[, c("sample1", "sample2")]
    # Find connected components in duplicate groups
    remaining <- unique(c(pairs$sample1, pairs$sample2))
    duplicate_remove <- character()
    while (length(remaining) > 0) {
      group <- remaining[1]
      previous_size <- 0
      while (length(group) > previous_size) {
        previous_size <- length(group)
        linked <- unique(c(pairs$sample2[pairs$sample1 %in% group],
                           pairs$sample1[pairs$sample2 %in% group]))
        group <- unique(c(group, linked))
      }
      
      # Order according .fam
      group <- sample_ids[sample_ids %in% group]
      
      # if group member is already being removed, keep another member
      diag_removed <- removal_table$sample[removal_table$reason == "diagonal_outlier"]
      available <- setdiff(group, diag_removed)
      
      if (length(available) > 0) {
        keep_sample <- available[1]
      } else {
        keep_sample <- group[1]
      }
      duplicate_remove <- c(duplicate_remove, setdiff(group, keep_sample))
      remaining <- setdiff(remaining, group)
    }
    
    # Add duplicates to remove table
    if (length(duplicate_remove) > 0) {
      removal_table <- rbind(removal_table,
                             data.frame(
                               sample = duplicate_remove,
                               reason = "near_duplicate",
                               stringsAsFactors = FALSE
                               )
      )
    }
  }

  if (nrow(removal_table) > 0) {
    removal_table <- aggregate(
      reason ~ sample,
      data = removal_table,
      FUN = function(x) paste(unique(x), collapse = ";")
    )
    # Preserve .fam order
    removal_table <- removal_table[
      match(sample_ids[sample_ids %in% removal_table$sample], removal_table$sample),
    ]
    rownames(removal_table) <- NULL
  }
  
  if (nrow(removal_table) > 0) {
    removal_table <- merge(
      removal_table,
      kin_diag$diagonal[ , c("sample", "diagonal", "robust_z")],
      by = "sample",
      all.x = TRUE,
      sort = FALSE
    )
    removal_table <- merge(
      removal_table, kin_diag$max_similarity[ , c("sample", "max_similarity")],
      by = "sample",
      all.x = TRUE,
      sort = FALSE
    )
  }
  
  if (nrow(removal_table) == 0) {
    message("No samples filtered according to criteria.")
    return(list(
      qc_prefix = qc_prefix,
      removal_table = removal_table,
      n_before = nrow(fam),
      n_removed = 0,
      n_after = nrow(fam),
      filtered = FALSE
    ))
  }
  
  print(removal_table, row.names = FALSE)

  remove_fam <- fam[fam[, 2] %in% removal_table$sample, c(1, 2)]
  if (nrow(remove_fam) != nrow(removal_table)) {
    stop("Not all samples selected for removal found in the .fam file.")
  }
  remove_file <- paste0(qc_prefix, "_kinship_remove.txt")
  
  write.table(remove_fam, remove_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save filtering record
  write.csv(removal_table, paste0(qc_prefix, "_kinship_removed_samples.csv"), row.names = FALSE)
 
  # Temporary QC prefix
  temp_prefix <- paste0(qc_prefix, "_kinship_tmp")
  
  args <- c(
    "--bfile", qc_prefix,
    "--allow-extra-chr",
    "--remove", remove_file,
    "--make-bed",
    "--out", temp_prefix
  )
  
  message("\nRunning PLINK kinship filtering:")
  message(plink, " ", paste(args, collapse = " "))
  
  status <- system2(
    plink,
    args = args
  )
  
  if (status != 0) {
    stop("PLINK kinship sample filtering failed.")
  }

  # Confirm output was created
  temp_files <- paste0(temp_prefix, c(".bed", ".bim", ".fam"))
  if (!all(file.exists(temp_files))) {
    stop("Temporary filtered PLINK dataset was not created correctly.")
  }
  
  filtered_fam <- read.table(paste0(temp_prefix, ".fam"), stringsAsFactors = FALSE)
  expected_n <- nrow(fam) - nrow(removal_table)
  
  if (nrow(filtered_fam) != expected_n) {
    stop("Unexpected sample count after filtering. Expected ", expected_n,
         ", found ", nrow(filtered_fam))
  }
 
  extensions <- c(".bed", ".bim", ".fam")
  
  for (ext in extensions) {
    success <- file.copy(
      from = paste0(temp_prefix, ext),
      to = paste0(qc_prefix, ext),
      overwrite = TRUE
    )
    if (!success) {
      stop("Failed to overwrite QC file: ", paste0(qc_prefix, ext))
    }
  }
  
  # Remove temporary PLINK files
  unlink(paste0(temp_prefix, c(".bed", ".bim", ".fam", ".log", ".nosex")))
  
  message("\nKinship filtering complete.")
  message("Samples before: ", nrow(fam))
  message("Samples removed: ", nrow(removal_table))
  message("Samples after: ", nrow(filtered_fam))
  
  return(list(
    qc_prefix = qc_prefix,
    removal_table = removal_table,
    remove_file = remove_file,
    n_before = nrow(fam),
    n_removed = nrow(removal_table),
    n_after = nrow(filtered_fam),
    filtered = TRUE
  ))
}
