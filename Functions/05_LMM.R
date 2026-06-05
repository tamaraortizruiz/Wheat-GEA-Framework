# prepare_gemma_phenotype(climate_data, fam_file, phenotype, sample_col, output_file)
# Creates a phenotype file in GEMMA required format from environmental data
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
  
  # If there already is a GEMMA kinship file and overwrite FALSE
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

# climate_data = Data frame with sample IDs and environmental variables
# fam_file = PLINK .fam file for GEMMA analysis
# phenotype = Environmental variable to analyze (e.g. "BIO1")
# sample_col = Column containing sample IDs
# output_file = Output path for the GEMMA phenotype file
# Output: Plain text file containing one phenotype value per sample in FAM order
# Returns: Path to the generated phenotype file
prepare_gemma_phenotype <- function(
    climate_data,
    fam_file,
    phenotype,
    sample_col = "SeedID",
    output_file
) {
  
  # Read PLINK .fam file
  fam <- read.table(fam_file, stringsAsFactors = FALSE)
  colnames(fam) <- c("family.ID", "sample.ID", "paternal.ID",
                     "maternal.ID", "sex", "phenotype")
  
  climate_data[[sample_col]] <- as.character(climate_data[[sample_col]])
  fam$sample.ID <- as.character(fam$sample.ID)
  
  # Match genotype to climate sample order
  climate_ordered <- climate_data[match(fam$sample.ID, climate_data[[sample_col]]), ]
  
  # If genotype samples do not match climate_data
  if (any(is.na(climate_ordered[[sample_col]]))) {
    stop("Some FAM samples are missing from climate_data.")
  }
  
  # Create data frame for selected phenotype
  pheno <- climate_ordered[[phenotype]]
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.table(pheno, output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("GEMMA phenotype file saved to: ", output_file)
  
  return(output_file)
}

# prepare_gemma_covariates(covariates_file, output_file, sample_col, n_pcs)
# Creates a ovariate file in GEMMA required format from PCA results
# covariates_file = .csv file containing PCA covariates
# output_file = Output path for the GEMMA covariate file
# sample_col = Sample ID column to remove before export
# n_pcs = Number of principal components to retain
# Output: Plain text file containing selected PCs in GEMMA format
# Returns: Path to the generated covariate file
prepare_gemma_covariates <- function(
    covariates_file,
    output_file,
    sample_col = "sample.ID",
    n_pcs = 3
) {
  
  # Read PCA covariates file
  covariates <- read.csv(covariates_file, check.names = FALSE)
  # Filter to selected PCs
  pc_cols <- paste0("PC", seq_len(n_pcs))
  covariates <- covariates[, pc_cols, drop = FALSE]
  
  # Match samples to covariates
  if (sample_col %in% colnames(covariates)) {
    covariates <- covariates[, setdiff(colnames(covariates), sample_col), drop = FALSE]
  }
  
  # Write covariates in requires GEMMA format
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  write.table(covariates, output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("GEMMA covariates saved to: ", output_file)
  return(output_file)
}

# run_gemma_lmm(gemma, bfile, phenotype_file, kinship_file, covariates_file,
#               output_prefix, output_dir)
# Runs a single GEMMA linear mixed model (LMM) genome-environment association analysis
# gemma = Path to the GEMMA executable
# bfile = PLINK dataset prefix (.bed/.bim/.fam)
# phenotype_file = GEMMA phenotype file
# kinship_file = GEMMA kinship matrix (.cXX.txt)
# covariates_file = Optional GEMMA covariate file
# output_prefix = Prefix for GEMMA output files
# output_dir = Directory where GEMMA results are written
# Output: GEMMA association results file (.assoc.txt)
# Returns: Path to the association results file
run_gemma_lmm <- function(
    gemma = "gemma",
    bfile,
    phenotype_file,
    kinship_file,
    covariates_file = NULL,
    output_prefix,
    output_dir = "Output/GEA/GEMMA"
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!dir.exists(output_dir)) {
    stop("Could not create GEMMA output directory: ", output_dir)
  }
  
  # Define arguments for GEMMA
  args <- c(
    "-bfile", bfile,
    "-p", phenotype_file,
    "-k", kinship_file,
    "-lmm", "4",
    "-o", output_prefix,
    "-outdir", output_dir
  )
  
  # If covariates file available, includes covariates as argument
  if (!is.null(covariates_file)) {
    args <- c(args, "-c", covariates_file)
  }
  
  message("\nRunning GEMMA:")
  message(gemma, " ", paste(args, collapse = " "))
  
  # Run GEMMA command
  status <- system2(gemma, args = args)
  
  if (status != 0) {
    stop("GEMMA failed. Check GEMMA log files.")
  }
  
  return(file.path(output_dir, paste0(output_prefix, ".assoc.txt")))
}

# run_gemma_lmm_strategies(config, qc_prefix, climate_data, phenotype, kinship_file)
# Runs and evaluates multiple GEMMA population structure correction strategies
# config = Configuration list loaded from YAML
# qc_prefix = PLINK dataset prefix after QC and environmental filtering
# climate_data = Data frame with environmental variables
# phenotype = Environmental variable to be analyzed by GEMMA
# strategies = Vector of number of PCs to use for structure correction (can be 0)
# kinship_file = GEMMA kinship matrix generated from the genotype dataset
# q_threshold = FDR threshold for significance
# Evaluation metrics:
#   λGC (genomic inflation factor)
#   Minimum p-value
#   Number of Bonferroni-significant SNPs
#   Number of FDR-significant SNPs
# Best strategy:
#   Selected automatically as the strategy with λGC closest to 1
# Returns: List of
#   results = Combined results from all strategies
#   evaluation = Evaluation summary table
#   best_strategy = Selected strategy
#   best_results = Results corresponding to the selected strategy
#   individual = Individual results for each strategy
run_gemma_lmm_strategies <- function(
    config,
    qc_prefix,
    climate_data,
    phenotype,
    strategies,
    kinship_file,
    q_threshold = 0.05
) {
  
  # Create lists to save output results
  all_results <- list()
  all_eval <- list()
  
  # For each structure correction strategy
  for (n_pcs in strategies) {
    if (n_pcs == 0) {
      strategy <- "kinship_only"
    } else {
      strategy <- paste0("kinship_", n_pcs, "PCs")
    }
    
    message("\nRunning GEMMA strategy: ", strategy)
    
    # Create GEMMA phenotype file path
    pheno_file <- file.path(config$gemma$output_dir,
                            paste0("pheno_", phenotype, "_", strategy, ".txt"))
    
    # Create GEMMA phenotype file
    gemma_pheno <- prepare_gemma_phenotype(
      climate_data = climate_data,
      fam_file = paste0(qc_prefix, ".fam"),
      phenotype = phenotype,
      sample_col = config$metadata$sample_col,
      output_file = pheno_file
    )
    
    if (n_pcs == 0) {
      gemma_covariates <- NULL
    } else {
      gemma_covariates <- prepare_gemma_covariates(
        covariates_file = config$pca$covariates_file,
        output_file = file.path(config$gemma$output_dir, paste0("covariates_", strategy, ".txt")),
        sample_col = "sample.ID",
        n_pcs = n_pcs
      )
    }
    
    # Run GEMMA LMM and outputs result file
    assoc_file <- run_gemma_lmm(
      gemma = config$gemma$path,
      bfile = qc_prefix,
      phenotype_file = gemma_pheno,
      kinship_file = kinship_file,
      covariates_file = gemma_covariates,
      output_prefix = paste0(phenotype, "_", strategy),
      output_dir = config$gemma$output_dir
    )
    
    # Read GEMMA output file
    result <- read.table(assoc_file, header = TRUE, stringsAsFactors = FALSE)
    
    # Construct structured results table
    result <- result %>%
      mutate(
        method = "GEMMA",
        phenotype = phenotype,
        strategy = strategy,
        chr = chr,
        marker = rs,
        position = ps,
        p_value = p_wald,
        q_value = p.adjust(p_value, method = "fdr"),
        bonferroni_threshold = 0.05 / n(),
        bonferroni_significant = p_value < bonferroni_threshold,
        fdr_significant = q_value < q_threshold
      )
    
    # Genomic inflation factor
    lambda_gc <- median(qchisq(1 - result$p_value, df = 1), na.rm = TRUE) /
      qchisq(0.5, df = 1)
    
    # Construct evaluation data frame
    evaluation <- data.frame(
      method = "GEMMA",
      phenotype = phenotype,
      strategy = strategy,
      n_snps = nrow(result),
      lambda_gc = lambda_gc,
      min_p = min(result$p_value, na.rm = TRUE),
      n_bonferroni = sum(result$bonferroni_significant, na.rm = TRUE),
      n_fdr = sum(result$fdr_significant, na.rm = TRUE),
      bonferroni_threshold = unique(result$bonferroni_threshold)
    )
    
    all_results[[strategy]] <- result
    all_eval[[strategy]] <- evaluation
  }
  
  results_all <- bind_rows(all_results)
  
  evaluation_all <- bind_rows(all_eval) %>%
    mutate(lambda_distance = abs(lambda_gc - 1)) %>%
    arrange(
      lambda_distance,
      desc(n_bonferroni),
      desc(n_fdr)
    )

  # Best strategy selected by:
  # 1. λGC closest to 1
  # 2. highest number of Bonferroni-significant SNPs
  # 3. highest number of FDR-significant SNPs
  best_strategy <- evaluation_all$strategy[1]
  
  return(list(
    results = results_all,
    evaluation = evaluation_all,
    best_strategy = best_strategy,
    best_results = results_all %>% filter(strategy == best_strategy),
    individual = all_results
  ))
}

# run_gemma_lmm_all_variables(config, qc_prefix, climate_data, kinship_file,
#                             phenotypes, strategies)
# Runs GEMMA linear mixed model (LMM) genome-environment association analyses
# across multiple environmental variables and population structure correction strategies.
# config = Configuration list loaded from YAML
# qc_prefix = PLINK dataset prefix after QC and environmental filtering
# climate_data = Data frame with environmental variables
# kinship_file = GEMMA kinship matrix generated from the genotype dataset
# phenotypes = Character vector of environmental variables
# strategies = Numeric vector indicating number of PCs to use for structure correction
# Evaluation metrics:
#   λGC (genomic inflation factor)
#   Minimum p-value
#   Number of Bonferroni-significant SNPs
#   Number of FDR-significant SNPs
# Best strategy:
#   Selected independently for each environmental variable by λGC closest to 1
# Output files:
#   gemma_all_results.csv
#   gemma_all_evaluation.csv
#   gemma_best_strategy_by_var.csv
# Returns: List of
#   by_variable = Nested list with GEMMA outputs for each environmental variable
#   results = Combined results across all variables and strategies
#   evaluation = Combined evaluation table across all variables and strategies
#   best_strategy = Selected strategy
#   best_by_variable = Summary table containing the best strategy for each variable
run_gemma_lmm_all_variables <- function(
    config,
    qc_prefix,
    climate_data,
    kinship_file,
    phenotypes = "BIO1",
    strategies = 0
) {
  
  # Create list for results
  all_gemma <- list()
  
  # For each variable
  for (phenotype in phenotypes) {
    message("\n==============================")
    message("Running GEMMA for: ", phenotype)
    message("==============================")
    
    # Run all strategies
    all_gemma[[phenotype]] <- run_gemma_lmm_strategies(
      config = config,
      qc_prefix = qc_prefix,
      climate_data = climate_data,
      phenotype = phenotype,
      strategies = strategies,
      kinship_file = kinship_file
    )
  }
  
  # Bind all results
  combined_results <- bind_rows(lapply(all_gemma, function(x) x$results))
  combined_evaluation <- bind_rows(lapply(all_gemma, function(x) x$evaluation))
  
  # Write .csv of results
  write_csv(combined_results, file.path(config$gemma$output_dir, "gemma_all_results.csv"))
  write_csv(combined_evaluation, file.path(config$gemma$output_dir, "gemma_all_evaluation.csv"))
  
  # Best strategy for each variable
  best_by_variable <- combined_evaluation %>%
    group_by(phenotype) %>%
    arrange(
      lambda_distance,
      desc(n_bonferroni),
      desc(n_fdr),
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup()
  
  write_csv(best_by_variable, file.path(config$gemma$output_dir, "gemma_best_strategy_by_var.csv"))
  
  return(list(
    by_variable = all_gemma,
    results = combined_results,
    evaluation = combined_evaluation,
    best_by_variable = best_by_variable
  ))
}

# plot_gemma_qq(results, strategy, p_col)
# Creates a quantile-quantile (QQ) plot for a selected GEMMA strategy
# results = Combined results table produced by run_gemma_lmm_strategies()
# phenotype = Environmental variable to plot
# strategy = Strategy to visualize
# p_col = Column containing p-values
# Output: QQ plot showing expected vs observed -log10(p-values)
# Returns: ggplot object
plot_gemma_qq <- function(
    results,
    phenotype,
    strategy,
    p_col = "p_value"
) {
  
  # Extract data for plotting
  plot_data <- results %>%
    filter(
      phenotype == !!phenotype,
      strategy == !!strategy,
      !is.na(.data[[p_col]]),
      .data[[p_col]] > 0,
      .data[[p_col]] <= 1
    )
  
  # Expected vs observed
  observed <- -log10(sort(plot_data[[p_col]]))
  expected <- -log10(ppoints(length(observed)))
  
  # Plot df
  qq_df <- data.frame(
    expected = expected,
    observed = observed
  )
  
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = paste("QQ plot:", phenotype, strategy),
         x = "Expected -log10(p)",
         y = "Observed -log10(p)") +
    theme_classic()
}
