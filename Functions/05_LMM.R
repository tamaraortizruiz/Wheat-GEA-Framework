# ---- Linear Mixed Model ----

# prepare_gemma_phenotype(climate_data, fam_file, phenotype, sample_col, output_file)
# Creates a phenotype file in GEMMA required format from environmental data
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
  
  # Import .fam file
  fam <- read.table(fam_file, stringsAsFactors = FALSE)
  # Set column names
  colnames(fam) <- c("family.ID",
                     "sample.ID",
                     "paternal.ID",
                     "maternal.ID",
                     "sex",
                     "phenotype")
  # Format
  climate_data[[sample_col]] <- as.character(climate_data[[sample_col]])
  fam$sample.ID <- as.character(fam$sample.ID)
  
  # Match genotype to climate sample order
  climate_ordered <- climate_data[match(fam$sample.ID, climate_data[[sample_col]]), ]
  
  # if genotype samples do not match climate_data
  if (any(is.na(climate_ordered[[sample_col]]))) {
    stop("Some FAM samples are missing from climate_data.")
  }
  
  # Create data frame for selected phenotype
  pheno <- climate_ordered[[phenotype]]
  
  # Write prepared file
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.table(pheno, output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("GEMMA phenotype file saved to: ", output_file)
  
  return(output_file)
}

# select_best_strategy(evaluation, lambda_min, lambda_max)
# Selects the optimal population structure correction strategy
# Criteria:
# 1. Prioritize strategies with acceptable λGC values (0.8-1.2)
#    and at least one significant SNP (Bonferroni or FDR)
# 2. Among these, select the strategy with λGC closest to 1
# 3. If multiple strategies have similar λGC, prioritize the strategy with the
#    most Bonferroni-significant SNPs, followed by the most FDR-significant SNPs
# 4. If no acceptable strategy contains significant SNPs,
#    select the best-calibrated strategy within the acceptable λGC range
# 5. If no strategy falls within the acceptable λGC range,
#    select the strategy with λGC closest to 1 across all strategies
# evaluation = Evaluation table containing λGC and significance metrics
# lambda_min = Lower acceptable λGC threshold
# lambda_max = Upper acceptable λGC threshold
# Returns: One-row data frame containing the selected strategy
select_best_strategy <- function(evaluation, lambda_min = 0.8, lambda_max = 1.2) {
  
  # Formats evaluation table
  evaluation <- evaluation %>%
    mutate(
      lambda_distance = abs(lambda_gc - 1),
      lambda_acceptable = lambda_gc >= lambda_min & lambda_gc <= lambda_max,
      has_signal = n_bonferroni > 0 | n_fdr > 0
    )
  
  # if lambda_gc falls between acceptable lambda_gc range and there is an adaptive signal
  if (any(evaluation$lambda_acceptable & evaluation$has_signal)) {
    evaluation %>%
      filter(lambda_acceptable, has_signal) %>%
      arrange(
        lambda_distance,
        desc(n_bonferroni),
        desc(n_fdr)
      ) %>%
      slice(1)
  } else if (any(evaluation$lambda_acceptable)) {
    evaluation %>%
      filter(lambda_acceptable) %>%
      arrange(
        lambda_distance,
        desc(n_bonferroni),
        desc(n_fdr)
      ) %>%
      slice(1)
  } else {
    evaluation %>%
      arrange(
        lambda_distance,
        desc(n_bonferroni),
        desc(n_fdr)
      ) %>%
      slice(1)
  }
}

# prepare_gemma_covariates(covariates_file, output_file, sample_col, n_pcs)
# Creates a covariate file in GEMMA required format from PCA results
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
# Runs a single GEMMA LMM
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
  
  # if covariates file available, includes covariates as argument
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
# Runs GEMMA across multiple population structure correction strategies
# config = Configuration list loaded from YAML
# qc_prefix = PLINK dataset prefix after QC and environmental filtering
# climate_data = Data frame with environmental variables
# phenotype = Environmental variable to be analyzed by GEMMA
# strategies = Vector of number of PCs to use for structure correction (can be 0)
# kinship_file = GEMMA kinship matrix generated from the genotype dataset
# q_threshold = FDR threshold for significance
# Returns: List containing results, evaluation, best_strategy, best_results and
#          strategy results
run_gemma_lmm_strategies <- function(
    config,
    qc_prefix,
    climate_data,
    phenotype,
    kinship_file,
    q_threshold = config$gemma$q_threshold
) {
  
  # Create lists to save output results
  all_results <- list()
  all_eval <- list()
  
  # For each structure correction strategy
  for (n_pcs in config$gemma$pc_strategies) {
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
    
    # Prepare covariates
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
    
    # Run GEMMA LMM
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
    lambda_gc <- median(qchisq(1 - result$p_value, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
    
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
    
    # Save to all results
    all_results[[strategy]] <- result
    all_eval[[strategy]] <- evaluation
  }
  
  # Bind all results
  results_all <- bind_rows(all_results)
  
  # Bind all result evaluations
  evaluation_all <- bind_rows(all_eval)
  
  # Best strategy
  best_strategy <- select_best_strategy(evaluation_all)$strategy[1]
  
  return(list(
    results = results_all,
    evaluation = evaluation_all,
    best_strategy = best_strategy,
    best_results = results_all %>% filter(strategy == best_strategy),
    individual = all_results
  ))
}

# run_gemma_lmm_all_variables(config, qc_prefix, climate_data, kinship_file, phenotypes)
# Runs GEMMA LMM across selected environmental variables and population structure correction strategies
# config = Configuration list loaded from YAML
# qc_prefix = PLINK dataset prefix after QC and environmental filtering
# climate_data = Data frame with environmental variables
# kinship_file = GEMMA kinship matrix generated from the genotype dataset
# phenotypes = Character vector of environmental variables=
# Output: .csv files of all results, all evaluation results and best strategy per variable
# Returns: List containing results by variable, combined results,
#          combined evaluation and best_by_variable table
run_gemma_lmm_all_variables <- function(
    config,
    qc_prefix,
    climate_data,
    kinship_file,
    phenotypes = "BIO1"
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
      kinship_file = kinship_file
    )
  }
  
  # Bind all results
  combined_results <- bind_rows(lapply(all_gemma, function(x) x$results))
  combined_evaluation <- bind_rows(lapply(all_gemma, function(x) x$evaluation))
  
  # Best strategy for each variable
  best_by_variable <- combined_evaluation %>%
    group_by(phenotype) %>%
    group_modify(~ select_best_strategy(.x)) %>%
    ungroup()
  
  # Write results
  write_csv(combined_results, file.path(config$gemma$output_dir, "gemma_all_results.csv"))
  write_csv(combined_evaluation, file.path(config$gemma$output_dir, "gemma_all_evaluation.csv"))
  write_csv(best_by_variable, file.path(config$gemma$output_dir, "gemma_best_strategy_by_var.csv"))
  
  return(list(
    by_variable = all_gemma,
    results = combined_results,
    evaluation = combined_evaluation,
    best_by_variable = best_by_variable
  ))
}

# plot_gemma_qq(results, phenotype, p_col)
# Creates a quantile-quantile (QQ) plot for all GEMMA strategies
# results = LMM results table
# phenotype = Environmental variable to plot
# p_col = Column containing p-values
# Output: QQ plots showing expected vs observed -log10(p-values)
# Returns: ggplot object
plot_gemma_qq <- function(
    results,
    phenotype,
    p_col = "p_value"
) {
  
  # Filter plot data for phenotype and valid p-value
  plot_data <- results %>%
    filter(phenotype == !!phenotype,
           !is.na(.data[[p_col]]), .data[[p_col]] > 0, .data[[p_col]] <= 1)
  
  # Create qq data frame grouped by strategy
  qq_df <- plot_data %>%
    group_by(strategy) %>%
    # For each strategy, sorts p-values and creates expected vs observed values
    group_modify(~{
      pvals <- sort(.x[[p_col]])
      data.frame(
        expected = -log10(ppoints(length(pvals))),
        observed = -log10(pvals)
      )}) %>%
    ungroup()
  
  # Plot
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ strategy, scales = "free") +
    labs(title = paste("QQ plots:", phenotype),
         x = "Expected -log10(p)",
         y = "Observed -log10(p)") +
    theme_classic()
}