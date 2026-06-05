# impute_geno_mean()
# Performs mean imputation of missing genotype values for each SNP
# geno = Genotype matrix or bigsnpr genotype object with samples as rows and SNPs as columns
# Output: Imputed genotype matrix
# Returns: Imputed genotype matrix
impute_geno_mean <- function(geno) {
  G <- as.matrix(geno[])
  
  snp_means <- colMeans(G, na.rm = TRUE)
  
  for (j in seq_len(ncol(G))) {
    missing <- is.na(G[, j])
    if (any(missing)) {
      G[missing, j] <- snp_means[j]
    }
  }
  G
}

# prepare_lfmm_inputs()
# Prepares genotype, marker, sample, and environmental data for LFMM analysis
# geno = Genotype matrix or bigsnpr genotype object
# map = Marker map corresponding to genotype columns
# fam = Sample information table corresponding to genotype rows
# climate_data = Environmental data frame containing sample IDs and climate variables
# phenotype = Environmental variable to test
# sample_col = Column name containing sample IDs in climate_data
# Output: Aligned genotype matrix, environmental matrix, marker map, fam table, and climate data
# Returns: List containing geno, env, map, fam, and climate data
prepare_lfmm_inputs <- function(
    geno,
    map,
    fam,
    climate_data,
    phenotype,
    sample_col = "SeedID"
) {
  # Convert sample IDs to characters for matching
  fam$sample.ID <- as.character(fam$sample.ID)
  climate_data[[sample_col]] <- as.character(climate_data[[sample_col]])
  
  # Reorder to match climate data
  climate_ordered <- climate_data[match(fam$sample.ID, climate_data[[sample_col]]),]
  
  # Ensure all samples are in climate data
  if (any(is.na(climate_ordered[[sample_col]]))) {
    stop("Some FAM samples are missing from climate_data")
  }
  
  # Select climatic variable
  env <- climate_ordered[[phenotype]]
  # Remove NAs
  keep <- !is.na(env)
  geno_lfmm <- geno[keep, , drop = FALSE]
  fam_lfmm <- fam[keep, , drop = FALSE]
  env_lfmm <- scale(as.matrix(env[keep]))
  
  colnames(env_lfmm) <- phenotype
  
  # Impute missing values
  geno_lfmm <- impute_geno_mean(geno_lfmm)
  
  # Remove SNPs with 0 variance
  snp_var <- apply(geno_lfmm, 2, var)
  keep_snps <- !is.na(snp_var) & snp_var > 0
  
  message("Removing ", sum(!keep_snps), " zero-variance SNPs before LFMM.")
  
  geno_lfmm <- geno_lfmm[, keep_snps, drop = FALSE]
  map_lfmm <- map[keep_snps, , drop = FALSE]
  
  # Aligned genotype matrix, environmental matrix, marker map, fam table and climate data
  list(
    geno = geno_lfmm,
    env = env_lfmm,
    map = map_lfmm,
    fam = fam_lfmm,
    climate = climate_ordered[keep, ]
  )
}

# run_lfmm_single()
# Runs LFMM2 GEA model
# geno = LFMM prepared genotype matrix
# env = Scaled environmental matrix for one phenotype
# map = Marker map matching genotype columns
# phenotype = Name of climatic variable tested
# K = Number of latent factors
# output_prefix = Prefix for output files
# output_dir = Directory to write LFMM outputs
# genomic_control = Whether to apply genomic control correction in lfmm2.test()
# q_threshold = FDR threshold for significance
# Output:SNP-level LFMM results table and model evaluation table saved as .csv
# Returns: List containing results, evaluation, model, and test object
run_lfmm_single <- function(
    geno,
    env,
    map,
    phenotype,
    K,
    output_prefix,
    output_dir = "Output/GEA/LFMM",
    genomic_control = TRUE,
    q_threshold = 0.05
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\nRunning LFMM2: ", phenotype, " | K = ", K)
  
  # Fit LFMM model with K latent factors
  mod <- lfmm2(
    input = geno,
    env = env,
    K = K
  )
  
  # Test association
  test <- lfmm2.test(
    object = mod,
    input = geno,
    env = env,
    genomic.control = genomic_control,
    linear = TRUE
  )
  
  # Extract p values and z scores
  p_values <- as.numeric(test$pvalues[, 1])
  z_scores <- as.numeric(test$zscores[, 1])
  
  # Construct structured results table
  result <- map %>%
    mutate(
      method = "LFMM",
      phenotype = phenotype,
      strategy = paste0("K = ", K),
      marker = if ("marker.ID" %in% colnames(map)) marker.ID else rownames(map),
      p_value = p_values,
      z_score = z_scores,
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
    method = "LFMM",
    phenotype = phenotype,
    strategy = paste0("K = ", K),
    K = K,
    n_snps = nrow(result),
    lambda_gc = lambda_gc,
    min_p = min(result$p_value, na.rm = TRUE),
    n_bonferroni = sum(result$bonferroni_significant, na.rm = TRUE),
    n_fdr = sum(result$fdr_significant, na.rm = TRUE),
    bonferroni_threshold = unique(result$bonferroni_threshold)
  )
  
  result_file <- file.path(output_dir, paste0(output_prefix, "_results.csv"))
  eval_file <- file.path(output_dir, paste0(output_prefix, "_evaluation.csv"))
  
  write_csv(result, result_file)
  write_csv(evaluation, eval_file)
  
  list(
    results = result,
    evaluation = evaluation,
    model = mod,
    test = test
  )
}

# run_lfmm_strategies()
# Runs LFMM2 across multiple K values for one phenotype
# config = Configuration list loaded from YAML
# geno = Genotype matrix or bigsnpr genotype object
# map = Marker map corresponding to genotype columns
# fam = Sample information table
# climate_data = Environmental data frame
# phenotype = Environmental variable to test
# Output: Combined LFMM results and evaluation tables for one phenotype
# Returns: List containing results, evaluation, best_strategy, best_results,
#          strategy results, and models
run_lfmm_strategies <- function(
    config,
    geno,
    map,
    fam,
    climate_data,
    phenotype
) {
  # Prepare the aligned genotype and environmental data
  lfmm_input <- prepare_lfmm_inputs(
    geno = geno,
    map = map,
    fam = fam,
    climate_data = climate_data,
    phenotype = phenotype,
    sample_col = config$metadata$sample_col
  )
  
  all_results <- list()
  all_eval <- list()
  all_models <- list()
  
  # for each K in K values to evaluate
  for (K in config$lfmm$K_values) {
    # Name strategy
    strategy <- paste0("K", K)
    
    # Name output prefix
    output_prefix <- paste0(phenotype, "_", strategy)
    
    # Run lfmm for K
    lfmm_run <- run_lfmm_single(
      geno = lfmm_input$geno,
      env = lfmm_input$env,
      map = lfmm_input$map,
      phenotype = phenotype,
      K = K,
      output_prefix = output_prefix,
      output_dir = config$lfmm$output_dir,
      genomic_control = config$lfmm$genomic_control,
      q_threshold = config$lfmm$q_threshold
    )
    
    # Save results for K
    all_results[[strategy]] <- lfmm_run$results
    all_eval[[strategy]] <- lfmm_run$evaluation
    all_models[[strategy]] <- lfmm_run$model
  }
  
  # Bind all results
  results_all <- bind_rows(all_results)
  
  # Bind all result evaluation
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
  
  list(
    results = results_all,
    evaluation = evaluation_all,
    best_strategy = best_strategy,
    best_results = results_all %>% filter(strategy == best_strategy),
    individual = all_results,
    models = all_models
  )
}

# run_lfmm_all_variables()
# Runs LFMM2 for all selected environmental variables
# config = Configuration list loaded from YAML
# geno = Genotype matrix or bigsnpr genotype object
# map = Marker map corresponding to genotype columns
# fam = Sample information table
# climate_data = Environmental data frame
# phenotypes = Vector of environmental variables to analyze
# Output: .csv files of all results, all evaluation results and best strategy per variable
# Returns: List containing results by variable, combined results,
#          combined evaluation, and best_by_variable table
run_lfmm_all_variables <- function(
    config,
    geno,
    map,
    fam,
    climate_data,
    phenotypes = config$climate$vars
) {
  
  lfmm_all <- list()
  
  # For each selected climatic variable
  for (phenotype in phenotypes) {
    
    message("\n==============================")
    message("Running LFMM for: ", phenotype)
    message("==============================")
    
    # Run all correction strategies
    lfmm_all[[phenotype]] <- run_lfmm_strategies(
      config = config,
      geno = geno,
      map = map,
      fam = fam,
      climate_data = climate_data,
      phenotype = phenotype
    )
  }
  
  # Bind all results
  combined_results <- bind_rows(lapply(lfmm_all, function(x) x$results))
  
  # Bind all results evaluation
  combined_evaluation <- bind_rows(lapply(lfmm_all, function(x) x$evaluation))
  combined_evaluation <- combined_evaluation %>%
    mutate(lambda_distance = abs(lambda_gc - 1))
  
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
  
  # Write results
  write_csv(combined_results, file.path(config$lfmm$output_dir, "lfmm_all_results.csv"))
  write_csv(combined_evaluation, file.path(config$lfmm$output_dir, "lfmm_all_evaluation.csv"))
  write_csv(best_by_variable, file.path(config$lfmm$output_dir, "lfmm_best_strategy_by_var.csv"))
  
  list(
    by_variable = lfmm_all,
    results = combined_results,
    evaluation = combined_evaluation,
    best_by_variable = best_by_variable
  )
}

# plot_lfmm_qq()
# Creates a QQ plot for LFMM p-values
# results = LFMM results table
# phenotype = Environmental variable to plot
# strategy = LFMM strategy to plot (by name)
# p_col = Column name containing p-values
# Output: QQ plot comparing expected and observed -log10(p-values)
# Returns: ggplot object
plot_lfmm_qq <- function(
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
  
  # Plot
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = paste("LFMM QQ plot:", phenotype, strategy),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_classic()
}
