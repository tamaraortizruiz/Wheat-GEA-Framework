# ---- Latent Factor Mixed Model ----
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
  
  return(G)
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
  
  fam$sample.ID <- as.character(fam$sample.ID)
  climate_data[[sample_col]] <- as.character(climate_data[[sample_col]])
  climate_ordered <- climate_data[match(fam$sample.ID, climate_data[[sample_col]]),]
  
  # Missing samples check
  if (any(is.na(climate_ordered[[sample_col]]))) {
    stop("Some FAM samples are missing from climate_data")
  }
  
  env <- climate_ordered[[phenotype]]
  keep <- !is.na(env)
  geno_lfmm <- geno[keep, , drop = FALSE]
  fam_lfmm <- fam[keep, , drop = FALSE]
  env_lfmm <- scale(as.matrix(env[keep]))
  
  colnames(env_lfmm) <- phenotype
  
  # Impute NAs
  geno_lfmm <- impute_geno_mean(geno_lfmm)
  
  # Remove SNPs with 0 variance
  snp_var <- apply(geno_lfmm, 2, stats::var)
  keep_snps <- !is.na(snp_var) & snp_var > 0
  
  message("Removing ", sum(!keep_snps), " zero-variance SNPs before LFMM.")
  
  geno_lfmm <- geno_lfmm[, keep_snps, drop = FALSE]
  map_lfmm <- map[keep_snps, , drop = FALSE]
  
  list(
    geno = geno_lfmm,
    env = env_lfmm,
    map = map_lfmm,
    fam = fam_lfmm,
    climate = climate_ordered[keep, ]
  )
}

# run_lfmm_single()
# Runs a single LFMM2
# geno = LFMM prepared genotype matrix
# env = Scaled environmental matrix for one phenotype
# map = Marker map matching genotype columns
# phenotype = Name of climatic variable tested
# K = Number of latent factors
# output_prefix = Prefix for output files
# output_dir = Directory to write LFMM outputs
# genomic_control = Whether to apply genomic control correction in lfmm2.test()
# q_threshold = FDR threshold for significance
# Output: SNP-level LFMM results table and model evaluation table saved as .csv
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
    q_threshold = config$lfmm$q_threshold
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\nRunning LFMM2: ", phenotype, " | K = ", K)
  
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
  
  p_values <- as.numeric(test$pvalues[, 1])
  z_scores <- as.numeric(test$zscores[, 1])
  
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
  lambda_gc <- median(qchisq(1 - result$p_value, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  
  # Evaluation data frame
  evaluation <- data.frame(
    method = "LFMM",
    phenotype = phenotype,
    strategy = paste0("K = ", K),
    K = K,
    n_snps = nrow(result),
    lambda_gc = lambda_gc,
    min_p = min(result$p_value, na.rm = TRUE),
    min_q = min(result$q_value, na.rm = TRUE),
    n_bonferroni = sum(result$bonferroni_significant, na.rm = TRUE),
    n_fdr = sum(result$fdr_significant, na.rm = TRUE),
    bonferroni_threshold = unique(result$bonferroni_threshold)
  )
  
  write_csv(result, file.path(output_dir, paste0(output_prefix, "_results.csv")))
  write_csv(evaluation, file.path(output_dir, paste0(output_prefix, "_evaluation.csv")))
  
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
#          strategy results and models
run_lfmm_strategies <- function(
    config,
    geno,
    map,
    fam,
    climate_data,
    phenotype
) {
  
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
  
  # for each K
  for (K in config$lfmm$K_values) {
    strategy <- paste0("K", K)
    output_prefix <- paste0(phenotype, "_", strategy)
    
    # Run LFMM for K
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
    
    all_results[[strategy]] <- lfmm_run$results
    all_eval[[strategy]] <- lfmm_run$evaluation
    all_models[[strategy]] <- lfmm_run$model
  }
  
  results_all <- bind_rows(all_results)
  evaluation_all <- bind_rows(all_eval)
  
  best_strategy <- select_best_strategy(evaluation_all)$strategy[1]
  
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
# Runs LFMM across selected environmental variables and population structure correction strategies
# config = Configuration list loaded from YAML
# geno = Genotype matrix or bigsnpr genotype object
# map = Marker map corresponding to genotype columns
# fam = Sample information table
# climate_data = Environmental data frame
# phenotypes = Vector of environmental variables to analyze
# Output: .csv files of all results, all evaluation results and best strategy per variable
# Returns: List containing results by variable, combined results, combined evaluation and best_by_variable table
run_lfmm_all_variables <- function(
    config,
    geno,
    map,
    fam,
    climate_data,
    phenotypes = config$climate$vars
) {
  
  lfmm_all <- list()
  
  for (phenotype in phenotypes) {
    
    message("\n==============================")
    message("Running LFMM for: ", phenotype)
    message("==============================")
    
    lfmm_all[[phenotype]] <- run_lfmm_strategies(
      config = config,
      geno = geno,
      map = map,
      fam = fam,
      climate_data = climate_data,
      phenotype = phenotype
    )
  }
  
  combined_results <- bind_rows(lapply(lfmm_all, function(x) x$results))
  combined_evaluation <- bind_rows(lapply(lfmm_all, function(x) x$evaluation))
  
  best_by_variable <- combined_evaluation %>%
    group_by(phenotype) %>%
    group_modify(~ select_best_strategy(.x)) %>%
    ungroup()
  
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
# Creates a quantile-quantile (QQ) plot for all LFMM strategies
# results = LFMM results table
# phenotype = Environmental variable to plot
# p_col = Column containing p-values
# Output: QQ plots comparing expected and observed -log10(p-values)
# Returns: ggplot object
plot_lfmm_qq <- function(
    results,
    phenotype,
    p_col = "p_value"
) {
  
  plot_data <- results %>%
    filter(phenotype == !!phenotype,
           !is.na(.data[[p_col]]), .data[[p_col]] > 0, .data[[p_col]] <= 1)
  
  qq_df <- plot_data %>%
    group_by(strategy) %>%
    group_modify(~{
      pvals <- sort(.x[[p_col]])
      data.frame(
        expected = -log10(ppoints(length(pvals))),
        observed = -log10(pvals)
      )
    }) %>%
    ungroup()
  
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    facet_wrap(~ strategy) +
    labs(
      title = paste("LFMM QQ plots:", phenotype),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_classic()
}

# load_lfmm_results()
# Helper function to load LFMM results from files
load_lfmm_results <- function(config) {
  
  results_file <- file.path(config$lfmm$output_dir, "lfmm_all_results.csv")
  evaluation_file <- file.path(config$lfmm$output_dir, "lfmm_all_evaluation.csv")
  best_file <- file.path(config$lfmm$output_dir, "lfmm_best_strategy_by_var.csv")
  
  if (!all(file.exists(c(results_file, evaluation_file, best_file)))) {
    stop("Saved LFMM result files not found.")
  }
  
  list(
    results = read.csv(results_file, check.names = FALSE),
    evaluation = read.csv(evaluation_file, check.names = FALSE),
    best_by_variable = read.csv(best_file, check.names = FALSE)
  )
}
