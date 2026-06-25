# prepare_rda_inputs()
# Prepares genotype, marker, sample, environmental, and PCA covariate data for RDA analysis
# geno = Genotype matrix with samples as rows and SNPs as columns
# map = Marker map corresponding to genotype columns
# fam = Sample metadata from PLINK .fam file
# climate_data = Environmental data frame containing climatic variables
# phenotype = Name of climatic variable to analyze
# covariates_file = Optional PCA covariate file
# n_pcs = Number of PCs to include as conditioning variables
# sample_col = Column containing sample IDs
# Output: Filtered and aligned genotype, environmental, and covariate data
# Returns: List containing genotype matrix, environmental variable, optional covariates,
# filtered map, fam, and climate data
prepare_rda_inputs <- function(
    geno,
    map,
    fam,
    climate_data,
    phenotype,
    covariates_file = NULL,
    n_pcs = 0,
    sample_col = "SeedID"
) {
  
  fam$sample.ID <- as.character(fam$sample.ID)
  climate_data[[sample_col]] <- as.character(climate_data[[sample_col]])
  climate_ordered <- climate_data[match(fam$sample.ID, climate_data[[sample_col]]), ]
  
  if (any(is.na(climate_ordered[[sample_col]]))) {
    stop("Some FAM samples are missing from climate_data.")
  }
  
  env <- climate_ordered[[phenotype]]
  keep <- !is.na(env)
  geno_rda <- geno[keep, , drop = FALSE]
  fam_rda <- fam[keep, , drop = FALSE]
  
  # Scale environmental data
  env_rda <- data.frame(env = as.numeric(scale(env[keep])))
  
  # Impute NAs
  geno_rda <- impute_geno_mean(geno_rda)
  # Remove SNPs with 0 variance
  snp_var <- apply(geno_rda, 2, stats::var)
  keep_snps <- !is.na(snp_var) & snp_var > 0
  
  message("Removing ", sum(!keep_snps), " zero-variance SNPs before RDA.")
  
  geno_rda <- geno_rda[, keep_snps, drop = FALSE]
  map_rda <- map[keep_snps, , drop = FALSE]
  
  covariates <- NULL
  
  if (!is.null(covariates_file) && n_pcs > 0) {
    pc_df <- read.csv(covariates_file, check.names = FALSE)
    pc_df$sample.ID <- as.character(pc_df$sample.ID)
    pc_ordered <- pc_df[match(fam_rda$sample.ID, pc_df$sample.ID), ]
    
    if (any(is.na(pc_ordered$sample.ID))) {
      stop("Some RDA samples are missing from PCA covariates.")
    }
    
    pc_cols <- paste0("PC", seq_len(n_pcs))
    covariates <- pc_ordered[, pc_cols, drop = FALSE]
  }
  
  list(
    geno = geno_rda,
    env = env_rda,
    covariates = covariates,
    map = map_rda,
    fam = fam_rda,
    climate = climate_ordered[keep, ]
  )
}

# run_rda_single()
# Runs RDA model
# geno = Prepared genotype matrix
# env = Scaled environmental variable data frame
# map = Marker map matching genotype columns
# covariates = Optional PCA covariates for partial RDA
# phenotype = Name of climatic variable tested
# strategy = Population structure correction strategy name
# output_prefix = Prefix for output files
# output_dir = Directory to write RDA outputs
# q_threshold = FDR threshold for significance
# Output: SNP-level RDA results table and model evaluation table saved as .csv
# Returns: List containing results, evaluation, and fitted RDA model
run_rda_single <- function(
    geno,
    env,
    map,
    covariates = NULL,
    phenotype,
    strategy,
    output_prefix,
    output_dir = "Output/GEA/RDA",
    q_threshold = 0.05
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\nRunning RDA: ", phenotype, " | ", strategy)
  
  # Scale genotype matrix
  geno_scaled <- as.data.frame(scale(geno))
  
  # if no covariates
  if (is.null(covariates)) {
    model_data <- env
    # Run model without covariates
    rda_model <- rda(geno_scaled ~ env, data = model_data)
  } else {
    covariates <- as.data.frame(covariates)
    model_data <- cbind(env, covariates)
    condition_term <- paste(colnames(covariates), collapse = " + ")
    # Define model with covariates
    rda_formula <- as.formula(paste0("geno_scaled ~ env + Condition(", condition_term, ")"))
    rda_model <- rda(rda_formula, data = model_data)
  }
  
  # Scores
  snp_scores <- scores(
    rda_model,
    display = "species",
    choices = 1
  )
  
  # Loadings
  rda_loading <- as.numeric(snp_scores[, 1])
  
  z_scores <- as.numeric(scale(rda_loading))
  p_values <- 2 * pnorm(-abs(z_scores))
  
  result <- map %>%
    mutate(
      method = "RDA",
      phenotype = phenotype,
      strategy = strategy,
      marker = if ("marker.ID" %in% colnames(map)) marker.ID else rownames(map),
      rda_axis = "RDA1",
      rda_loading = rda_loading,
      z_score = z_scores,
      p_value = p_values,
      q_value = p.adjust(p_value, method = "fdr"),
      bonferroni_threshold = 0.05 / n(),
      bonferroni_significant = p_value < bonferroni_threshold,
      fdr_significant = q_value < q_threshold
    )
  
  # Genomic inflation factor
  lambda_gc <- median(qchisq(1 - result$p_value, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  
  # Evaluation data frame
  evaluation <- data.frame(
    method = "RDA",
    phenotype = phenotype,
    strategy = strategy,
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
    model = rda_model
  )
}

# run_rda_strategies()
# Runs RDA across multiple population structure correction strategies
# config = Configuration list loaded from YAML
# geno = Genotype matrix
# map = Marker map
# fam = Sample metadata from PLINK .fam file
# climate_data = Environmental data frame
# phenotype = Climatic variable to analyze
# Output: Combined SNP-level results and evaluation tables for all strategies
# Returns: List containing results, evaluation, best_strategy, best_results and
#          strategy results and models
run_rda_strategies <- function(
    config,
    geno,
    map,
    fam,
    climate_data,
    phenotype
) {
  
  all_results <- list()
  all_eval <- list()
  all_models <- list()
  
  # For each structure correction strategy (number of PCs)
  for (n_pcs in config$rda$pc_strategies) {
    
    if (n_pcs == 0) {
      strategy <- "no_PCs"
    } else {
      strategy <- paste0(n_pcs, "PCs")
    }
    
    rda_input <- prepare_rda_inputs(
      geno = geno,
      map = map,
      fam = fam,
      climate_data = climate_data,
      phenotype = phenotype,
      covariates_file = config$pca$covariates_file,
      n_pcs = n_pcs,
      sample_col = config$metadata$sample_col
    )
    
    rda_run <- run_rda_single(
      geno = rda_input$geno,
      env = rda_input$env,
      map = rda_input$map,
      covariates = rda_input$covariates,
      phenotype = phenotype,
      strategy = strategy,
      output_prefix = paste0(phenotype, "_", strategy),
      output_dir = config$rda$output_dir,
      q_threshold = config$rda$q_threshold
    )
    
    all_results[[strategy]] <- rda_run$results
    all_eval[[strategy]] <- rda_run$evaluation
    all_models[[strategy]] <- rda_run$model
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

# run_rda_all_variables()
# Runs RDA across selected environmental variables and population structure correction strategies
# config = Configuration list loaded from YAML
# geno = Genotype matrix
# map = Marker map
# fam = Sample metadata from PLINK .fam file
# climate_data = Environmental data frame
# phenotypes = Character vector of climatic variables to analyze
# Output: .csv files of all results, all evaluation results and best strategy per variable
# Returns: List containing results by variable, combined results, combined evaluation and best_by_variable table
run_rda_all_variables <- function(
    config,
    geno,
    map,
    fam,
    climate_data,
    phenotypes = config$env$vars
) {
  
  rda_all <- list()
  
  # For each climatic variable
  for (phenotype in phenotypes) {
    
    message("\n==============================")
    message("Running RDA for: ", phenotype)
    message("==============================")
    
    # Run all correction strategies
    rda_all[[phenotype]] <- run_rda_strategies(
      config = config,
      geno = geno,
      map = map,
      fam = fam,
      climate_data = climate_data,
      phenotype = phenotype
    )
  }
  
  combined_results <- bind_rows(lapply(rda_all, function(x) x$results))
  combined_evaluation <- bind_rows(lapply(rda_all, function(x) x$evaluation))
  
  best_by_variable <- combined_evaluation %>%
    group_by(phenotype) %>%
    group_modify(~ select_best_strategy(.x)) %>%
    ungroup()
  
  write_csv(combined_results, file.path(config$rda$output_dir, "rda_all_results.csv"))
  write_csv(combined_evaluation, file.path(config$rda$output_dir, "rda_all_evaluation.csv"))
  write_csv(best_by_variable, file.path(config$rda$output_dir, "rda_best_strategy_by_var.csv"))
  
  list(
    by_variable = rda_all,
    results = combined_results,
    evaluation = combined_evaluation,
    best_by_variable = best_by_variable
  )
}

# plot_rda_qq()
# Creates a quantile-quantile (QQ) plot for all RDA strategies
# results = RDA results table
# phenotype = Environmental variable to plot
# p_col = Column containing p-values
# Output: QQ plots showing expected vs observed -log10(p-values)
# Returns: ggplot object
plot_rda_qq <- function(
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
    geom_point(alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ strategy) +
    labs(
      title = paste("RDA QQ plots:", phenotype),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_classic()
}

# load_rda_results()
# Helper function to load RDA results from files
load_rda_results <- function(config) {
  
  results_file <- file.path(config$rda$output_dir, "rda_all_results.csv")
  evaluation_file <- file.path(config$rda$output_dir, "rda_all_evaluation.csv")
  best_file <- file.path(config$rda$output_dir, "rda_best_strategy_by_var.csv")
  
  if (!all(file.exists(c(results_file, evaluation_file, best_file)))) {
    stop("Saved RDA result files not found.")
  }
  
  list(
    results = read.csv(results_file, check.names = FALSE),
    evaluation = read.csv(evaluation_file, check.names = FALSE),
    best_by_variable = read.csv(best_file, check.names = FALSE)
  )
}