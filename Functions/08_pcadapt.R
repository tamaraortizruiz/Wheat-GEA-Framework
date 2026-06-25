# ---- pcadapt ----

# prepare_pcadapt_inputs()
# Prepares PLINK genotype data and marker map for pcadapt analysis
# plink_prefix = PLINK dataset prefix
# Output: pcadapt genotype object and marker map
# Returns: List containing pcadapt genotype object and marker map
prepare_pcadapt_inputs <- function(plink_prefix) {
  
  bed_file <- paste0(plink_prefix, ".bed")
  bim_file <- paste0(plink_prefix, ".bim")
  
  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file)
  }
  if (!file.exists(bim_file)) {
    stop("BIM file not found: ", bim_file)
  }
  
  # bed genotype data to required pcadapt format
  bed <- read.pcadapt(
    bed_file,
    type = "bed"
  )
  
  map <- read.table(bim_file, stringsAsFactors = FALSE)
  colnames(map) <- c(
    "chr",
    "marker",
    "genetic_dist",
    "position",
    "allele1",
    "allele2"
  )
  
  list(
    bed = bed,
    map = map
  )
}

# run_pcadapt_single()
# Runs pcadapt for one K value
# bed = pcadapt genotype object
# map = Aligned marker map
# K = Number of principal components for analysis
# output_prefix = Prefix for output files
# output_dir = Directory to write pcadapt outputs
# q_threshold = FDR threshold for significance
# Output: SNP-level pcadapt results table and evaluation table saved as .csv
# Returns: List containing results, evaluation, and fitted pcadapt model
run_pcadapt_single <- function(
    bed,
    map,
    K,
    output_prefix,
    output_dir = "Output/GEA/pcadapt",
    q_threshold = 0.05
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\nRunning pcadapt: K = ", K)
  
  model <- pcadapt(
    input = bed,
    K = K
  )
  
  p_values <- model$pvalues
  stat <- model$stat
  
  result <- map %>%
    mutate(
      method = "pcadapt",
      phenotype = NA,
      strategy = paste0("K = ", K),
      K = K,
      marker = marker,
      p_value = p_values,
      statistic = stat,
      q_value = p.adjust(p_value, method = "fdr"),
      bonferroni_threshold = 0.05 / n(),
      bonferroni_significant = p_value < bonferroni_threshold,
      fdr_significant = q_value < q_threshold
    )
  
  # Genomic inflation factor
  lambda_gc <- median(qchisq(1 - result$p_value, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  
  # Evaluation data frame
  evaluation <- data.frame(
    method = "pcadapt",
    phenotype = NA,
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
    model = model
  )
}

# run_pcadapt_strategies()
# Runs pcadapt across multiple K values
# config = Configuration list loaded from YAML
# qc_prefix = PLINK dataset prefix after QC and environmental filtering
# Output: SNP-level pcadapt results and evaluation tables for all K values
# Returns: List containing combined results, evaluation, best_strategy,
#          best_results, individual strategy outputs, and fitted models
run_pcadapt_strategies <- function(
    config,
    qc_prefix
) {
  
  pcadapt_input <- prepare_pcadapt_inputs(
    plink_prefix = qc_prefix
  )
  
  all_results <- list()
  all_eval <- list()
  all_models <- list()
  
  # for each K in K values to evaluate
  for (K in config$pcadapt$K_values) {
    strategy <- paste0("K", K)
    
    pcadapt_run <- run_pcadapt_single(
      bed = pcadapt_input$bed,
      map = pcadapt_input$map,
      K = K,
      output_prefix = strategy,
      output_dir = config$pcadapt$output_dir,
      q_threshold = config$pcadapt$q_threshold
    )
    
    all_results[[strategy]] <- pcadapt_run$results
    all_eval[[strategy]] <- pcadapt_run$evaluation
    all_models[[strategy]] <- pcadapt_run$model
  }
  
  results_all <- bind_rows(all_results)
  evaluation_all <- bind_rows(all_eval)
  
  best_strategy <- select_best_strategy(evaluation_all)$strategy[1]
  best_eval <- select_best_strategy(evaluation_all)
  
  write_csv(results_all, file.path(config$pcadapt$output_dir, "pcadapt_all_results.csv"))
  write_csv(evaluation_all, file.path(config$pcadapt$output_dir, "pcadapt_all_evaluation.csv"))
  write_csv(best_eval, file.path(config$pcadapt$output_dir, "pcadapt_best_strategy.csv"))
  
  list(
    results = results_all,
    evaluation = evaluation_all,
    best_strategy = best_strategy,
    best_results = results_all %>% filter(strategy == best_strategy),
    individual = all_results,
    models = all_models,
    best_evaluation = best_eval
  )
}

# plot_pcadapt_qq()
# Creates a quantile-quantile (QQ) plot for pcadapt p-values
# results = pcadapt results table
# p_col = Column containing p-values
# Output: QQ plots showing expected vs observed -log10(p-values)
# Returns: ggplot object
plot_pcadapt_qq <- function(
    results,
    p_col = "p_value"
) {
  
  plot_data <- results %>%
    filter(!is.na(.data[[p_col]]), .data[[p_col]] > 0, .data[[p_col]] <= 1)
  
  qq_df <- plot_data %>%
    group_by(strategy) %>%
    group_modify(~{
      pvals <- sort(.x[[p_col]])
      data.frame(
        expected = -log10(ppoints(length(pvals))),
        observed = -log10(pvals)
      )}) %>%
    ungroup()
  
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ strategy) +
    labs(title = "pcadapt QQ plots",
         x = "Expected -log10(p)",
         y = "Observed -log10(p)") +
    theme_classic()
}

# load_pcadapt_results()
# Helper function to load pcadapt results from files
load_pcadapt_results <- function(config) {
  
  results_file <- file.path(config$pcadapt$output_dir, "pcadapt_all_results.csv")
  evaluation_file <- file.path(config$pcadapt$output_dir, "pcadapt_all_evaluation.csv")
  best_file <- file.path(config$pcadapt$output_dir, "pcadapt_best_strategy.csv")
  
  if (!all(file.exists(c(results_file, evaluation_file, best_file)))) {
    stop("Saved pcadapt result files not found.")
  }
  
  results <- read.csv(results_file, check.names = FALSE)
  best_eval <- read.csv(best_file, check.names = FALSE)
  
  list(
    results = results,
    evaluation = read.csv(evaluation_file, check.names = FALSE),
    best_evaluation = best_eval,
    best_strategy = best_eval$strategy[1],
    best_results = results %>% filter(strategy == best_eval$strategy[1])
  )
}