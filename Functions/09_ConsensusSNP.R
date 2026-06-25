# ---- Consensus SNP Sets ----

# standardize_consensus_results()
# Standardizes best-result tables from GEA method
# results = SNP-level result table from one method
# method_name = Name of GEA method
# phenotype_name = Environmental variable to retain
# q_threshold = q_value to utilize as threshold
# Output: Standardized universal SNP-level table
# Returns: Data frame with method, phenotype, marker, chr, position, p_value, q_value
standardize_consensus_results <- function(
    results,
    method_name,
    phenotype_name = NULL,
    q_threshold = 0.1
) {
  
  # Data frame for no significant SNPs
  empty_out <- data.frame(
    method = character(),
    phenotype = character(),
    marker = character(),
    chr = character(),
    position = numeric(),
    p_value = numeric(),
    q_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (is.null(results) || nrow(results) == 0) {
    return(empty_out)
  }
  
  if (!is.null(phenotype_name) && "phenotype" %in% colnames(results)) {
    results <- results %>%
      filter(phenotype == phenotype_name)
  }
  
  # if phenotype has no significant SNPs
  if (nrow(results) == 0) {
    return(empty_out)
  }
  
  # SNPs passing selected q-value threshold
  results <- results %>%
    filter(!is.na(q_value), q_value <= q_threshold)
  
  # if no significant SNPs
  if (nrow(results) == 0) {
    return(empty_out)
  }
  
  if (!"chr" %in% colnames(results)) {
    results <- results %>%
      mutate(chr = NA_character_)
  }
  
  if (!"position" %in% colnames(results)) {
    results <- results %>%
      mutate(position = NA_real_)
  }
  
  results %>%
    transmute(
      method = method_name,
      phenotype = ifelse(is.null(phenotype_name), NA_character_, phenotype_name),
      marker = as.character(marker),
      chr = as.character(chr),
      position = as.numeric(position),
      p_value = as.numeric(p_value),
      q_value = as.numeric(q_value)
    ) %>%
    filter(!is.na(marker), marker != "", !is.na(p_value), p_value > 0, p_value <= 1, !is.na(q_value), q_value <= q_threshold) %>%
    # Keeps all distinct marker method combinations
    distinct(method, marker, .keep_all = TRUE)
}

# prepare_consensus_inputs()
# Prepares consensus results for best LMM, LFMM, RDA and pcadapt results
# gemma_results = GEMMA result object
# lfmm_results = LFMM result object
# rda_results = RDA result object
# pcadapt_results = pcadapt result object
# phenotype = Environmental variable for consensus construction
# q_threshold = q_value to utilize as threshold
# Output: Combined standardized SNP table
# Returns: Data frame of candidate SNPs across methods
prepare_consensus_inputs <- function(
    gemma_results,
    lfmm_results,
    rda_results,
    pcadapt_results,
    phenotype,
    q_threshold = 0.1
) {
  
  # Best GEA method results
  gemma_best <- gemma_results$results %>%
    inner_join(gemma_results$best_by_variable %>% select(phenotype, strategy),
               by = c("phenotype", "strategy"))
  lfmm_best <- lfmm_results$results %>%
    inner_join(lfmm_results$best_by_variable %>% select(phenotype, strategy),
               by = c("phenotype", "strategy"))
  rda_best <- rda_results$results %>%
    inner_join(rda_results$best_by_variable %>% select(phenotype, strategy),
               by = c("phenotype", "strategy"))
  
  # Standardize GEA results
  gemma <- standardize_consensus_results(
    gemma_best,
    "GEMMA",
    phenotype,
    q_threshold = q_threshold
    )
  lfmm  <- standardize_consensus_results(
    lfmm_best,
    "LFMM",
    phenotype,
    q_threshold = q_threshold
    )
  rda   <- standardize_consensus_results(
    rda_best,
    "RDA",
    phenotype,
    q_threshold = q_threshold
    )
  pcadapt <- standardize_consensus_results(
    pcadapt_results$best_results,
    "pcadapt",
    phenotype_name = NULL,
    q_threshold = q_threshold
    ) %>%
    mutate(
      phenotype = phenotype
      )
  
  bind_rows(gemma, lfmm, rda, pcadapt)
}

# build_consensus_categories()
# Builds consensus SNP sets according to method support type
# broad_2methods = q-significant in at least 2 methods (includes all methods)
# env_2methods = q-significant in at least 2 environment-specific GEA methods
# high_confidence = q-significant in at least 2 environmental methods plus pcadapt support
# Output: List of consensus SNP sets
build_consensus_categories <- function(consensus_input) {
  
  env_methods <- c("GEMMA", "LFMM", "RDA")
  
  # if empty consensus input
  empty_out <- data.frame()
  if (is.null(consensus_input) || nrow(consensus_input) == 0) {
    return(list(
      broad_2methods = empty_out,
      env_2methods = empty_out,
      high_confidence = empty_out
    ))
  }
  
  consensus_summary <- consensus_input %>%
    mutate(
      is_env_method = method %in% env_methods,
      is_pcadapt = method == "pcadapt"
    ) %>%
    group_by(phenotype, marker) %>%
    summarise(
      chr = first(na.omit(chr), default = NA_character_),
      position = first(na.omit(position), default = NA_real_),
      n_methods = n_distinct(method),
      methods = paste(sort(unique(method)), collapse = ";"),
      n_env_methods = n_distinct(method[is_env_method]),
      env_methods = paste(sort(unique(method[is_env_method])), collapse = ";"),
      pcadapt_support = any(is_pcadapt),
      min_p = min(p_value, na.rm = TRUE),
      min_q = min(q_value, na.rm = TRUE),
      min_env_p = ifelse(any(is_env_method), min(p_value[is_env_method], na.rm = TRUE), NA_real_),
      min_env_q = ifelse(any(is_env_method), min(q_value[is_env_method], na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      min_env_p = ifelse(is.infinite(min_env_p), NA_real_, min_env_p),
      min_env_q = ifelse(is.infinite(min_env_q), NA_real_, min_env_q)
    )
  
  # Broad 2 method support set
  broad_2methods <- consensus_summary %>%
    filter(n_methods >= 2) %>%
    mutate(consensus_set = "broad_2methods") %>%
    arrange(desc(n_methods), min_q, min_p) %>%
    as.data.frame()
  
  # Environmental method support set
  env_2methods <- consensus_summary %>%
    filter(n_env_methods >= 2) %>%
    mutate(consensus_set = "env_2methods") %>%
    arrange(desc(n_env_methods), min_env_q, min_env_p) %>%
    as.data.frame()
  
  # Environmental + pcadapt method support set
  high_confidence <- consensus_summary %>%
    filter(n_env_methods >= 2, pcadapt_support) %>%
    mutate(consensus_set = "high_confidence") %>%
    arrange(desc(n_env_methods), min_env_q, min_env_p) %>%
    as.data.frame()
  
  list(
    broad_2methods = broad_2methods,
    env_2methods = env_2methods,
    high_confidence = high_confidence
  )
}

# evaluate_consensus_set()
# Summarizes one consensus SNP set
# consensus_df = Consensus SNP table
# consensus_name = Name of consensus set
# phenotype = Environmental variable
# Output: Singlerow summary data frame
evaluate_consensus_set <- function(
    consensus_df,
    consensus_name,
    phenotype
) {

  # if no signal
    if (is.null(consensus_df) || nrow(consensus_df) == 0) {
    return(data.frame(
      phenotype = phenotype,
      consensus_set = consensus_name,
      n_snps = 0,
      mean_methods = NA_real_,
      median_methods = NA_real_,
      max_methods = NA_real_,
      mean_env_methods = NA_real_,
      pcadapt_supported_snps = 0
    ))
  }
  
  # Evaluation data frame
  data.frame(
    phenotype = phenotype,
    consensus_set = consensus_name,
    n_snps = nrow(consensus_df),
    mean_methods = mean(consensus_df$n_methods, na.rm = TRUE),
    median_methods = median(consensus_df$n_methods, na.rm = TRUE),
    max_methods = max(consensus_df$n_methods, na.rm = TRUE),
    mean_env_methods = mean(consensus_df$n_env_methods, na.rm = TRUE),
    pcadapt_supported_snps = sum(consensus_df$pcadapt_support, na.rm = TRUE),
    median_qval = median(consensus_df$min_q , na.rm = TRUE)
  )
}

# run_consensus_single_variable(gemma_results, lfmm_results, rda_results, pcadapt_results
# phenotype, output_dir, q_threshold)
# Builds consensus SNP sets for a single environmental variable
# gemma_results = GEMMA result object
# lfmm_results = LFMM result object
# rda_results = RDA result object
# pcadapt_results = pcadapt result object
# phenotype = Environmental variable
# output_dir = Consensus output directory
# q_threshold = q-value threshold for candidate SNP inclusion
# Output: Consensus SNP sets and evaluation table saved as CSV
run_consensus_single_variable <- function(
    gemma_results,
    lfmm_results,
    rda_results,
    pcadapt_results,
    phenotype,
    output_dir = "Output/ConsensusSNP",
    q_threshold = 0.1
) {
  
  message("\nBuilding consensus sets for: ", phenotype)
  
  phenotype_dir <- file.path(output_dir, phenotype)
  dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)
  
  consensus_input <- prepare_consensus_inputs(
    gemma_results = gemma_results,
    lfmm_results = lfmm_results,
    rda_results = rda_results,
    pcadapt_results = pcadapt_results,
    phenotype = phenotype,
    q_threshold = q_threshold
  )
  
  consensus_sets <- build_consensus_categories(consensus_input)
  
  evaluation <- bind_rows(
    evaluate_consensus_set(
      consensus_sets$broad_2methods,
      "broad_2methods",
      phenotype
    ),
    evaluate_consensus_set(
      consensus_sets$env_2methods,
      "env_2methods",
      phenotype
    ),
    evaluate_consensus_set(
      consensus_sets$high_confidence,
      "high_confidence",
      phenotype
    )
  )
  
  write_csv(consensus_input, file.path(phenotype_dir, "consensus_input.csv"))
  write_csv(consensus_sets$broad_2methods, file.path(phenotype_dir, "broad_2methods.csv"))
  write_csv(consensus_sets$env_2methods, file.path(phenotype_dir, "env_2methods.csv"))
  write_csv(consensus_sets$high_confidence, file.path(phenotype_dir, "high_confidence.csv"))
  write_csv(evaluation, file.path(phenotype_dir, "consensus_evaluation.csv"))
  
  list(
    input = consensus_input,
    broad_2methods = consensus_sets$broad_2methods,
    env_2methods = consensus_sets$env_2methods,
    high_confidence = consensus_sets$high_confidence,
    evaluation = evaluation
  )
}

# run_consensus_all_variables()
# Builds consensus SNP sets for all selected environmental variables
# config = Configuration list loaded from YAML
# gemma_results = GEMMA result object
# lfmm_results = LFMM result object
# rda_results = RDA result object
# pcadapt_results = pcadapt result object
# phenotypes = Vector of environmental variables
# Output: Consensus files per variable and combined evaluation table
run_consensus_all_variables <- function(
    config,
    gemma_results,
    lfmm_results,
    rda_results,
    pcadapt_results,
    phenotypes = config$env$vars
) {
  
  consensus_all <- list()
  
  # For each environmental variable
  for (phenotype in phenotypes) {
    consensus_all[[phenotype]] <- run_consensus_single_variable(
      gemma_results = gemma_results,
      lfmm_results = lfmm_results,
      rda_results = rda_results,
      pcadapt_results = pcadapt_results,
      phenotype = phenotype,
      output_dir = config$consensus$output_dir,
      q_threshold = config$consensus$q_threshold
    )
  }
  
  combined_evaluation <- bind_rows(lapply(consensus_all, function(x) x$evaluation))
  
  write_csv(combined_evaluation, file.path(config$consensus$output_dir, "consensus_all_evaluation.csv"))
  
  list(
    by_variable = consensus_all,
    evaluation = combined_evaluation
  )
}

# read_csv_safe()
# Reads .csv file if it exists and is not empty
# file = .csv file path
# Output: Data frame
read_csv_safe <- function(file) {
  
  # if file empty
  if (!file.exists(file) || file.info(file)$size == 0) {
    return(data.frame())
  }
  
  read_csv(file, show_col_types = FALSE) %>%
    as.data.frame()
}

# consensus_files_exist()
# Checks whether consensus results exist for all phenotypes
# config = YAML configuration list
# phenotypes = List of bioclimatic variables
consensus_files_exist <- function(
    config,
    phenotypes = config$env$vars
) {
  
  output_dir <- config$consensus$output_dir
  required_files <- c(
    file.path(output_dir, "consensus_all_evaluation.csv"),
    unlist(lapply(phenotypes, function(phenotype) {
      file.path(output_dir, phenotype, c("consensus_input.csv", "broad_2methods.csv",
                                         "env_2methods.csv", "high_confidence.csv",
                                         "consensus_evaluation.csv"))
    }))
  )
  
  all(file.exists(required_files))
}

# load_consensus_results()
# Loads saved consensus results from .csv files if available
# config = YAML configuration list
# phenotypes = List of bioclimatic variables
load_consensus_results <- function(
    config,
    phenotypes = config$env$vars
) {
  
  output_dir <- config$consensus$output_dir
  consensus_all <- list()
  
  for (phenotype in phenotypes) {
    phenotype_dir <- file.path(output_dir, phenotype)
    consensus_all[[phenotype]] <- list(
      input = read_csv_safe(file.path(phenotype_dir, "consensus_input.csv")),
      broad_2methods = read_csv_safe(file.path(phenotype_dir, "broad_2methods.csv")),
      env_2methods = read_csv_safe(file.path(phenotype_dir, "env_2methods.csv")),
      high_confidence = read_csv_safe(file.path(phenotype_dir, "high_confidence.csv")),
      evaluation = read_csv_safe(file.path(phenotype_dir, "consensus_evaluation.csv"))
    )
  }
  
  combined_evaluation <- read_csv_safe(file.path(output_dir, "consensus_all_evaluation.csv"))
  
  list(
    by_variable = consensus_all,
    evaluation = combined_evaluation
  )
}

