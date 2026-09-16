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
  
  # Check marker map data
  pick_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  marker_col <- pick_col(results, c("marker", "rs", "marker.ID"))
  chr_col <- pick_col(results, c("chr", "chromosome"))
  pos_col <- pick_col(results, c("position", "ps", "physical.pos"))
  
  if (is.na(marker_col)) {
    stop("No marker column found for ", method_name)
  }
  
  results %>%
    transmute(
      method = method_name,
      phenotype = ifelse(is.null(phenotype_name), NA_character_, phenotype_name),
      marker = as.character(.data[[marker_col]]),
      chr = if (!is.na(chr_col)) as.character(.data[[chr_col]]) else NA_character_,
      position = if (!is.na(pos_col)) as.numeric(.data[[pos_col]]) else NA_real_,
      p_value = as.numeric(p_value),
      q_value = as.numeric(q_value)
    ) %>%
    filter(
      !is.na(marker),
      marker != "",
      !is.na(p_value),
      p_value > 0,
      p_value <= 1,
      !is.na(q_value),
      q_value <= q_threshold
    ) %>%
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
  
  target_phenotype <- as.character(phenotype)
  
  # Best GEA method results
  gemma_best <- gemma_results$results %>%
    inner_join(gemma_results$best_by_variable %>% dplyr::select(phenotype, strategy),
               by = c("phenotype", "strategy"))
  lfmm_best <- lfmm_results$results %>%
    inner_join(lfmm_results$best_by_variable %>% dplyr::select(phenotype, strategy),
               by = c("phenotype", "strategy"))
  rda_best <- rda_results$results %>%
    inner_join(rda_results$best_by_variable %>% dplyr::select(phenotype, strategy),
               by = c("phenotype", "strategy"))
  
  # Standardize GEA results
  gemma <- standardize_consensus_results(
    gemma_best,
    "GEMMA",
    phenotype_name = target_phenotype,
    q_threshold = q_threshold
    )
  lfmm  <- standardize_consensus_results(
    lfmm_best,
    "LFMM",
    phenotype_name = target_phenotype,
    q_threshold = q_threshold
    )
  rda   <- standardize_consensus_results(
    rda_best,
    "RDA",
    phenotype_name = target_phenotype,
    q_threshold = q_threshold
    )
  pcadapt <- standardize_consensus_results(
    pcadapt_results$best_results,
    "pcadapt",
    phenotype_name = NULL,
    q_threshold = q_threshold
  ) %>%
    mutate(
      phenotype = .env$target_phenotype
    )
  
  bind_rows(gemma, lfmm, rda, pcadapt) %>%
    mutate(
      phenotype = as.character(phenotype),
      marker = trimws(as.character(marker)),
      method = trimws(as.character(method))
    )
}

# prepare_snp_overlap_data()
# Converts standardized significant SNP results into method support
# consensus_input = Output from prepare_consensus_inputs()
# Output: One row per phenotype-SNP with TRUE/FALSE membership for each method
# Note: This function does not construct or label consensus sets
prepare_snp_overlap_data <- function(consensus_input) {
  
  method_names <- c("GEMMA", "LFMM", "RDA", "pcadapt")
  
  # Empty output with expected columns
  if (is.null(consensus_input) || nrow(consensus_input) == 0) {
    return(data.frame(
      phenotype = character(),
      marker = character(),
      chr = character(),
      position = numeric(),
      GEMMA = logical(),
      LFMM = logical(),
      RDA = logical(),
      pcadapt = logical(),
      n_methods = integer(),
      intersection = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Keep one chromosome and position per SNP
  marker_coordinates <- consensus_input %>%
    group_by(phenotype, marker) %>%
    summarise(
      chr = dplyr::first(chr[!is.na(chr) & chr != ""], default = NA_character_),
      position = dplyr::first(position[!is.na(position)], default = NA_real_),
      .groups = "drop"
    )
  
  # Convert to binary membership columns
  overlap_data <- consensus_input %>%
    filter(method %in% method_names) %>%
    distinct(phenotype, marker, method) %>%
    mutate(supported = TRUE) %>%
    pivot_wider(
      id_cols = c(phenotype, marker),
      names_from = method,
      values_from = supported,
      values_fill = FALSE
    )
  
  # If method found no significant SNPs
  for (method_name in method_names) {
    if (!method_name %in% names(overlap_data)) {
      overlap_data[[method_name]] <- FALSE
    }
  }
  
  # Ensure logical values
  overlap_data <- overlap_data %>%
    mutate(across(all_of(method_names), ~ tidyr::replace_na(as.logical(.x), FALSE)))
  
  # Method combination for each SNP
  method_matrix <- as.matrix(overlap_data[, method_names, drop = FALSE])
  overlap_data$n_methods <- rowSums(method_matrix)
  overlap_data$intersection <- apply(method_matrix, 1, function(supported) {
      paste(method_names[as.logical(supported)], collapse = " + ")
    }
  )
  
  overlap_data %>%
    left_join(marker_coordinates, by = c("phenotype", "marker")
    ) %>%
    dplyr::select(
      phenotype,
      marker,
      chr,
      position,
      all_of(method_names),
      n_methods,
      intersection
    ) %>%
    arrange(
      desc(n_methods),
      intersection,
      chr,
      position
    ) %>%
    as.data.frame()
}

# summarize_snp_intersections()
# Counts SNPs with each exact method support combination
# overlap_data = Output from prepare_snp_overlap_data()
# Output: One row per exact intersection
summarize_snp_intersections <- function(overlap_data) {
  
  method_names <- c("GEMMA", "LFMM", "RDA", "pcadapt")
  
  if (is.null(overlap_data) || nrow(overlap_data) == 0) {
    return(data.frame())
  }
  
  intersection_summary <- overlap_data %>%
    group_by(
      phenotype,
      across(all_of(method_names))
    ) %>%
    summarise(
      n_snps = n(),
      .groups = "drop"
    )
  
  method_matrix <- as.matrix(intersection_summary[, method_names, drop = FALSE])
  intersection_summary$n_methods <- rowSums(method_matrix)
  intersection_summary$intersection <- apply(method_matrix, 1, function(supported) {
      paste(method_names[as.logical(supported)], collapse = " + ")
    }
  )
  
  intersection_summary %>%
    dplyr::select(
      phenotype,
      intersection,
      all_of(method_names),
      n_methods,
      n_snps
    ) %>%
    arrange(desc(n_methods), desc(n_snps)) %>%
    as.data.frame()
}

# plot_snp_overlap()
# Creates an UpSet plot using ggplot2
plot_snp_overlap <- function(
    overlap_data,
    phenotype
) {
  
  method_names <- c("GEMMA", "LFMM", "RDA", "pcadapt")
  
  if (is.null(overlap_data) || nrow(overlap_data) == 0) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0,
          y = 0,
          label = paste("No significant SNPs were found for", phenotype),
          size = 5
        ) +
        xlim(-1, 1) +
        ylim(-1, 1) +
        theme_void() +
        labs(title = paste("SNP method support overlap:", phenotype))
    )
  }
  
  # Exclude only pcadapt
  plot_data <- overlap_data %>%
    dplyr::filter(GEMMA | LFMM | RDA)
  
  # Count SNPs in every exact method intersection
  intersection_summary <- summarize_snp_intersections(plot_data)
  
  # Order intersections by number of supporting methods and SNP count
  intersection_order <- intersection_summary %>%
    arrange(
      desc(n_methods),
      desc(n_snps),
      intersection
    ) %>%
    pull(intersection)
  
  intersection_summary <- intersection_summary %>%
    mutate(
      intersection = factor(intersection, levels = intersection_order)
    )
  
  # Long format for dot matrix
  matrix_data <- intersection_summary %>%
    dplyr::select(
      intersection,
      all_of(method_names)
    ) %>%
    pivot_longer(
      cols = all_of(method_names),
      names_to = "method",
      values_to = "supported"
    ) %>%
    mutate(
      method = factor(method, levels = rev(method_names)),
      supported = as.logical(supported)
    )
  
  # Data used to connect supported methods within intersections
  connection_data <- matrix_data %>%
    dplyr::filter(supported)
  
  # Intersection bar chart
  intersection_bars <- ggplot(intersection_summary,
                              aes(x = intersection,
                                  y = n_snps)) +
    geom_col(width = 0.7, fill = "#2C7FB8") +
    geom_text(aes(label = n_snps), vjust = -0.3, size = 3.5) +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    labs(
      y = "Number of significant SNPs",
      x = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(t = 5, r = 10, b = 0, l = 10)
    )
  
  # Method-membership matrix
  intersection_matrix <- ggplot(matrix_data, aes(x = intersection, y = method)) +
    # Connect supported methods vertically
    geom_line(data = connection_data, aes(group = intersection), colour = "#252525", linewidth = 0.8) +
    # Unsupported methods in light grey
    geom_point(shape = 16, size = 3.5, colour = "#D9D9D9") +
    # Supported methods in black
    geom_point(data = connection_data, shape = 16, size = 4, colour = "#252525") +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Exact method intersection",
      y = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.4),
      plot.margin = margin(t = 0, r = 10, b = 5, l = 10)
    )
  
  # Stack and align top and bottom of plot
  patchwork::wrap_plots(
    intersection_bars,
    intersection_matrix,
    ncol = 1,
    heights = c(2.3, 1.3)
  ) +
    patchwork::plot_annotation(
      title = paste("SNP-level method support overlap:", phenotype)
      )
}

# run_snp_overlap_single_variable()
# Performs exploratory SNP-level overlap for one environmental variable
run_snp_overlap_single_variable <- function(
    gemma_results,
    lfmm_results,
    rda_results,
    pcadapt_results,
    phenotype,
    output_dir = "Output/ConsensusSNP/SNPOverlap",
    q_threshold = 0.1
) {
  
  phenotype_dir <- file.path(output_dir, phenotype)
  dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)
  
  standardized_input <- prepare_consensus_inputs(
    gemma_results = gemma_results,
    lfmm_results = lfmm_results,
    rda_results = rda_results,
    pcadapt_results = pcadapt_results,
    phenotype = phenotype,
    q_threshold = q_threshold
  )
  
  overlap_data <- prepare_snp_overlap_data(standardized_input)
  
  intersection_summary <- summarize_snp_intersections(overlap_data)
  
  overlap_plot <- plot_snp_overlap(
    overlap_data = overlap_data,
    phenotype = phenotype
  )

  
  # Save one row per SNP with method membership
  write_csv(overlap_data, file.path(phenotype_dir, "snp_method_membership.csv"))
  
  # Save counts for each exact intersection
  write_csv(intersection_summary, file.path(phenotype_dir, "snp_intersection_summary.csv"))
  
  # Save plot
  ggsave(
    filename = file.path(phenotype_dir, "snp_overlap_upset.png"),
    plot = overlap_plot,
    width = 11,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  list(
    overlap_data = overlap_data,
    intersection_summary = intersection_summary,
    overlap_plot = overlap_plot
  )
}

# run_snp_overlap_all_variables()
# Runs SNP-level overlap for all environmental variables
# Output: Per-variable tables, plots and combined intersection summary
run_snp_overlap_all_variables <- function(
    config,
    gemma_results,
    lfmm_results,
    rda_results,
    pcadapt_results,
    phenotypes = config$env$vars,
    output_dir = "Output/ConsensusSNP/SNPOverlap"
) {
  
  overlap_all <- list()
  
  for (phenotype in phenotypes) {
    overlap_all[[phenotype]] <- run_snp_overlap_single_variable(
      gemma_results = gemma_results,
      lfmm_results = lfmm_results,
      rda_results = rda_results,
      pcadapt_results = pcadapt_results,
      phenotype = phenotype,
      output_dir = output_dir,
      q_threshold = config$consensus$q_threshold
    )
  }
  
  combined_intersection_summary <- bind_rows(
    lapply(overlap_all, function(result) result$intersection_summary
    )
  )
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(combined_intersection_summary, file.path(output_dir, "all_variables_intersection_summary.csv"))
  
  list(
    by_variable = overlap_all,
    intersection_summary = combined_intersection_summary
  )
}


# build_consensus_categories()
# Builds tiered consensus SNP sets according to method support type
# env_2methods = q-significant in at least 2 environmental GEA methods
# high_confidence = all 3 environmental methods OR at least 2 environmental
#                   methods + pcadapt
# exploratory_support = exactly 1 environmental method + pcadapt
# Output: Full support summary and fixed consensus evidence sets
build_consensus_categories <- function(consensus_input) {
  
  env_methods <- c("GEMMA", "LFMM", "RDA")
  
  # Return empty outputs when no candidate SNPs are available
  empty_out <- data.frame()
  
  if (is.null(consensus_input) || nrow(consensus_input) == 0) {
    return(list(
      summary = empty_out,
      env_2methods = empty_out,
      high_confidence = empty_out,
      exploratory_support = empty_out
    ))
  }
  
  # Summarize method support for every phenotype-SNP combination
  consensus_summary <- consensus_input %>%
    dplyr::mutate(
      is_env_method = method %in% env_methods,
      is_pcadapt = method == "pcadapt"
    ) %>%
    dplyr::group_by(phenotype, marker) %>%
    dplyr::summarise(
      chr = dplyr::first(chr[!is.na(chr) & chr != ""], default = NA_character_),
      position = dplyr::first(position[!is.na(position)], default = NA_real_),
      n_methods = dplyr::n_distinct(method),
      methods = paste(sort(unique(method)), collapse = ";"),
      n_env_methods = dplyr::n_distinct(method[is_env_method]),
      env_methods = paste(sort(unique(method[is_env_method])), collapse = ";"),
      pcadapt_support = any(is_pcadapt),
      min_p = min(p_value, na.rm = TRUE),
      min_q = min(q_value, na.rm = TRUE),
      min_env_p = ifelse(any(is_env_method), min(p_value[is_env_method], na.rm = TRUE), NA_real_),
      min_env_q = ifelse(any(is_env_method), min(q_value[is_env_method], na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      min_env_p = ifelse(is.infinite(min_env_p), NA_real_, min_env_p),
      min_env_q = ifelse(is.infinite(min_env_q), NA_real_, min_env_q)
    )
  
  # Environmental consensus
  # SNPs supported by at least 2 environmental methods
  env_2methods <- consensus_summary %>%
    dplyr::filter(
      n_env_methods >= 2
    ) %>%
    dplyr::mutate(
      consensus_set = "env_2methods"
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_env_methods),
      dplyr::desc(pcadapt_support),
      min_env_q,
      min_env_p
    ) %>%
    as.data.frame()
  
  # High-confidence subset
  # All 3 environmental methods OR at least 2 environmental methods + pcadapt
  high_confidence <- consensus_summary %>%
    dplyr::filter(
      n_env_methods >= 3 | (n_env_methods >= 2 & pcadapt_support)
    ) %>%
    dplyr::mutate(
      consensus_set = "high_confidence"
    ) %>%
    dplyr::arrange(
      dplyr::desc(n_env_methods),
      min_env_q,
      min_env_p
    ) %>%
    as.data.frame()
  
  # Exploratory complementary support
  # Exactly 1 environmental method + pcadapt
  exploratory_support <- consensus_summary %>%
    dplyr::filter(
      n_env_methods == 1,
      pcadapt_support
    ) %>%
    dplyr::mutate(
      consensus_set = "exploratory_support"
    ) %>%
    dplyr::arrange(
      min_env_q,
      min_env_p
    ) %>%
    as.data.frame()
  
  list(
    summary = as.data.frame(consensus_summary),
    env_2methods = env_2methods,
    high_confidence = high_confidence,
    exploratory_support = exploratory_support
  )
}

# validate_primary_set()
# Validates the configured name of the primary consensus set
# primary_set = "env_2methods" (default) or "high_confidence"
# Returns: Validated primary set name
validate_primary_set <- function(
    primary_set = "env_2methods"
) {

  allowed_primary_sets <- c("env_2methods", "high_confidence")

  if (length(primary_set) != 1 || is.na(primary_set) || !primary_set %in% allowed_primary_sets) {
    stop("Primary set must be one of: ", paste(allowed_primary_sets, collapse = ", "))
  }

  primary_set
}

# evaluate_consensus_set()
# Summarizes one consensus SNP set
# consensus_df = Consensus SNP table
# consensus_name = Name of consensus set
# phenotype = Environmental variable
# Output: Single row summary data frame
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
      pcadapt_supported_snps = 0,
      median_qval = NA_real_
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

# run_consensus_single_variable()
# Builds consensus SNP sets for a single environmental variable
# gemma_results = GEMMA result object
# lfmm_results = LFMM result object
# rda_results = RDA result object
# pcadapt_results = pcadapt result object
# phenotype = Environmental variable
# output_dir = Consensus output directory
# q_threshold = q-value threshold for candidate SNP inclusion
# primary_set = Evidence set used as primary consensus; defaults to env_2methods
# Output: Consensus SNP sets and evaluation table saved as CSV
run_consensus_single_variable <- function(
    gemma_results,
    lfmm_results,
    rda_results,
    pcadapt_results,
    phenotype,
    output_dir = "Output/ConsensusSNP",
    q_threshold = 0.1,
    primary_set = "env_2methods"
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
  primary_set <- validate_primary_set(primary_set)
  
  evaluation <- bind_rows(
    evaluate_consensus_set(
      consensus_sets$env_2methods,
      "env_2methods",
      phenotype
    ),
    evaluate_consensus_set(
      consensus_sets$high_confidence,
      "high_confidence",
      phenotype
    ),
    evaluate_consensus_set(
      consensus_sets$exploratory_support,
      "exploratory_support",
      phenotype
    )
  ) %>%
    dplyr::mutate(
      selected_as_primary = consensus_set == primary_set
    )

  write_csv(consensus_input, file.path(phenotype_dir, "consensus_input.csv"))
  write_csv(consensus_sets$summary, file.path(phenotype_dir, "consensus_support_summary.csv"))
  write_csv(consensus_sets$env_2methods, file.path(phenotype_dir, "env_2methods.csv"))
  write_csv(consensus_sets$high_confidence, file.path(phenotype_dir, "high_confidence.csv"))
  write_csv(consensus_sets$exploratory_support, file.path(phenotype_dir, "exploratory_support.csv"))
  write_csv(evaluation, file.path(phenotype_dir, "consensus_evaluation.csv"))
  
  list(
    input = consensus_input,
    summary = consensus_sets$summary,
    env_2methods = consensus_sets$env_2methods,
    high_confidence = consensus_sets$high_confidence,
    exploratory_support = consensus_sets$exploratory_support,
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

  # Default to the inclusive environmental consensus when YAML omits the key
  primary_set <- config$consensus$primary_set
  
  if (is.null(primary_set) || length(primary_set) == 0 ||
    (length(primary_set) == 1 && (is.na(primary_set) || primary_set == ""))
  ) {
    primary_set <- "env_2methods"
  }
  
  primary_set <- validate_primary_set(primary_set)
  
  # For each environmental variable
  for (phenotype in phenotypes) {
    consensus_all[[phenotype]] <- run_consensus_single_variable(
      gemma_results = gemma_results,
      lfmm_results = lfmm_results,
      rda_results = rda_results,
      pcadapt_results = pcadapt_results,
      phenotype = phenotype,
      output_dir = config$consensus$output_dir,
      q_threshold = config$consensus$q_threshold,
      primary_set = primary_set
    )
  }
  
  combined_evaluation <- bind_rows(lapply(consensus_all, function(x) x$evaluation))
  
  write_csv(combined_evaluation, file.path(config$consensus$output_dir, "consensus_all_evaluation.csv"))
  
  list(
    by_variable = consensus_all,
    primary_set = primary_set,
    evaluation = combined_evaluation
  )
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
      file.path(output_dir, phenotype, c("consensus_input.csv",
                                         "consensus_support_summary.csv",
                                         "env_2methods.csv",
                                         "high_confidence.csv",
                                         "exploratory_support.csv",
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

  primary_set <- config$consensus$primary_set
  if (
    is.null(primary_set) ||
    length(primary_set) == 0 ||
    (length(primary_set) == 1 && (is.na(primary_set) || primary_set == ""))
  ) {
    primary_set <- "env_2methods"
  }
  primary_set <- validate_primary_set(primary_set)
  
  for (phenotype in phenotypes) {
    phenotype_dir <- file.path(output_dir, phenotype)
    consensus_all[[phenotype]] <- list(
      input = read.csv(file.path(phenotype_dir, "consensus_input.csv"), check.names = FALSE),
      summary = read.csv(file.path(phenotype_dir, "consensus_support_summary.csv"), check.names = FALSE),
      env_2methods = read.csv(file.path(phenotype_dir, "env_2methods.csv"), check.names = FALSE),
      high_confidence = read.csv(file.path(phenotype_dir, "high_confidence.csv"), check.names = FALSE),
      exploratory_support = read.csv(file.path(phenotype_dir, "exploratory_support.csv"), check.names = FALSE),
      evaluation = read.csv(file.path(phenotype_dir, "consensus_evaluation.csv"), check.names = FALSE) %>%
        dplyr::mutate(selected_as_primary = consensus_set == primary_set)
    )
  }
  
  combined_evaluation <- read.csv(file.path(output_dir, "consensus_all_evaluation.csv"), check.names = FALSE) %>%
    dplyr::mutate(selected_as_primary = consensus_set == primary_set)
  
  list(
    by_variable = consensus_all,
    primary_set = primary_set,
    evaluation = combined_evaluation
  )
}
