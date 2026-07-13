# ---- Accession-Level Adaptive Germplasm Scoring Module ----

# make_direction_table()
# Make direction table for one method
# method_results = Method specific results
# value_col = Method specific direction value column
# source_name = Direction source method
# Returns: Formatted method results with direction information
make_direction_table <- function(
    method_results,
    value_col,
    source_name
) {
  
  results <- method_results$results
  
  results %>%
    filter(!is.na(.data[[value_col]])) %>%
    transmute(
      phenotype = as.character(phenotype),
      marker = as.character(marker),
      direction_value = as.numeric(.data[[value_col]]),
      p_value = p_value,
      # define direction
      direction = case_when(
        direction_value > 0 ~ 1,
        direction_value < 0 ~ -1,
        TRUE ~ NA_real_
      ),
      direction_source = source_name
    ) %>%
    filter(!is.na(marker), !is.na(direction)) %>%
    group_by(phenotype, marker) %>%
    # smallest p-value
    arrange(p_value, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()
}


# infer_adaptive_snp_direction()
# Infer adaptive direction for primary lead SNPs
# primary_lead_snps = Selected primary lead SNPs
# gemma_results = GEMMA GEA results
# rda_results = RDA GEA results
# lfmm_results = LFMM GEA results
# Returns: Primary lead SNPs with direction information
infer_adaptive_snp_direction <- function(
    primary_lead_snps,
    gemma_results,
    rda_results,
    lfmm_results
) {
  
  primary_snps <- primary_lead_snps %>%
    dplyr::select(phenotype, marker) %>%
    mutate(
      phenotype = as.character(phenotype),
      marker = as.character(marker)
    ) %>%
    distinct()
  
  # GEMMA -> direction from beta
  gemma_direction <- make_direction_table(
    method_results = gemma_results,
    value_col = "beta",
    source_name = "GEMMA"
  )
  
  # RDA -> direction from oriented_rda_loading
  rda_direction <- make_direction_table(
    method_results = rda_results,
    value_col = "oriented_rda_loading",
    source_name = "RDA"
  )
  
  # LFMM -> direction from z_score
  lfmm_direction <- make_direction_table(
    method_results = lfmm_results,
    value_col = "z_score",
    source_name = "LFMM"
  )
  
  # Add direction to primary lead SNPs
  primary_snps %>%
    left_join(
      gemma_direction %>%
        dplyr::select(
          phenotype,
          marker,
          gemma_direction = direction,
          gemma_value = direction_value
        ),
      by = c("phenotype", "marker")
    ) %>%
    left_join(
      rda_direction %>%
        dplyr::select(
          phenotype,
          marker,
          rda_direction = direction,
          rda_value = direction_value
        ),
      by = c("phenotype", "marker")
    ) %>%
    left_join(
      lfmm_direction %>%
        dplyr::select(
          phenotype,
          marker,
          lfmm_direction = direction,
          lfmm_value = direction_value
        ),
      by = c("phenotype", "marker")
    ) %>%
    mutate(
      adaptive_direction = case_when(
        !is.na(gemma_direction) ~ gemma_direction,
        is.na(gemma_direction) & !is.na(rda_direction) ~ rda_direction,
        is.na(gemma_direction) & is.na(rda_direction) & !is.na(lfmm_direction) ~ lfmm_direction,
        TRUE ~ NA_real_
      ),
      direction_value = case_when(
        !is.na(gemma_direction) ~ gemma_value,
        is.na(gemma_direction) & !is.na(rda_direction) ~ rda_value,
        is.na(gemma_direction) & is.na(rda_direction) & !is.na(lfmm_direction) ~ lfmm_value,
        TRUE ~ NA_real_
      ),
      direction_source = case_when(
        !is.na(gemma_direction) ~ "GEMMA",
        is.na(gemma_direction) & !is.na(rda_direction) ~ "RDA",
        is.na(gemma_direction) & is.na(rda_direction) & !is.na(lfmm_direction) ~ "LFMM",
        TRUE ~ NA_character_
      )
    )
}



# score_accessions_one_variable()
# Score accessions for one environmental variable
# geno = Genotype data matrix
# map = Genotype marker map
# fam = Sample information
# direction_table = Output data frame from infer_adaptive_snp_direction()
# phenotype = Environmental variable
# Returns: Returns directional scores ranging from -1 to +1
score_accessions_one_variable <- function(
    geno,
    map,
    fam,
    direction_table,
    phenotype
) {
  
  # Filter direction table
  snp_direction <- direction_table %>%
    filter(
      .data$phenotype == .env$phenotype,
      !is.na(adaptive_direction)
    )
  
  if (nrow(snp_direction) == 0) {
    return(data.frame())
  }
  
  marker_index <- match(
    snp_direction$marker,
    map$marker.ID
  )
  
  keep <- !is.na(marker_index)
  snp_direction <- snp_direction[keep, ]
  marker_index <- marker_index[keep]
  
  if (length(marker_index) == 0) {
    return(data.frame())
  }
  
  # Genotype dosage matrix filtered to LD-pruned lead SNPs
  G <- as.matrix(geno[, marker_index])
  G <- apply(G, 2, as.numeric)
  
  if (is.null(dim(G))) {
    G <- matrix(G, ncol = 1)
  }
  
  colnames(G) <- snp_direction$marker
  
  # Orient dosage:
  # 2 = two alleles associated with higher environmental values
  # 1 = one allele associated with higher environmental values
  # 0 = no alleles associated with higher environmental values
  oriented_G <- G
  negative_snps <- snp_direction$adaptive_direction < 0
  
  if (any(negative_snps)) {
    oriented_G[, negative_snps] <-
      2 - oriented_G[, negative_snps]
  }
  
  n_total_snps <- ncol(oriented_G)
  n_scored_snps <- rowSums(!is.na(oriented_G))
  adaptive_dosage_sum <- rowSums(oriented_G, na.rm = TRUE)
  
  # Raw score: 0 to 1
  adaptive_score <- ifelse(
    n_scored_snps > 0,
    adaptive_dosage_sum / (2 * n_scored_snps),
    NA_real_
  )
  
  # Centered directional score: -1 to +1
  directional_score <- ifelse(
    !is.na(adaptive_score),
    2 * adaptive_score - 1,
    NA_real_
  )
  
  # Magnitude of directional differentiation: 0 to 1
  absolute_directional_score <- abs(directional_score)
  
  # Label effect direction
  direction_label <- case_when(
    directional_score > 0 ~ "higher_environmental_values",
    directional_score < 0 ~ "lower_environmental_values",
    !is.na(directional_score) ~ "balanced",
    TRUE ~ NA_character_
  )
  
  data.frame(
    sample_id = as.character(fam$sample.ID),
    phenotype = phenotype,
    adaptive_score = adaptive_score,
    directional_score = directional_score,
    absolute_directional_score = absolute_directional_score,
    direction_label = direction_label,
    n_total_snps = n_total_snps,
    n_scored_snps = n_scored_snps,
    scored_snp_fraction = n_scored_snps / n_total_snps,
    adaptive_dosage_sum = adaptive_dosage_sum,
    max_possible_dosage = 2 * n_scored_snps
  ) %>%
    arrange(desc(absolute_directional_score)) %>%
    mutate(
      extremeness_rank = row_number(),
      extremeness_percentile =
        100 * percent_rank(absolute_directional_score)
    )
}

# run_adaptive_germplasm_scoring()
# Run accession-level directional germplasm scoring
# robustness_results = Consensus robustness results
# gemma_results = GEMMA GEA results
# rda_results = RDA GEA results
# lfmm_results = LFMM GEA results
# qc_prefix = QC-filtered PLINK prefix
# metadata = Aligned metadata data frame
# output_dir = Adaptive scoring output directory
# sample_col = Sample identifier column in metadata
# overwrite = Defaults to FALSE; uses existing PLINK conversion
# Returns: SNP direction table, SNP direction summary, accession directional scores,
# top 50 directionally extreme accessions per variable, directional score summary
run_adaptive_germplasm_scoring <- function(
    robustness_results,
    gemma_results,
    rda_results,
    lfmm_results,
    qc_prefix,
    metadata,
    output_dir = "Output/AdaptiveScoring",
    sample_col = "SeedID",
    overwrite = FALSE
) {
  
  message("\nRunning accession-level directional germplasm scoring")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract selected LD-pruned primary lead SNPs
  primary_lead_snps <- robustness_results$primary_lead_snps
  
  if (is.null(primary_lead_snps) || nrow(primary_lead_snps) == 0
  ) {
    stop("No primary lead SNPs found.")
  }
  
  primary_lead_snps <- primary_lead_snps %>%
    filter(
      selected_primary == TRUE,
      is_lead == TRUE
    ) %>%
    dplyr::select(
      phenotype,
      marker
    ) %>%
    mutate(
      phenotype = as.character(phenotype),
      marker = as.character(marker)
    ) %>%
    distinct()
  
  if (nrow(primary_lead_snps) == 0) {
    stop("No SNPs remained after filtering for selected primary lead SNPs.")
  }
  
  # Direction of each selected SNP
  direction_table <- infer_adaptive_snp_direction(
    primary_lead_snps = primary_lead_snps,
    gemma_results = gemma_results,
    rda_results = rda_results,
    lfmm_results = lfmm_results
  )
  
  write_csv(direction_table, file.path(output_dir, "adaptive_snp_direction_table.csv"))
  
  # Summarize direction sources
  direction_summary <- direction_table %>%
    count(
      phenotype,
      direction_source,
      name = "n_snps"
    )
  
  write_csv(direction_summary, file.path(output_dir, "adaptive_snp_direction_summary.csv"))
  
  # Read PLINK genotype data
  qc_obj <- plink_to_bigSNP(
    bed_file = paste0(qc_prefix, ".bed"),
    overwrite = overwrite
  )
  
  geno <- qc_obj$genotypes
  map <- qc_obj$map
  fam <- qc_obj$fam
  
  # Score accessions separately for each environmental variable
  score_list <- list()
  for (
    selected_phenotype in
    unique(primary_lead_snps$phenotype)
  ) {
    score_list[[selected_phenotype]] <-
      score_accessions_one_variable(
        geno = geno,
        map = map,
        fam = fam,
        direction_table = direction_table,
        phenotype = selected_phenotype
      )
  }
  
  adaptive_scores <- bind_rows(score_list)
  
  if (nrow(adaptive_scores) == 0) {
    warning("No directional scores were calculated")
    return(
      list(
        snp_directions = direction_table,
        direction_summary = direction_summary,
        adaptive_scores = adaptive_scores
      )
    )
  }
  
  message(
    "\nDirectional score interpretation:",
    "\n  -1 = strong enrichment for alleles associated with lower",
    " environmental values",
    "\n   0 = balanced or intermediate allele dosage",
    "\n  +1 = strong enrichment for alleles associated with higher",
    " environmental values",
    "\n",
    "\nSign indicates direction.",
    "\nAbsolute value indicates the strength of the directional",
    " genetic profile.",
    "\nAccessions are ranked using the absolute directional score.",
    "\nThese scores represent genotype-environment associations,",
    " not direct measures of fitness or yield."
  )
  
  # Add accession metadata and rank by directional extremeness
  adaptive_scores <- adaptive_scores %>%
    left_join(metadata, by = setNames(sample_col, "sample_id")) %>%
    group_by(phenotype) %>%
    arrange(desc(absolute_directional_score), .by_group = TRUE
    ) %>%
    mutate(
      extremeness_rank = row_number(),
      extremeness_percentile =
        100 * percent_rank(
          absolute_directional_score
        )
    ) %>%
    ungroup()
  
  # Extract the 50 strongest directional profiles
  top_50_extreme_accessions <- adaptive_scores %>%
    group_by(phenotype) %>%
    slice_max(
      order_by = absolute_directional_score,
      n = 50,
      with_ties = FALSE
    ) %>%
    arrange(
      phenotype,
      desc(absolute_directional_score)
    ) %>%
    ungroup()
  
  # Extract high values direction
  top_50_higher_direction <- adaptive_scores %>%
    group_by(phenotype) %>%
    slice_max(
      order_by = directional_score,
      n = 50,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  # Extract low values direction
  top_50_lower_direction <- adaptive_scores %>%
    group_by(phenotype) %>%
    slice_min(
      order_by = directional_score,
      n = 50,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  # directional score summary
  directional_score_summary <- adaptive_scores %>%
    group_by(phenotype) %>%
    summarise(
      n_accessions = n(),
      n_selected_snps = max(n_total_snps, na.rm = TRUE),
      median_scored_snp_fraction = median(scored_snp_fraction, na.rm = TRUE),
      directional_score_min = min(directional_score, na.rm = TRUE),
      directional_score_median = median(directional_score, na.rm = TRUE),
      directional_score_max = max(directional_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save output tables
  write_csv(adaptive_scores, file.path(output_dir, "accession_directional_scores.csv"))
  write_csv(top_50_extreme_accessions,file.path(output_dir, paste0("top_50_directionally_extreme_",
                                                                   "accessions_by_variable.csv")))
  write_csv(top_50_higher_direction, file.path(output_dir, 
                                               "top_50_higher_direction_accessions_by_variable.csv"))
  write_csv(top_50_lower_direction, file.path(output_dir, "top_50_lower_direction_accessions_by_variable.csv"))
  write_csv(directional_score_summary, file.path(output_dir, "directional_score_summary_by_variable.csv"))
  
  # Return all output objects
  list(
    snp_directions = direction_table,
    direction_summary = direction_summary,
    adaptive_scores = adaptive_scores,
    top_50_extreme_accessions = top_50_extreme_accessions,
    top_50_higher_direction = top_50_higher_direction,
    top_50_lower_direction = top_50_lower_direction,
    directional_score_summary = directional_score_summary
  )
}

