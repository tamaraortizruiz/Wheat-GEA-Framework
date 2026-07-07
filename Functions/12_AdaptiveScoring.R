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
# Returns: Ranked accessions scores for a single phenotype
score_accessions_one_variable <- function(
    geno,
    map,
    fam,
    direction_table,
    phenotype
) {
  
  snp_direction <- direction_table %>%
    filter(
      .data$phenotype == .env$phenotype,
      !is.na(adaptive_direction)
    )
  
  # Empty direction check
  if (nrow(snp_direction) == 0) {
    return(data.frame())
  }
  
  marker_index <- match(snp_direction$marker, map$marker.ID)
  
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
  
  # Recode genotype dosage -> higher values = more adaptive dosage
  adaptive_G <- G
  negative_snps <- snp_direction$adaptive_direction < 0
  
  if (any(negative_snps)) {
    adaptive_G[, negative_snps] <- 2 - adaptive_G[, negative_snps]
  }
  
  n_scored_snps <- rowSums(!is.na(adaptive_G))
  adaptive_dosage_sum <- rowSums(adaptive_G, na.rm = TRUE)
  
  adaptive_score <- ifelse(
    n_scored_snps > 0, # if there are lead adaptive SNPs
    # adaptive_score = adaptive_dosage_sum / max_possible_dosage
    adaptive_dosage_sum / (2 * n_scored_snps),
    NA_real_ # no adaptive SNPs
  )
  
  data.frame(
    sample_id = as.character(fam$sample.ID),
    phenotype = phenotype,
    adaptive_score = adaptive_score,
    n_scored_snps = n_scored_snps,
    adaptive_dosage_sum = adaptive_dosage_sum,
    max_possible_dosage = 2 * n_scored_snps
  ) %>%
    arrange(desc(adaptive_score)) %>%
    mutate(
      rank = row_number(),
      percentile_rank = 100 * percent_rank(adaptive_score)
    )
}


# run_adaptive_germplasm_scoring()
# robustness_results = Consensus robustness results
# gemma_results = GEMMA GEA results
# rda_results = RDA GEA results
# lfmm_results = LFMM GEA results
# qc_prefix = QC filtered PLINK prefix
# metadata = Aligned metadata data frame
# output_dir = Adaptive scoring output directory
# sample_col = Sample identifier column in metadata
# overwrite = defaults to FALSE, uses existing results
# Returns: List of marker effect direction table, direction summary, sample adaptive scores,
#         top 50 adaptive accessions per variable, and adaptive score summary
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
  
  message("\nRunning accession-level adaptive germplasm scoring")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  primary_lead_snps <- robustness_results$primary_lead_snps
  
  if (is.null(primary_lead_snps) || nrow(primary_lead_snps) == 0) {
    stop("No primary lead SNPs found.")
  }
  
  primary_lead_snps <- primary_lead_snps %>%
    filter(selected_primary == TRUE, is_lead == TRUE) %>%
    dplyr::select(phenotype, marker) %>%
    mutate(
      phenotype = as.character(phenotype),
      marker = as.character(marker)
    ) %>%
    distinct()
  
  # SNP effect direction
  direction_table <- infer_adaptive_snp_direction(
    primary_lead_snps = primary_lead_snps,
    gemma_results = gemma_results,
    rda_results = rda_results,
    lfmm_results = lfmm_results
  )
  
  write_csv(direction_table, file.path(output_dir, "adaptive_snp_direction_table.csv"))
  
  # Direction effect summary
  direction_summary <- direction_table %>%
    count(phenotype, direction_source, name = "n_snps")
  
  write_csv(direction_summary, file.path(output_dir, "adaptive_snp_direction_summary.csv"))
  
  qc_obj <- plink_to_bigSNP(
    bed_file = paste0(qc_prefix, ".bed"),
    overwrite = overwrite
  )
  
  geno <- qc_obj$genotypes
  map <- qc_obj$map
  fam <- qc_obj$fam
  
  score_list <- list()
  
  for (selected_phenotype in unique(primary_lead_snps$phenotype)) {
    score_list[[selected_phenotype]] <- score_accessions_one_variable(
      geno = geno,
      map = map,
      fam = fam,
      direction_table = direction_table,
      phenotype = selected_phenotype
    )
  }
  
  adaptive_scores <- bind_rows(score_list)
  
  if (nrow(adaptive_scores) == 0) {
    warning("No adaptive scores were calculated.")
    return(list(
      snp_directions = direction_table,
      direction_summary = direction_summary,
      adaptive_scores = adaptive_scores
    ))
  }
  
  message("\n A high adaptive score means the accession carries many copies of alleles associated with higher values of the environmental variable
  \nA low adaptive score means the accession carries fewer copies of those alleles, or more copies of alleles associated with the opposite end of the environmental gradient")
  
  # Produce adaptive scores for each accessions (ranked order)
  adaptive_scores <- adaptive_scores %>%
    left_join(
      metadata,
      by = setNames(sample_col, "sample_id")
    ) %>%
    group_by(phenotype) %>%
    arrange(desc(adaptive_score), .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      percentile_rank = 100 * percent_rank(adaptive_score)
    ) %>%
    ungroup()
  
  # Extract top 50 ranked accessions (per variable)
  top_50_adaptive_accessions <- adaptive_scores %>%
    group_by(phenotype) %>%
    slice_max(
      order_by = adaptive_score,
      n = 50,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  # Produce adaptive score summary
  adaptive_score_summary <- adaptive_scores %>%
    group_by(phenotype) %>%
    summarise(
      n_accessions = n(),
      n_scored_snps_min = min(n_scored_snps, na.rm = TRUE),
      n_scored_snps_median = median(n_scored_snps, na.rm = TRUE),
      n_scored_snps_max = max(n_scored_snps, na.rm = TRUE),
      adaptive_score_min = min(adaptive_score, na.rm = TRUE),
      adaptive_score_mean = mean(adaptive_score, na.rm = TRUE),
      adaptive_score_median = median(adaptive_score, na.rm = TRUE),
      adaptive_score_max = max(adaptive_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(adaptive_scores, file.path(output_dir, "accession_adaptive_scores.csv"))
  write_csv(top_50_adaptive_accessions, file.path(output_dir, "top_50_adaptive_accessions_by_variable.csv"))
  write_csv(adaptive_score_summary, file.path(output_dir, "adaptive_score_summary_by_variable.csv"))
  
  list(
    snp_directions = direction_table,
    direction_summary = direction_summary,
    adaptive_scores = adaptive_scores,
    top_50_adaptive_accessions = top_50_adaptive_accessions,
    adaptive_score_summary = adaptive_score_summary
  )
}
