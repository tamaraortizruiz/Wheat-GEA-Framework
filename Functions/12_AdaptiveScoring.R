# ---- Accession-Level Adaptive Germplasm Scoring Module ----

# make_direction_table()
# Make direction table for one method
# method_results = Method specific results
# value_col = Method specific direction value column
# source_name = Direction source method
# scoring_map = Marker map and allele data
# q_threshold = Configurable q-value threshold
# Returns: Formatted method results with direction information
make_direction_table <- function(
    method_results,
    value_col,
    source_name,
    scoring_map,
    q_threshold
) {
  
  # Selected strategy per environmental variable
  selected <- method_results$best_by_variable %>%
    dplyr::select(phenotype, strategy) %>%
    dplyr::distinct()
  
  if (anyDuplicated(selected$phenotype)) {
    stop("More than one selected strategy per variable: ", source_name)
  }
  
  # Retain significant supporting results
  results <- method_results$results %>%
    dplyr::semi_join(selected, by = c("phenotype", "strategy")) %>%
    dplyr::filter(
      !is.na(marker),
      is.finite(q_value),
      q_value <= q_threshold,
      is.finite(.data[[value_col]]),
      .data[[value_col]] != 0
    )
  
  # GEMMA calls second allele allele0
  # LFMM and RDA retain allele2
  other_col <- if (source_name == "GEMMA") {
    "allele0"
  } else {
    "allele2"
  }
  
  required <- c("allele1", other_col)
  
  if (!all(required %in% names(results))) {
    stop("Missing allele columns in ", source_name, " results")
  }
  
  if (!all(c("marker.ID", "allele1", "allele2") %in% names(scoring_map))) {
    stop("Missing marker or allele columns in scoring_map")
  }
  
  if (anyDuplicated(scoring_map$marker.ID)) {
    stop("Scoring map contains duplicate marker IDs")
  }
  
  # Match result to scoring genotype map
  marker_index <- match(as.character(results$marker), as.character(scoring_map$marker.ID))
  
  if (anyNA(marker_index)) {
    stop("Some ", source_name, " direction markers are absent from the scoring map")
  }
  
  clean_allele <- function(x) {
    toupper(trimws(as.character(x)))
  }
  
  method_a1 <- clean_allele(results$allele1)
  method_a2 <- clean_allele(results[[other_col]])
  
  scoring_a1 <- clean_allele(scoring_map$allele1[marker_index])
  scoring_a2 <- clean_allele(scoring_map$allele2[marker_index])
  
  valid <- !is.na(method_a1) & !is.na(method_a2) &
    !is.na(scoring_a1) & !is.na(scoring_a2) &
    !method_a1 %in% c("", "0") &
    !method_a2 %in% c("", "0") &
    !scoring_a1 %in% c("", "0") &
    !scoring_a2 %in% c("", "0") &
    method_a1 != method_a2 &
    scoring_a1 != scoring_a2
  
  same <- valid &
    method_a1 == scoring_a1 &
    method_a2 == scoring_a2
  
  swapped <- valid &
    method_a1 == scoring_a2 &
    method_a2 == scoring_a1
  
  same[is.na(same)] <- FALSE
  swapped[is.na(swapped)] <- FALSE
  
  if (any(!(same | swapped))) {
    bad_markers <- unique(results$marker[!(same | swapped)])
    stop(
      source_name, " allele mismatch for: ",
      paste(utils::head(bad_markers, 10), collapse = ", ")
    )
  }
  
  # Express every direction relative to scoring allele1
  multiplier <- ifelse(same, 1, -1)
  
  direction_value <- as.numeric(results[[value_col]]) * multiplier
  
  direction_table <- data.frame(
    phenotype = as.character(results$phenotype),
    marker = as.character(results$marker),
    direction_value = direction_value,
    p_value = results$p_value,
    direction = sign(direction_value),
    direction_source = source_name,
    direction_allele = scoring_a1,
    other_allele = scoring_a2,
    allele_alignment = ifelse(same, "matched", "swapped"),
    stringsAsFactors = FALSE
  )
  
  if (anyDuplicated(direction_table[c("phenotype", "marker")])) {
    stop("Duplicate selected results for ", source_name)
  }
  
  message(
    source_name, ": ",
    sum(same), " matched allele pairs; ",
    sum(swapped), " swapped pairs."
  )
  
  direction_table
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
    lfmm_results,
    scoring_map,
    q_threshold
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
    source_name = "GEMMA",
    scoring_map = scoring_map,
    q_threshold = q_threshold
  )
  
  # RDA -> direction from oriented_rda_loading
  rda_direction <- make_direction_table(
    method_results = rda_results,
    value_col = "oriented_rda_loading",
    source_name = "RDA",
    scoring_map = scoring_map,
    q_threshold = q_threshold
  )
  
  # LFMM -> direction from z_score
  lfmm_direction <- make_direction_table(
    method_results = lfmm_results,
    value_col = "z_score",
    source_name = "LFMM",
    scoring_map = scoring_map,
    q_threshold = q_threshold
  )
  
  # Add direction to primary lead SNPs
  final_directions <- primary_snps %>%
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
  
  direction_audit <- dplyr::bind_rows(gemma_direction, rda_direction, lfmm_direction)
  direction_audit <- direction_audit %>%
    dplyr::semi_join(primary_snps, by = c("phenotype", "marker"))
  
  attr(final_directions, "direction_audit") <- direction_audit
  
  return(final_directions)
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
  
  # Numeric genotype dosage matrix filtered to lead SNPs
  oriented_G <- as.matrix(geno[, marker_index])
  
  if (!is.numeric(oriented_G)) {
    stop("Extracted genotype dosage matrix is not numeric.")
  }
  
  colnames(oriented_G) <- snp_direction$marker
  
  # Orient dosage:
  # 2 = two alleles associated with higher environmental values
  # 1 = one allele associated with higher environmental values
  # 0 = no alleles associated with higher environmental values
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
# primary_lead_snps = Final LD-pruned SNP table from the configured primary set
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
    primary_lead_snps,
    gemma_results,
    rda_results,
    lfmm_results,
    qc_prefix,
    metadata,
    output_dir = "Output/AdaptiveScoring",
    sample_col = "SeedID",
    overwrite = FALSE,
    q_threshold = config$consensus$q_threshold
) {
  
  message("\nRunning accession-level directional germplasm scoring")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(primary_lead_snps) || nrow(primary_lead_snps) == 0
  ) {
    stop("No primary lead SNPs found.")
  }
  
  # Read PLINK genotype data
  qc_obj <- plink_to_bigSNP(
    bed_file = paste0(qc_prefix, ".bed"),
    overwrite = overwrite
  )
  
  geno <- qc_obj$genotypes
  map <- qc_obj$map
  fam <- qc_obj$fam

  # Direction of each selected SNP
  direction_table <- infer_adaptive_snp_direction(
    primary_lead_snps = primary_lead_snps,
    gemma_results = gemma_results,
    rda_results = rda_results,
    lfmm_results = lfmm_results,
    scoring_map = map,
    q_threshold = q_threshold
  )
  
  write_csv(direction_table, file.path(output_dir, "adaptive_snp_direction_table.csv"))
  write_csv(attr(direction_table, "direction_audit"), file.path(output_dir, "direction_allele_audit.csv"))
  
  # Summarize direction sources
  direction_summary <- direction_table %>%
    count(
      phenotype,
      direction_source,
      name = "n_snps"
    )
  
  write_csv(direction_summary, file.path(output_dir, "adaptive_snp_direction_summary.csv"))
  
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

# load_adaptive_scoring_results()
# Loads saved accession-level scoring tables for reporting or downstream stages
load_adaptive_scoring_results <- function(config) {
  output_dir <- config$adaptive_scoring$output_dir
  files <- c(
    snp_directions = "adaptive_snp_direction_table.csv",
    direction_summary = "adaptive_snp_direction_summary.csv",
    adaptive_scores = "accession_directional_scores.csv",
    top_50_extreme_accessions = "top_50_directionally_extreme_accessions_by_variable.csv",
    top_50_higher_direction = "top_50_higher_direction_accessions_by_variable.csv",
    top_50_lower_direction = "top_50_lower_direction_accessions_by_variable.csv",
    directional_score_summary = "directional_score_summary_by_variable.csv"
  )
  paths <- file.path(output_dir, unname(files))
  
  if (!all(file.exists(paths))) {
    stop(
      "Saved adaptive-scoring results were not found"
    )
  }
  
  setNames(lapply(paths, read.csv, check.names = FALSE), names(files))
}
