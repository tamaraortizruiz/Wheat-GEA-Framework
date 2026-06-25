# ---- Consensus Robustness Evaluation Module ----

# summarize_robustness_all_sets(ld_results)
# Builds robustness metrics from LD pruned consensus results
# ld_results = Output from LD pruning module
summarize_robustness_all_sets <- function(ld_results) {
  
  robustness_list <- list()
  
  # For each variable
  for (phenotype in names(ld_results$by_variable)) {
    # Extract consensus sets
    consensus_sets <- ld_results$by_variable[[phenotype]]$by_consensus_set
    
    # For each consensus set
    for (consensus_set in names(consensus_sets)) {
      # Extract LD result
      set_result <- consensus_sets[[consensus_set]]
      # Extract results before and after pruning for LD retention
      pre_ld <- set_result$pre_ld
      lead_snps <- set_result$lead_snps
      # Count numver of markers before and after pruning
      n_pre <- nrow(pre_ld)
      n_leads <- nrow(lead_snps)
      
      # Check for if there are no lead SNPs (no signal)
      if (n_leads == 0) {
        # Empty robustness data frame
        robustness_list[[paste(phenotype, consensus_set, sep = "__")]] <- data.frame(
          phenotype = phenotype,
          consensus_set = consensus_set,
          n_snps_pre_ld = n_pre,
          n_lead_snps_post_ld = 0,
          n_ld_blocks = 0,
          ld_retention_rate = ifelse(n_pre > 0, 0, NA_real_),
          mean_method_support = NA_real_,
          mean_env_method_support = NA_real_,
          evidence_score = NA_real_
        )
      } else {
        # Extract q values
        q_values <- lead_snps$min_q
        # Retain only valid q values
        q_values <- q_values[!is.na(q_values) & q_values > 0 & q_values <= 1]
        # Create robustness data frame
        robustness_list[[paste(phenotype, consensus_set, sep = "__")]] <- data.frame(
          phenotype = phenotype,
          consensus_set = consensus_set,
          n_snps_pre_ld = n_pre,
          n_lead_snps_post_ld = n_leads,
          n_ld_blocks = n_distinct(lead_snps$ld_block),
          ld_retention_rate = n_leads / n_pre,
          mean_method_support = mean(lead_snps$n_methods, na.rm = TRUE),
          mean_env_method_support = mean(lead_snps$n_env_methods, na.rm = TRUE),
          evidence_score = ifelse(length(q_values) > 0, median(-log10(q_values), na.rm = TRUE), NA_real_)
        )
      }
    }
  }
  
  # Return global robustness table
  bind_rows(robustness_list)
}

# scale_component()
# Scales metric between 0 and 1 within a phenotype while keeping empty consensus as NA
scale_component <- function(x, has_lead_snps) {
  x <- as.numeric(x)
  
  # Prevent empty sets from scoring
  x[!has_lead_snps] <- NA_real_

  # If all values are missing
  if (all(is.na(x))) {
    # Return numerical NA
    return(rep(NA_real_, length(x)))
  }
  
  # If all valid values are the same
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    # Score as 1
    out <- rep(NA_real_, length(x))
    out[!is.na(x)] <- 1
    return(out)
  }
  
  # Else, min max scale
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# score_consensus_robustness(robustness_summary, method_weight, ld_weight, evidence_weight)
# Scales metrics within each phenotype and calculates robustness score
score_consensus_robustness <- function(
    robustness_summary,
    method_weight = 0.4,
    ld_weight = 0.3,
    evidence_weight = 0.3
) {
  
  robustness_summary %>%
    # Group by variable
    group_by(phenotype) %>%
    mutate(
      has_lead_snps = n_lead_snps_post_ld > 0,
      # Scaled method support score
      method_support_score = scale_component(mean_method_support, has_lead_snps),
      # Scaled LD retention score
      ld_retention_score = scale_component(ld_retention_rate, has_lead_snps),
      # Scaled evidence score
      evidence_score_scaled = scale_component(evidence_score, has_lead_snps),
      # Scaled robustness score, NA treated as 0
      robustness_score =
        method_weight * coalesce(method_support_score, 0) +
        ld_weight * coalesce(ld_retention_score, 0) +
        evidence_weight * coalesce(evidence_score_scaled, 0),
    ) %>%
    ungroup()
}

# run_consensus_robustness(ld_results, output_dir)
# Evaluates LD pruned consensus sets and selects a representative primary set
run_consensus_robustness <- function(
    ld_results,
    output_dir = "Output/Consensus_Robustness"
) {
  
  message("\nRunning consensus robustness evaluation")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate robustness metrics
  robustness_summary <- summarize_robustness_all_sets(ld_results)
  # Calculate robustness score
  robustness_scored <- score_consensus_robustness(
    robustness_summary = robustness_summary,
    method_weight = 0.4,
    ld_weight = 0.3,
    evidence_weight = 0.3
  )
  
  # Create lead SNP signatures to detect identical LD-pruned sets
  signatures <- list()
  
  # For each variable
  for (phenotype in names(ld_results$by_variable)) {
    # Extract consensus sets
    consensus_sets <- ld_results$by_variable[[phenotype]]$by_consensus_set
    # For each consensus set
    for (consensus_set in names(consensus_sets)) {
      # Extract lead SNPs
      lead_snps <- consensus_sets[[consensus_set]]$lead_snps
      # if no signal
      if (is.null(lead_snps) || nrow(lead_snps) == 0) {
        lead_snp_signature <- NA_character_
      } else {
        # Extract lead SNPs (only unique)
        lead_snp_signature <- paste(sort(unique(as.character(lead_snps$marker))), collapse = ";")
      }
      # Create data frame per consensus set of variable, set and lead SNP signature (unique lead SNPs)
      signatures[[paste(phenotype, consensus_set, sep = "__")]] <- data.frame(
        phenotype = phenotype,
        consensus_set = consensus_set,
        lead_snp_signature = lead_snp_signature
      )
    }
  }
  
  # Bind into single data frame
  signatures <- bind_rows(signatures)
  
  # Evaluate empty sets, identical sets and primary consensus set
  robustness_evaluation <- robustness_scored %>%
    # Add lead SNP signature to metrics data frame
    left_join(signatures, by = c("phenotype", "consensus_set")) %>%
    mutate(
      # if no adaptive signal
      is_empty = is.na(n_lead_snps_post_ld) | n_lead_snps_post_ld == 0,
      # Prioritize sets by criteria for identical set cases (most strict set is preferred)
      consensus_priority = case_when(consensus_set == "high_confidence" ~ 3,
                                     consensus_set == "env_2methods" ~ 2,
                                     consensus_set == "broad_2methods" ~ 1,
                                     TRUE ~ 0)
    ) %>%
    # Group by phenotype (selection is done by variable)
    group_by(phenotype) %>%
    # Order
    arrange(
      is_empty, # non-empty
      desc(robustness_score), # robustness score
      desc(consensus_priority), # consensus priority type (for identical sets)
      desc(mean_env_method_support), # environmental method support
      desc(mean_method_support), # method support
      desc(ld_retention_rate), # ld retention rate
      desc(evidence_score), # evidence score
      .by_group = TRUE
    ) %>%
    mutate(
      # First row as primary consensus set
      primary_consensus_set = row_number() == 1 & !is_empty,
      primary_signature = lead_snp_signature[primary_consensus_set][1],
      # Check for non empty sets identical to selected primary set
      equivalent_to_primary = !is_empty & !is.na(lead_snp_signature) & !is.na(primary_signature) &
        lead_snp_signature == primary_signature,
      # Selection criteria
      selection_reason = case_when(
        is_empty ~ "Discarded; consensus set was empty",
        primary_consensus_set & sum(equivalent_to_primary, na.rm = TRUE) > 1 ~
          "Selected; representative of identical non-empty consensus sets with highest robustness score",
        primary_consensus_set ~
          "Selected; highest robustness score among non-empty consensus sets",
        equivalent_to_primary ~
          "Equivalent to selected primary consensus set",
        TRUE ~
          "Not selected"
      )
    ) %>%
    # All results
    ungroup() %>%
    # Columns for final evaluation table
    select(
      phenotype,
      consensus_set,
      robustness_score,
      n_ld_blocks,
      n_snps_pre_ld,
      n_lead_snps_post_ld,
      ld_retention_rate,
      mean_method_support,
      mean_env_method_support,
      evidence_score,
      is_empty,
      equivalent_to_primary,
      primary_consensus_set,
      selection_reason,
      lead_snp_signature
    )
  
  # Extract primary summary
  primary_summary <- robustness_evaluation %>%
    filter(primary_consensus_set)
  
  # Extract selected primary LD pruned lead SNPs
  primary_lead_snps <- list()
  
  for (i in seq_len(nrow(primary_summary))) {
    
    phenotype_i <- primary_summary$phenotype[i]
    consensus_i <- primary_summary$consensus_set[i]
    
    lead_snps_i <- ld_results$by_variable[[phenotype_i]]$by_consensus_set[[consensus_i]]$lead_snps
    
    if (!is.null(lead_snps_i) && nrow(lead_snps_i) > 0) {
      primary_lead_snps[[phenotype_i]] <- lead_snps_i %>%
        mutate(
          marker = as.character(marker),
          chr = as.character(chr),
          phenotype = as.character(phenotype),
          consensus_set = as.character(consensus_set),
          ld_block = as.character(ld_block),
          selected_primary = TRUE,
          primary_consensus_set = as.character(consensus_i)
        )
    }
  }
  
  primary_lead_snps <- if (length(primary_lead_snps) > 0) {
    bind_rows(primary_lead_snps)
  } else {
    data.frame()
  }
  
  # Write outputs
  write_csv(robustness_scored, file.path(output_dir, "consensus_robustness_metrics.csv"))
  write_csv(robustness_evaluation, file.path(output_dir, "consensus_robustness_evaluation.csv"))
  write_csv(primary_summary, file.path(output_dir, "primary_consensus_summary.csv"))
  write_csv(primary_lead_snps, file.path(output_dir, "primary_ld_pruned_lead_snps.csv"))
  
  list(
    metrics = robustness_scored,
    evaluation = robustness_evaluation,
    primary_summary = primary_summary,
    primary_lead_snps = primary_lead_snps
  )
}
