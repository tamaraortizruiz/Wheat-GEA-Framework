# ---- Post hoc LD Block Pruning ----

# safe_name()
# Ensures valid characters for folders or files
safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

# check_ld_consensus_input()
# Validates one consensus SNP table before LD block definition according to expected columns
# consensus_df = Consensus data frame from past module
# phenotype = Bioclimatic variable
# consensus_name = Consensus set name
# Returns: Reformated and filtered consensus data frame
check_ld_consensus_input <- function(consensus_df,
                                     phenotype,
                                     consensus_name) {
  

  if (is.null(consensus_df) || nrow(consensus_df) == 0) {
    return(NULL)
  }
  
  required_cols <- c(
    "marker",
    "chr",
    "position",
    "n_methods",
    "methods",
    "n_env_methods",
    "env_methods",
    "pcadapt_support",
    "min_p",
    "min_q",
    "min_env_p",
    "min_env_q"
  )
  missing_cols <- setdiff(required_cols, colnames(consensus_df))
  if (length(missing_cols) > 0) {
    stop("Consensus table for ", phenotype, " / ", consensus_name, " is missing columns: ",
         paste(missing_cols, collapse = ", ")
    )
  }
  
  # Format consensus data frame
  # Filter invalid rows and duplicates
  consensus_df %>%
    mutate(
      phenotype = phenotype,
      consensus_set = consensus_name,
      marker = as.character(marker),
      chr = as.character(chr),
      position = as.numeric(position),
      n_methods = as.numeric(n_methods),
      methods = as.character(methods),
      n_env_methods = as.numeric(n_env_methods),
      env_methods = as.character(env_methods),
      pcadapt_support = as.logical(pcadapt_support),
      min_p = as.numeric(min_p),
      min_q = as.numeric(min_q),
      min_env_p = as.numeric(min_env_p),
      min_env_q = as.numeric(min_env_q)
    ) %>%
    filter(!is.na(marker), marker != "", !is.na(min_p), min_p > 0, min_p <= 1, !is.na(min_q), min_q >= 0, min_q <= 1) %>%
    distinct(marker, .keep_all = TRUE)
}

# plink_pairwise_ld()
# Uses PLINK to calculate pairwise LD among candidate consensus SNPs
# bfile = PLINK dataset prefix (.bed/.bim/.fam)
# snp_file = Text file containing the candidate SNP IDs (one per line)
# out_prefix = Prefix for output PLINK .ld and .log files
# plink = Path to PLINK executable
# ld_window_kb = LD window size in kb
# ld_window = Number of SNPs per window
# ld_window_r2 = LD threshold
# allow_extra_chr = Allows non-standard chromosomes (A, B, D)
# Returns: Path to .ld file
plink_pairwise_ld <- function(
    bfile,
    snp_file,
    out_prefix,
    plink = "plink",
    ld_window_kb = 10000,
    ld_window = 99999,
    ld_window_r2 = 0,
    allow_extra_chr = TRUE
) {
  
  # PLINK arguments
  args <- c(
    "--bfile", bfile,
    "--extract", snp_file,
    "--r2",
    "--ld-window-kb", ld_window_kb,
    "--ld-window", ld_window,
    "--ld-window-r2", ld_window_r2,
    "--out", out_prefix,
    "--allow-no-sex"
  )
  
  if (allow_extra_chr) {
    args <- c(args, "--allow-extra-chr")
  }
  
  log_file <- paste0(out_prefix, "_plink.log")
  
  # Run PLINK
  status <- system2(
    command = plink,
    args = args,
    stdout = log_file,
    stderr = log_file
  )
  
  if (!identical(status, 0L)) {
    stop("PLINK LD calculation failed: ", log_file)
  }
  
  ld_file <- paste0(out_prefix, ".ld")
  if (!file.exists(ld_file)) {
    stop("PLINK did not create LD file: ", ld_file)
  }
  
  ld_file
}

# read_plink_ld()
# Reads PLINK .ld output and standardizes
# ld_file: Path to .ld file output
# Returns: LD pairs data frame
read_plink_ld <- function(ld_file) {
  
  if (!file.exists(ld_file)) {
    stop("LD file not found: ", ld_file)
  }
  
  if (file.info(ld_file)$size == 0) {
    return(data.frame(
      snp_a = character(),
      snp_b = character(),
      r2 = numeric()
    ))
  }
  
  # Read PLINK file (safe)
  ld <- tryCatch(
    read.table(file = ld_file, header = FALSE, skip = 1, stringsAsFactors = FALSE, fill = TRUE),
    error = function(e) {
      if (grepl("no lines available in input", conditionMessage(e))) {
        return(data.frame())
      } else {
        stop(e)
      }
    }
  )
  
  # No LD pairs
  if (nrow(ld) == 0) {
    return(data.frame(
      snp_a = character(),
      snp_b = character(),
      r2 = numeric()
    ))
  }
  
  # PLINK output formaat check 
  if (ncol(ld) < 7) {
    stop("PLINK LD file has fewer than 7 columns: ", ld_file)
  }
  
  ld <- ld[, 1:7]
  colnames(ld) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
  
  # Format LD pairs
  ld %>%
    transmute(
      snp_a = as.character(SNP_A),
      snp_b = as.character(SNP_B),
      r2 = as.numeric(R2)
    ) %>%
    filter(
      !is.na(snp_a),
      !is.na(snp_b),
      !is.na(r2),
      snp_a != snp_b
    ) %>%
    as.data.frame()
}

# define_ld_blocks()
# Defines LD blocks as connected components of SNPs with r2 >= threshold
# Graph theory where:
#   SNP = node
#   LD relationship = edge
#   r² >= threshold = connected SNPs
#   connected group = LD block
# consensus_df = Standardized consensus SNP data frame
# ld_pairs = LD pairs data frame
# r2_threshold = r2 value to use as threshold
define_ld_blocks <- function(
    consensus_df,
    ld_pairs,
    r2_threshold = 0.2
) {
  
  markers <- unique(consensus_df$marker)
  
  if (length(markers) == 0) {
    return(data.frame())
  }
  
  # LD pairs above r2 threshold
  ld_edges <- ld_pairs %>%
    filter(
      snp_a %in% markers,
      snp_b %in% markers,
      r2 >= r2_threshold
    )
  
  # no LD edges = each SNP is its own single block
  if (nrow(ld_edges) == 0) {
    return(consensus_df %>%
             arrange(chr, position) %>%
             transmute(
               marker = marker,
               ld_block = paste0("LD", row_number()),
               block_size = 1,
               max_r2_in_block = NA_real_ # numeric NA
             ) %>%
             as.data.frame())
  }
  
  # Graph vertices
  vertices <- data.frame(name = markers)
  
  # Undirected graph
  graph <- graph_from_data_frame(
    d = ld_edges %>% dplyr::select(snp_a, snp_b),
    directed = FALSE,
    vertices = vertices)
  
  # Find connected components
  components <- igraph::components(graph)$membership
  
  block_df <- data.frame(
    marker = names(components),
    component_id = as.integer(components)
  )
  block_df <- block_df %>%
    left_join(consensus_df %>% dplyr::select(marker, chr, position), by = "marker")
  
  block_order <- block_df %>%
    group_by(component_id) %>%
    summarise(
      chr = dplyr::first(na.omit(chr)),
      min_position = suppressWarnings(min(position, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      min_position = ifelse(is.infinite(min_position), NA_real_, min_position)
    ) %>%
    arrange(chr, min_position) %>%
    mutate(ld_block = paste0("LD", row_number()))
  
  block_df <- block_df %>%
    left_join(block_order %>% dplyr::select(component_id, ld_block), by = "component_id")
  
  # Maximum R2 per block
  edge_blocks <- ld_edges %>%
    # Attach component ID to SNP A
    left_join(block_df %>% dplyr::select(marker, component_id), by = c("snp_a" = "marker")) %>%
    dplyr::rename(component_a = component_id) %>%
    # Attach component ID to SNP B
    left_join(block_df %>% dplyr::select(marker, component_id), by = c("snp_b" = "marker")) %>%
    dplyr::rename(component_b = component_id) %>%
    filter(component_a == component_b) %>%
    group_by(component_a) %>%
    # Max R2 in block
    summarise(
      max_r2_in_block = max(r2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(component_id = component_a)
  
  block_df %>%
    group_by(ld_block) %>%
    mutate(block_size = n()) %>%
    ungroup() %>%
    left_join(edge_blocks, by = "component_id") %>%
    dplyr::select(marker, ld_block, block_size, max_r2_in_block) %>%
    arrange(ld_block) %>%
    as.data.frame()
}

# choose_lead_snps()
# Chooses one lead SNP per LD block
# consensus_df = Standardized consensus SNP data frame
# ld_blocks = LD blocks data frame
# Returns: list with data frame of all SNPs with LD block and lead SNP status and lead SNP data frame
choose_lead_snps <- function(consensus_df, ld_blocks) {
  
  annotated <- consensus_df %>%
    dplyr::left_join(ld_blocks, by = "marker")
  
  # Add missing optional columns if needed
  if (!"n_env_methods" %in% colnames(annotated)) {
    annotated$n_env_methods <- NA_real_
  }
  
  if (!"pcadapt_support" %in% colnames(annotated)) {
    annotated$pcadapt_support <- FALSE
  }
  
  if (!"min_q" %in% colnames(annotated)) {
    annotated$min_q <- NA_real_
  }
  
  if (!"min_env_q" %in% colnames(annotated)) {
    annotated$min_env_q <- NA_real_
  }
  
  if (!"min_env_p" %in% colnames(annotated)) {
    annotated$min_env_p <- NA_real_
  }
  
  # Lead SNP selection
  lead_snps <- annotated %>%
    group_by(.data$ld_block) %>%
    arrange(
      desc(coalesce(.data$n_methods, -Inf)),
      desc(coalesce(.data$n_env_methods, -Inf)),
      desc(coalesce(as.integer(.data$pcadapt_support), 0L)),
      coalesce(.data$min_env_q, Inf),
      coalesce(.data$min_q, Inf),
      coalesce(.data$min_env_p, Inf),
      coalesce(.data$min_p, Inf),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    mutate(is_lead = TRUE)
  
  all_candidates <- annotated %>%
    mutate(is_lead = .data$marker %in% lead_snps$marker) %>%
    arrange(.data$chr, .data$position, .data$ld_block, desc(.data$is_lead))
  
  list(
    all_candidates = all_candidates,
    lead_snps = lead_snps
  )
}

# summarize_ld_blocks()
# Summarizes LD block pruning results
# consensus_df = Standardized consensus SNP data frame
# lead_snps = Lead SNPs data frame
# phenotype = Bioclimatic variable
# consensus_name = Consensus set name
summarize_ld_blocks <- function(
    consensus_df,
    lead_snps,
    phenotype,
    consensus_name
) {
  
  data.frame(
    phenotype = phenotype,
    consensus_set = consensus_name,
    n_snps_pre_ld = nrow(consensus_df),
    n_lead_snps_post_ld = nrow(lead_snps),
    n_ld_blocks = n_distinct(lead_snps$ld_block),
    ld_retention_rate = nrow(lead_snps) / nrow(consensus_df),
    mean_method_support = mean(lead_snps$n_methods, na.rm = TRUE),
    evidence_score = median(-log10(lead_snps$min_q), na.rm = TRUE)
  )
}

# run_ld_blocks_one_set()
# Runs PLINK + R LD block definition for one consensus SNP set
# consensus_df = Standardized consensus SNP data frame
# bfile = PLINK dataset prefix (.bed/.bim/.fam)
# phenotype = Bioclimatic variable
# consensus_name = Consensus set name
# output_dir = Output directory
# r2_threshold = R2 threshold value
# ld_window_kb = LD window in kb
# ld_window = Number of SNPs in window
# allow_extra_chr = TRUE by default to allow non-standard chromosomes
# Returns: list of results
run_ld_blocks_one_set <- function(
    consensus_df,
    bfile,
    phenotype,
    consensus_name,
    output_dir = "Output/Consensus_LD",
    plink = "plink",
    r2_threshold = 0.2,
    ld_window_kb = 10000,
    ld_window = 99999,
    allow_extra_chr = TRUE
) {
  
  message("LD blocks: ", phenotype, " | ", consensus_name)
  
  set_dir <- file.path(output_dir, safe_name(phenotype), safe_name(consensus_name))
  dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)
  
  consensus_df <- check_ld_consensus_input(
    consensus_df = consensus_df,
    phenotype = phenotype,
    consensus_name = consensus_name
  )
  
  # Empty consensus set check
  if (is.null(consensus_df) || nrow(consensus_df) == 0) {
    empty_summary <- data.frame(
      phenotype = phenotype,
      consensus_set = consensus_name,
      n_snps_pre_ld = 0,
      n_lead_snps_post_ld = 0,
      n_ld_blocks = 0,
      ld_retention_rate = NA_real_,
      mean_method_support = NA_real_,
      evidence_score = NA_real_
    )
    
    write_csv(empty_summary, file.path(set_dir, "ld_block_summary.csv"))
    
    return(list(
      pre_ld = data.frame(),
      ld_pairs = data.frame(),
      ld_blocks = data.frame(),
      all_candidates = data.frame(),
      lead_snps = data.frame(),
      summary = empty_summary
    ))
  }
  
  write_csv(consensus_df, file.path(set_dir, "consensus_pre_ld.csv"))
  
  # Single SNP case
  if (nrow(consensus_df) == 1) {
    # Single block for SNP
    ld_blocks <- consensus_df %>%
      transmute(
        marker = marker,
        ld_block = "LD1",
        block_size = 1,
        max_r2_in_block = NA_real_
      )
    
    lead_results <- choose_lead_snps(consensus_df, ld_blocks)
    
    summary <- summarize_ld_blocks(
      consensus_df = consensus_df,
      lead_snps = lead_results$lead_snps,
      phenotype = phenotype,
      consensus_name = consensus_name
    )
    
    write_csv(ld_blocks, file.path(set_dir, "ld_blocks.csv"))
    write_csv(lead_results$all_candidates, file.path(set_dir, "ld_annotated_candidates.csv"))
    write_csv(lead_results$lead_snps, file.path(set_dir, "ld_pruned_lead_snps.csv"))
    write_csv(summary, file.path(set_dir, "ld_block_summary.csv"))
    
    return(list(
      pre_ld = consensus_df,
      ld_pairs = data.frame(),
      ld_blocks = ld_blocks,
      all_candidates = lead_results$all_candidates,
      lead_snps = lead_results$lead_snps,
      summary = summary
    ))
  }
  
  # For all other >1 SNP cases
  snp_file <- file.path(set_dir, "candidate_snps.txt")
  write_lines(consensus_df$marker, snp_file)
  
  plink_out <- file.path(set_dir, "candidate_pairwise_ld")
  ld_file <- plink_pairwise_ld(
    bfile = bfile,
    snp_file = snp_file,
    out_prefix = plink_out,
    plink = plink,
    ld_window_kb = ld_window_kb,
    ld_window = ld_window,
    ld_window_r2 = r2_threshold,
    allow_extra_chr = allow_extra_chr
  )
  
  ld_pairs <- read_plink_ld(ld_file)
  
  # Define LD blocks
  ld_blocks <- define_ld_blocks(
    consensus_df = consensus_df,
    ld_pairs = ld_pairs,
    r2_threshold = r2_threshold
  )
  
  # Choose lead SNPs
  lead_results <- choose_lead_snps(
    consensus_df = consensus_df,
    ld_blocks = ld_blocks
  )
  
  # LD block summary
  summary <- summarize_ld_blocks(
    consensus_df = consensus_df,
    lead_snps = lead_results$lead_snps,
    phenotype = phenotype,
    consensus_name = consensus_name
  )
  
  write_csv(ld_pairs, file.path(set_dir, "plink_pairwise_ld.csv"))
  write_csv(ld_blocks, file.path(set_dir, "ld_blocks.csv"))
  write_csv(lead_results$all_candidates, file.path(set_dir, "ld_annotated_candidates.csv"))
  write_csv(lead_results$lead_snps, file.path(set_dir, "ld_pruned_lead_snps.csv"))
  write_csv(summary, file.path(set_dir, "ld_block_summary.csv"))
  
  list(
    pre_ld = consensus_df,
    ld_pairs = ld_pairs,
    ld_blocks = ld_blocks,
    all_candidates = lead_results$all_candidates,
    lead_snps = lead_results$lead_snps,
    summary = summary
  )
}

# run_ld_blocks_single_variable()
# Runs LD block pruning for all four consensus sets of one environmental variable
# consensus_variable_result = Consensus SNP set creation output object for a single variable (all 4 consensus sets)
# bfile = PLINK dataset prefix (.bed/.bim/.fam)
# phenotype = Bioclimatic variable
# output_dir = Output directory
# r2_threshold = R2 threshold value
# ld_window_kb = LD window in kb
# ld_window = Number of SNPs in window
# allow_extra_chr = TRUE by default to allow non-standard chromosomes
run_ld_blocks_single_variable <- function(
    consensus_variable_result,
    bfile,
    phenotype,
    output_dir = "Output/Consensus_LD",
    plink = "plink",
    r2_threshold = 0.2,
    ld_window_kb = 10000,
    ld_window = 99999,
    allow_extra_chr = TRUE
) {
  
  message("\nRunning LD block pruning for: ", phenotype)
  
  consensus_sets <- list(
    broad_2methods = consensus_variable_result$broad_2methods,
    env_2methods = consensus_variable_result$env_2methods,
    high_confidence = consensus_variable_result$high_confidence
  )
  
  ld_results <- list()
  
  for (consensus_name in names(consensus_sets)) {
    ld_results[[consensus_name]] <- run_ld_blocks_one_set(
      consensus_df = consensus_sets[[consensus_name]],
      bfile = bfile,
      phenotype = phenotype,
      consensus_name = consensus_name,
      output_dir = output_dir,
      plink = plink,
      r2_threshold = r2_threshold,
      ld_window_kb = ld_window_kb,
      ld_window = ld_window,
      allow_extra_chr = allow_extra_chr
    )
  }
  
  combined_summary <- bind_rows(lapply(ld_results, function(x) x$summary))
  
  phenotype_dir <- file.path(output_dir, safe_name(phenotype))
  dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(combined_summary,file.path(phenotype_dir, "ld_block_summary_all_sets.csv"))
  
  list(
    by_consensus_set = ld_results,
    summary = combined_summary
  )
}

# run_ld_blocks_all_variables()
# Runs post hoc LD block pruning for all environmental variables
# config = YAML configuration list
# consensus_results = Consensus SNP set results object for all climatic variables
# phenotypes = List of phenotypes for analysis
run_ld_blocks_all_variables <- function(
    config,
    qc_prefix = config$qc_outputs$qc_prefix,
    consensus_results,
    phenotypes = config$env$vars
) {
  
  required_files <- paste0(qc_prefix, c(".bed", ".bim", ".fam"))
  if (!all(file.exists(required_files))) {
    stop("Post hoc LD pruning PLINK files were not found")
  }
  
  ld_all <- list()
  
  # For each phenotype
  for (phenotype in phenotypes) {
    if (!phenotype %in% names(consensus_results$by_variable)) {
      warning("Skipping ", phenotype, ": not found in consensus_results$by_variable.")
      next
    }
    ld_all[[phenotype]] <- run_ld_blocks_single_variable(
      consensus_variable_result = consensus_results$by_variable[[phenotype]],
      bfile = qc_prefix,
      phenotype = phenotype,
      output_dir = config$ld_pruning_ph$output_dir,
      plink = config$plink$path,
      r2_threshold = config$ld_pruning_ph$r2_threshold,
      ld_window_kb = config$ld_pruning_ph$ld_window_kb,
      ld_window = config$ld_pruning_ph$ld_window,
      allow_extra_chr = config$ld_pruning_ph$allow_extra_chr
    )
  }
  
  combined_summary <- bind_rows(lapply(ld_all, function(x) x$summary))
  
  dir.create(config$ld_pruning_ph$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(combined_summary, file.path(config$ld_pruning_ph$output_dir, "ld_block_summary_all_variables.csv"))
  
  list(
    by_variable = ld_all,
    summary = combined_summary
  )
}

# ld_block_files_exist()
# Checks whether LD block results already exist
# config = YAML configuration list
# phenotypes = List of bioclimatic variables
ld_block_files_exist <- function(
    config,
    phenotypes = config$env$vars
) {
  
  output_dir <- config$ld_pruning_ph$output_dir
  
  required_files <- c(
    file.path(output_dir, "ld_block_summary_all_variables.csv"),
    unlist(lapply(phenotypes, function(phenotype) {
      file.path(output_dir, safe_name(phenotype), "ld_block_summary_all_sets.csv")
    }))
  )
  
  all(file.exists(required_files))
}

# load_ld_block_set()
# Loads LD block results for one phenotype and one consensus set
# output_dir = LD block output directory
# phenotype = Bioclimatic variable
# consensus_name = Consensus set name
load_ld_block_set <- function(
    output_dir,
    phenotype,
    consensus_name
) {
  
  set_dir <- file.path(output_dir, safe_name(phenotype), safe_name(consensus_name))
  
  list(
    pre_ld = read_csv_safe(file.path(set_dir, "consensus_pre_ld.csv")),
    ld_pairs = read_csv_safe(file.path(set_dir, "plink_pairwise_ld.csv")),
    ld_blocks = read_csv_safe(file.path(set_dir, "ld_blocks.csv")),
    all_candidates = read_csv_safe(file.path(set_dir, "ld_annotated_candidates.csv")),
    lead_snps = read_csv_safe(file.path(set_dir, "ld_pruned_lead_snps.csv")),
    summary = read_csv_safe(file.path(set_dir, "ld_block_summary.csv"))
  )
}

# load_ld_block_results()
# Loads saved post-hoc LD block results
# config = YAML configuration list
# phenotypes = List of bioclimatic variables
load_ld_block_results <- function(
    config,
    phenotypes = config$env$vars
) {
  
  output_dir <- config$ld_pruning_ph$output_dir
  consensus_sets <- c("broad_2methods", "env_2methods", "high_confidence")
  
  ld_all <- list()
  
  # For each phenotype
  for (phenotype in phenotypes) {
    
    set_results <- list()
    
    for (consensus_name in consensus_sets) {
      set_results[[consensus_name]] <- load_ld_block_set(
        output_dir = output_dir,
        phenotype = phenotype,
        consensus_name = consensus_name
      )
    }
    
    phenotype_summary_file <- file.path(output_dir, safe_name(phenotype), "ld_block_summary_all_sets.csv")
    
    ld_all[[phenotype]] <- list(
      by_consensus_set = set_results,
      summary = read_csv_safe(phenotype_summary_file)
    )
  }
  
  combined_summary <- read_csv_safe(file.path(output_dir, "ld_block_summary_all_variables.csv"))
  
  list(
    by_variable = ld_all,
    summary = combined_summary
  )
}