# ---- Primary LD-Pruned Consensus SNPs ----

# create_primary_snp_results()
# Extracts LD-pruned lead SNPs from the configured primary consensus set
# config = Pipeline configuration object
# ld_results = Output from the LD processing module
# phenotypes = Environmental variables to include
# output_dir = Directory for combined primary-set tables
# Returns: List containing set, summary and lead_snps
create_primary_snp_results <- function(
    config,
    ld_results,
    phenotypes = config$env$vars,
    output_dir = config$ld_pruning_ph$output_dir
) {
  
  primary_set <- config$consensus$primary_set
  
  if (is.null(primary_set) || length(primary_set) == 0 ||
      (length(primary_set) == 1 && (is.na(primary_set) || primary_set == ""))) {
    primary_set <- "env_2methods"
  }
  
  primary_set <- validate_primary_set(primary_set)
  summary_list <- list()
  snp_list <- list()
  
  for (phenotype in phenotypes) {
    phenotype_result <- ld_results$by_variable[[phenotype]]
    
    if (is.null(phenotype_result)) {
      warning("Skipping ", phenotype, ": no LD results were found.")
      next
    }
    
    set_result <- phenotype_result$by_consensus_set[[primary_set]]
    
    if (is.null(set_result)) {
      stop(
        "LD results for ", phenotype,
        " do not contain the configured primary set: ", primary_set
      )
    }
    
    summary_list[[phenotype]] <- set_result$summary %>%
      mutate(
        phenotype = .env$phenotype,
        consensus_set = .env$primary_set
      )
    
    if (!is.null(set_result$lead_snps) && nrow(set_result$lead_snps) > 0) {
      snp_list[[phenotype]] <- set_result$lead_snps %>%
        mutate(
          phenotype = as.character(.env$phenotype),
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
          min_env_q = as.numeric(min_env_q),
          consensus_set = as.character(.env$primary_set),
          ld_block = as.character(ld_block),
          ld_block_size = as.integer(block_size),
          max_r2_in_block = as.numeric(max_r2_in_block)
        )
    }
  }
  
  primary_summary <- bind_rows(summary_list)
  primary_lead_snps <- bind_rows(snp_list)
  
  if (nrow(primary_lead_snps) == 0) {
    primary_lead_snps <- data.frame(
      phenotype = character(),
      marker = character(),
      chr = character(),
      position = numeric(),
      n_methods = numeric(),
      methods = character(),
      n_env_methods = numeric(),
      env_methods = character(),
      pcadapt_support = logical(),
      min_p = numeric(),
      min_q = numeric(),
      min_env_p = numeric(),
      min_env_q = numeric(),
      consensus_set = character(),
      ld_block = character(),
      ld_block_size = integer(),
      max_r2_in_block = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    primary_lead_snps <- primary_lead_snps %>%
      mutate(
        phenotype = as.character(phenotype),
        marker = as.character(marker),
        chr = as.character(chr),
        position = as.numeric(position),
        consensus_set = as.character(consensus_set),
        ld_block = as.character(ld_block)
      ) %>%
      distinct(phenotype, marker, .keep_all = TRUE) %>%
      arrange(phenotype, chr, position, marker)
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(primary_summary, file.path(output_dir, "primary_ld_summary.csv"))
  write_csv(primary_lead_snps, file.path(output_dir, "final_primary_ld_pruned_snps.csv"))
  
  list(
    set = primary_set,
    summary = primary_summary,
    lead_snps = primary_lead_snps
  )
}

plot_primary_qtl <- function(
    primary_lead_snps,
    bim_file,
    phenotypes = NULL,
    direction_table = NULL,
    output_file = NULL
) {
  
  if (nrow(primary_lead_snps) == 0) {
    stop("No primary LD-pruned SNPs available for plotting.")
  }
  
  # Wheat chromosome order
  wheat_chr_order <- c(
    "1A", "1B", "1D",
    "2A", "2B", "2D",
    "3A", "3B", "3D",
    "4A", "4B", "4D",
    "5A", "5B", "5D",
    "6A", "6B", "6D",
    "7A", "7B", "7D"
  )
  
  # Read BIM for chromosome sizes
  bim <- read_table(
    bim_file,
    col_names = c("chr", "marker", "cm", "position", "a1", "a2"),
    show_col_types = FALSE
  ) %>%
    mutate(
      chr = as.character(chr),
      chr = gsub("^chr", "", chr, ignore.case = TRUE),
      marker = as.character(marker),
      position = as.numeric(position)
    ) %>%
    filter(
      !is.na(chr),
      !is.na(position),
      chr != "",
      chr != "NA",
      chr != ".",
      chr != "0"
    )
  
  # Extract chromosome info
  chr_info <- bim %>%
    group_by(chr) %>%
    summarise(chr_len = max(position, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      chr = factor(chr, levels = wheat_chr_order)
    ) %>%
    filter(!is.na(chr)) %>%
    arrange(chr) %>%
    mutate(
      chr = as.character(chr),
      chr_start = lag(cumsum(chr_len), default = 0),
      chr_end = chr_start + chr_len,
      chr_mid = chr_start + chr_len / 2,
      chr_index = row_number()
    )
  
  # Prepare primary lead SNPs
  plot_df <- primary_lead_snps %>%
    mutate(
      phenotype = as.character(phenotype),
      marker = as.character(marker),
      chr = as.character(chr),
      chr = gsub("^chr", "", chr, ignore.case = TRUE),
      position = as.numeric(position)
    ) %>%
    filter(
      !is.na(phenotype),
      !is.na(marker),
      !is.na(chr),
      !is.na(position),
      chr != "",
      chr != "NA",
      chr != ".",
      chr != "0"
    )
  
  if (!is.null(phenotypes)) {
    plot_df <- plot_df %>%
      filter(phenotype %in% phenotypes)
  }
  
  if (nrow(plot_df) == 0) {
    stop("No primary lead SNPs available for plotting.")
  }
  
  # Add allelic direction if available
  if (!is.null(direction_table)) {
    direction_table <- direction_table %>%
      mutate(
        phenotype = as.character(phenotype),
        marker = as.character(marker),
        adaptive_direction = as.numeric(adaptive_direction)
      ) %>%
      distinct(phenotype, marker, .keep_all = TRUE) %>%
      dplyr::select(phenotype, marker, adaptive_direction)

    plot_df <- plot_df %>%
      left_join(direction_table, by = c("phenotype", "marker")) %>%
      mutate(
        allelic_effect = case_when(
          adaptive_direction > 0 ~ "positive",
          adaptive_direction < 0 ~ "negative",
          TRUE ~ "unknown"
        )
      )
  } else {
    plot_df <- plot_df %>%
      mutate(allelic_effect = "primary lead SNP")
  }
  
  # Add cumulative genomic position
  plot_df <- plot_df %>%
    inner_join(chr_info %>% dplyr::select(chr, chr_start), by = "chr") %>%
    mutate(
      genome_pos = chr_start + position,
      phenotype = factor(phenotype, levels = rev(unique(phenotype)))
    ) %>%
    distinct(phenotype, marker, chr, position, .keep_all = TRUE)
  
  # Chromosome background
  chr_bg <- chr_info %>%
    mutate(
      fill_group = chr_index %% 2
    )
  
  p <- ggplot() +
    geom_rect(
      data = chr_bg,
      aes(
        xmin = chr_start,
        xmax = chr_end,
        ymin = -Inf,
        ymax = Inf,
        fill = factor(fill_group)
      ),
      alpha = 0.35
    ) +
    geom_point(
      data = plot_df,
      aes(
        x = genome_pos,
        y = phenotype,
        color = allelic_effect,
        size = n_methods
      ),
      alpha = 0.85
    ) +
    scale_x_continuous(
      breaks = chr_info$chr_mid,
      labels = chr_info$chr,
      expand = c(0.01, 0.01)
    ) +
    scale_fill_manual(
      values = c("grey85", "grey70"),
      guide = "none"
    ) +
    labs(
      title = "Primary lead SNPs across chromosomes",
      x = "Chromosomes",
      y = "Environmental variables",
      color = "Allelic effect",
      size = "Method support"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.position = "right"
    )
  
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = 13,
      height = 7,
      dpi = 300
    )
  }
  
  return(p)
}
