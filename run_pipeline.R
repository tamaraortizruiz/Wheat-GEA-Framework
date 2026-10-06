#!/usr/bin/env Rscript

# Pipeline can be ran as selected computational stages
# E.x.
# full pipeline
#   Rscript run_pipeline.R
# as configured in config.yaml
#   Rscript run_pipeline.R --config config.yaml
# selection through command-line
#   Rscript run_pipeline.R --stages env,qc,kinship,ld_pruning,pca
#   Rscript run_pipeline.R --stages consensus,consensus_ld,primary_snps

# Set wd to the repository root
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  setwd(dirname(script_path))
}

# Recognizes: 
#   --config
#   --stages
#   --help
parse_pipeline_args <- function(args) {
  out <- list(config = "config.yaml", stages = NULL)
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (grepl("^--config=", arg)) {
      out$config <- sub("^--config=", "", arg)
    } else if (identical(arg, "--config") && i < length(args)) {
      i <- i + 1
      out$config <- args[[i]]
    } else if (grepl("^--stages=", arg)) {
      out$stages <- strsplit(sub("^--stages=", "", arg), ",", fixed = TRUE)[[1]]
    } else if (identical(arg, "--stages") && i < length(args)) {
      i <- i + 1
      out$stages <- strsplit(args[[i]], ",", fixed = TRUE)[[1]]
    } else if (arg %in% c("--help", "-h")) {
      cat(
        "Use:\n",
        "  Rscript run_pipeline.R\n",
        "  Rscript run_pipeline.R --config file\n",
        "  Rscript run_pipeline.R --stages a,b,c\n",
        "  Rscript run_pipeline.R --config FILE --stages a,b,c\n\n",
        "Stages: env, qc, kinship, ld_pruning, pca, lmm, lfmm, rda, pcadapt,\n",
        "        snp_overlap, consensus, consensus_ld, primary_snps,\n",
        "        adaptive_scoring, bio_interp\n"
      )
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg)
    }
    i <- i + 1
  }
  if (!is.null(out$stages)) {
    out$stages <- trimws(out$stages)
  }
  out
}

# Define arguments
args <- parse_pipeline_args(commandArgs(trailingOnly = TRUE))

# Load function .R files
function_files <- sort(list.files("Functions", pattern = "\\.R$", full.names = TRUE))
invisible(lapply(function_files, source))

# Load configuration YAML file
config <- yaml::read_yaml(args$config)
overwrite <- isTRUE(config$run$overwrite)

if (overwrite && !is.null(args$stages)) {
  stop(
    "overwrite is set to TRUE",
    "run the full pipeline without --stages or set overwrite to FALSE for a targeted run"
  )
}

# Align stage names with configuration switches
stage_flags <- c(
  env = "run_env",
  qc = "run_qc",
  kinship = "run_kinship",
  ld_pruning = "run_ld_pruning",
  pca = "run_pca",
  lmm = "run_lmm",
  lfmm = "run_lfmm",
  rda = "run_rda",
  pcadapt = "run_pcadapt",
  snp_overlap = "run_snp_overlap",
  consensus = "run_consensus",
  consensus_ld = "run_consensus_ld",
  primary_snps = "run_primary_snps",
  adaptive_scoring = "run_adaptive_scoring",
  bio_interp = "run_bio_interp"
)

if (!is.null(args$stages)) {
  unknown <- setdiff(args$stages, names(stage_flags))
  if (length(unknown) > 0) {
    stop("Unknown stage(s): ", paste(unknown, collapse = ", "))
  }
}

enabled <- function(stage) {
  if (!is.null(args$stages)) {
    return(stage %in% args$stages)
  }
  isTRUE(config$analysis[[stage_flags[[stage]]]])
}

stage_cache_exists <- function(stage) {
  file.exists(pipeline_stage_file(config, stage))
}

# Define function to run each stage
run_stage <- function(stage, code) {
  # Check whether stage is enabled
  if (!enabled(stage)) {
    message("[SKIP] ", stage)
    return(invisible(NULL))
  }

  started <- Sys.time()
  message("\n", strrep("=", 72), "\n[STAGE] ", stage, "\n", strrep("=", 72))
  status <- "completed"
  
  # Check for completion
  tryCatch(
    force(code),
    error = function(e) {
      status <<- "failed"
      stop(e)
    },
    finally = {
      finished <- Sys.time()
      # Update pipeline_manifest.csv
      update_pipeline_manifest(config, stage, status, started, finished)
    }
  )
  
  invisible(NULL)
}

# Project directory setup
create_project_dirs(overwrite = overwrite)
dir.create(pipeline_stage_dir(config), recursive = TRUE, showWarnings = FALSE)

# Environmental module stage
run_stage("env", {
  if (!overwrite && stage_cache_exists("env")) {
    message("Reusing existing environmental stage")
  } else {
    # Read metadata
    metadata_initial <- if (isTRUE(config$subsetting$use)) {
      read.csv(config$subsetting$sample_metadata, stringsAsFactors = FALSE)
    } else {
      read.csv(config$metadata$file, stringsAsFactors = FALSE)
    }

    # Extract WorldClim data
    climate_data <- extract_climate_variables(
      metadata = metadata_initial,
      lon_col = config$env$lon_col,
      lat_col = config$env$lat_col,
      sample_col = config$metadata$sample_col,
      var = config$env$var,
      res = config$env$res,
      climate_dir = config$env$path,
      output_file = config$env$output_file,
      overwrite = overwrite
    )
    
    # Run climate PCA
    climate_pca <- run_climate_pca(
      climate_data = climate_data,
      output_file = "Output/Structure/climate_pca.rds"
    )
    
    # Save stage
    save_pipeline_stage(
      list(
        metadata = metadata_initial,
        climate_data = climate_data,
        climate_pca = climate_pca
        ),
      config, "env"
    )
  }
})

# Quality control module stage
run_stage("qc", {
  if (!overwrite && stage_cache_exists("qc")) {
    message("Reusing existing QC stage")
  } else {
    # Load env stage results
    env_result <- load_environment_stage(config)
    climate_data <- env_result$climate_data
    qc_files <- paste0(config$qc_outputs$qc_prefix, c(".bed", ".bim", ".fam"))
    if (!overwrite && all(file.exists(qc_files))) {
      message("Reusing filtered PLINK files: ", config$qc_outputs$qc_prefix)
      qc_prefix <- config$qc_outputs$qc_prefix
    } else {
      # PLINK keep file for accessions with env data 
      keep_file <- create_keep_file(
        metadata = climate_data,
        fam_file = paste0(config$genotype$input_prefix, ".fam"),
        output_file = config$subsetting$keep_file,
        sample_col = config$metadata$sample_col
      )
      # PLINK filtering
      qc_prefix <- plink_filter(
        input = config$genotype$input_prefix,
        output = config$qc_outputs$qc_prefix,
        plink = config$plink$path,
        keep = keep_file,
        geno_na = config$qc$geno_na,
        maf = config$qc$maf,
        ind_na = config$qc$ind_na
      )
    }
    # if marker subset
    if (isTRUE(config$marker_subsetting$use_extract)) {
      qc_prefix <- plink_extract_markers(
        input = qc_prefix,
        output = config$marker_subsetting$output_prefix,
        plink = config$plink$path,
        extract_file = config$marker_subsetting$extract_file,
        overwrite = overwrite
      )
    }
    # Filter + order metadata
    metadata <- filter_metadata(
      metadata_file = config$metadata$file,
      fam_file = paste0(qc_prefix, ".fam"),
      sample_col = config$metadata$sample_col,
      output_file = config$metadata$filtered_file
    )
    # QC plots
    qc_plots <- plot_genotype_qc_overview(
      prefix = qc_prefix,
      output_dir = "Output/QC",
      dataset_name = "QC_genotype"
    )
    # Save stage
    save_pipeline_stage(
      list(
        qc_prefix = qc_prefix,
        metadata = metadata,
        climate_data = climate_data,
        dimensions = count_plink_dataset(qc_prefix),
        marker_summary = qc_plots$marker_summary
      ),
      config, "qc"
    )
  }
})

# Kinship estimation module stage
run_stage("kinship", {
  if (!overwrite && stage_cache_exists("kinship")) {
    message("Reusing existing kinship stage")
  } else {
    # Load QC stage
    qc_result <- load_qc_stage(config)
    qc_prefix <- qc_result$qc_prefix
    # Kinship estimation
    gemma_kinship_file <- run_gemma_kinship(
      gemma = config$gemma$path,
      bfile = qc_prefix,
      output_prefix = config$gemma$kinship_prefix,
      output_dir = config$gemma$output_dir,
      overwrite = overwrite
    )
    # Kinship summary
    kinship_stats <- summarize_gemma_kinship(gemma_kinship_file)
    kinship_diag <- diagnose_gemma_kinship(
      kinship_file = gemma_kinship_file,
      bfile = qc_prefix,
      duplicate_threshold = config$kinship_qc$duplicate_threshold,
      diag_mad_threshold = config$kinship_qc$diagonal_mad_threshold
    )
    # Kinship filtering
    kinship_filter <- filter_kinship_samples(
      kin_diag = kinship_diag,
      qc_prefix = qc_prefix,
      plink = config$plink$path,
      remove_duplicates = config$kinship_qc$remove_duplicates,
      remove_diagonal_outliers = config$kinship_qc$remove_diagonal_outliers
    )
    
    # Re-estimate kinship for filtered samples
    if (isTRUE(kinship_filter$filtered)) {
      gemma_kinship_file <- run_gemma_kinship(
        gemma = config$gemma$path,
        bfile = qc_prefix,
        output_prefix = config$gemma$kinship_prefix,
        output_dir = config$gemma$output_dir,
        overwrite = TRUE
      )
      kinship_stats <- summarize_gemma_kinship(gemma_kinship_file)
    }
    
    final_fam <- read.table(paste0(qc_prefix, ".fam"), stringsAsFactors = FALSE)
    final_ids <- final_fam[, 2]
    sample_col <- config$metadata$sample_col
    climate_match <- match(final_ids, qc_result$climate_data[[sample_col]])
    if (anyNA(climate_match)) {
      stop("Mismatch between final genotype samples and climate data")
    }
    
    # Refilter env data
    climate_data <- qc_result$climate_data[climate_match, , drop = FALSE]
    
    # Refilter metadata
    metadata <- filter_metadata(
      metadata_file = config$metadata$file,
      fam_file = paste0(qc_prefix, ".fam"),
      sample_col = sample_col,
      output_file = config$metadata$filtered_file
    )
    
    # Define downstream genotype matrix, metadata, climate data
    qc_result$metadata <- metadata
    qc_result$climate_data <- climate_data
    qc_result$dimensions <- count_plink_dataset(qc_prefix)
    
    save_pipeline_stage(qc_result, config, "qc")
    save_pipeline_stage(
      list(
        gemma_kinship_file = gemma_kinship_file,
        stats = kinship_stats,
        diagnostics = kinship_diag,
        filter = kinship_filter
      ),
      config, "kinship"
    )
  }
})

# LD pruning module stage
run_stage("ld_pruning", {
  if (!overwrite && stage_cache_exists("ld_pruning")) {
    message("Reusing cached LD-pruning stage.")
  } else {
    qc_result <- load_qc_stage(config)
    ld_files <- paste0(config$qc_outputs$ld_prefix, c(".bed", ".bim", ".fam"))
    if (!overwrite && all(file.exists(ld_files))) {
      message("Reusing LD-pruned PLINK files: ", config$qc_outputs$ld_prefix)
      ld_prefix <- config$qc_outputs$ld_prefix
    } else {
      # LD pruning
      ld_prefix <- plink_ld_prune(
        input = qc_result$qc_prefix,
        output = config$qc_outputs$ld_prefix,
        plink = config$plink$path,
        window = config$ld_pruning$window,
        step = config$ld_pruning$step,
        r2 = config$ld_pruning$r2
      )
    }
    save_pipeline_stage(list(ld_prefix = ld_prefix), config, "ld_pruning")
  }
})

# PCA module stage
run_stage("pca", {
  if (!overwrite && file.exists(config$pca$output_file)) {
    message("Reusing cached PCA stage.")
  } else {
    # Load QC and LD stages
    ld_result <- load_ld_pruning_stage(config)
    qc_result <- load_qc_stage(config)
    ld_obj <- plink_to_bigSNP(
      bed_file = paste0(ld_result$ld_prefix, ".bed"),
      overwrite = overwrite
    )
    # Calculate PCA
    pca_result <- run_pca_bigsnpr(
      obj = ld_obj,
      n_pcs = config$pca$n_pcs,
      output_file = config$pca$output_file,
      covariates_file = config$pca$covariates_file
    )
    rm(ld_obj)
    gc()
  }
})

# LMM module stage
run_stage("lmm", {
  # Define GEMMA outputs
  gemma_files <- file.path(
    config$gemma$output_dir,
    c("gemma_all_results.csv", "gemma_all_evaluation.csv", "gemma_best_strategy_by_var.csv")
  )
  if (!overwrite && all(file.exists(gemma_files))) {
    message("Reusing saved GEMMA results")
    save_pipeline_stage(list(output_files = gemma_files), config, "lmm")
  } else {
    # Load QC + kinship
    qc_result <- load_qc_stage(config)
    kinship_result <- load_kinship_stage(config)
    # Run GEMMA
    gemma_results <- run_gemma_lmm_all_variables(
      config = config,
      qc_prefix = qc_result$qc_prefix,
      climate_data = qc_result$climate_data,
      kinship_file = kinship_result$gemma_kinship_file,
      phenotypes = config$env$vars
    )
    save_pipeline_stage(list(output_files = gemma_files), config, "lmm")
  }
})

# LFMM module stage
run_stage("lfmm", {
  # Define LFMM outputs
  lfmm_files <- file.path(
    config$lfmm$output_dir,
    c("lfmm_all_results.csv", "lfmm_all_evaluation.csv", "lfmm_best_strategy_by_var.csv")
  )
  if (!overwrite && all(file.exists(lfmm_files))) {
    message("Reusing saved LFMM results")
    save_pipeline_stage(list(output_files = lfmm_files), config, "lfmm")
  } else {
    # Load QC
    qc_result <- load_qc_stage(config)
    qc_obj <- plink_to_bigSNP(
      bed_file = paste0(qc_result$qc_prefix, ".bed"),
      overwrite = overwrite
    )
    # Run LFMM
    lfmm_results <- run_lfmm_all_variables(
      config = config,
      geno = qc_obj$genotypes,
      map = qc_obj$map,
      fam = qc_obj$fam,
      climate_data = qc_result$climate_data,
      phenotypes = config$env$vars
    )
    save_pipeline_stage(list(output_files = lfmm_files), config, "lfmm")
    rm(qc_obj)
    gc()
  }
})

# RDA module stage
run_stage("rda", {
  # Define RDA outputs
  rda_files <- file.path(
    config$rda$output_dir,
    c("rda_all_results.csv", "rda_all_evaluation.csv", "rda_best_strategy_by_var.csv")
  )
  if (!overwrite && all(file.exists(rda_files))) {
    message("Reusing saved RDA results.")
    save_pipeline_stage(list(output_files = rda_files), config, "rda")
  } else {
    # Load QC + LD stages
    qc_result <- load_qc_stage(config)
    if (isTRUE(config$rda$use_full_qc)) {
      rda_prefix <- qc_result$qc_prefix
    } else {
      ld_result <- load_ld_pruning_stage(config)
      rda_prefix <- ld_result$ld_prefix
    }
    rda_obj <- plink_to_bigSNP(
      bed_file = paste0(rda_prefix, ".bed"),
      overwrite = overwrite
    )
    # Run RDA
    rda_results <- run_rda_all_variables(
      config = config,
      geno = rda_obj$genotypes,
      map = rda_obj$map,
      fam = rda_obj$fam,
      climate_data = qc_result$climate_data,
      phenotypes = config$env$vars
    )
    save_pipeline_stage(list(output_files = rda_files), config, "rda")
    rm(rda_obj)
    gc()
  }
})

# pcadapt module stage
run_stage("pcadapt", {
  # Define pcadapt module outputs
  pcadapt_files <- file.path(
    config$pcadapt$output_dir,
    c("pcadapt_all_results.csv", "pcadapt_all_evaluation.csv", "pcadapt_best_strategy.csv")
  )
  if (!overwrite && all(file.exists(pcadapt_files))) {
    message("Reusing saved pcadapt results.")
    save_pipeline_stage(list(output_files = pcadapt_files), config, "pcadapt")
  } else {
    # Load QC stage
    qc_result <- load_qc_stage(config)
    # Run pcadapt
    pcadapt_results <- run_pcadapt_strategies(
      config = config,
      qc_prefix = qc_result$qc_prefix
    )
    save_pipeline_stage(list(output_files = pcadapt_files), config, "pcadapt")
  }
})

# SNP overlap module stage
run_stage("snp_overlap", {
  if (!overwrite && stage_cache_exists("snp_overlap")) {
    message("Reusing cached SNP-overlap stage.")
  } else {
    # Run SNP overlap for all vars
    overlap_results <- run_snp_overlap_all_variables(
      config = config,
      gemma_results = load_gemma_results(config),
      lfmm_results = load_lfmm_results(config),
      rda_results = load_rda_results(config),
      pcadapt_results = load_pcadapt_results(config)
    )
    save_pipeline_stage(overlap_results, config, "snp_overlap")
  }
})

# Consensus set module stage
run_stage("consensus", {
  if (!overwrite && consensus_files_exist(config, config$env$vars)) {
    message("Reusing saved consensus results.")
    consensus_results <- list(primary_set = validate_primary_set(config$consensus$primary_set))
  } else {
    # Run consensus evaluation for all vars
    consensus_results <- run_consensus_all_variables(
      config = config,
      gemma_results = load_gemma_results(config),
      lfmm_results = load_lfmm_results(config),
      rda_results = load_rda_results(config),
      pcadapt_results = load_pcadapt_results(config),
      phenotypes = config$env$vars
    )
  }
  save_pipeline_stage(
    list(
      output_dir = config$consensus$output_dir,
      primary_set = consensus_results$primary_set
    ),
    config, "consensus"
  )
})

# LD processing module stage
run_stage("consensus_ld", {
  if (!overwrite && ld_block_files_exist(config, config$env$vars)) {
    message("Reusing saved post-hoc LD results.")
  } else {
    # Load QC + consensus set stage
    qc_result <- load_qc_stage(config)
    consensus_results <- load_consensus_results(config)
    # Run LD processing
    ld_results <- run_ld_blocks_all_variables(
      config = config,
      qc_prefix = qc_result$qc_prefix,
      consensus_results = consensus_results,
      phenotypes = config$env$vars
    )
  }
  save_pipeline_stage(
    list(output_dir = config$ld_pruning_ph$output_dir),
    config, "consensus_ld"
  )
})

# Primary SNP set module stage
run_stage("primary_snps", {
  primary_files <- file.path(
    config$ld_pruning_ph$output_dir,
    c("primary_ld_summary.csv", "final_primary_ld_pruned_snps.csv")
  )
  if (!overwrite && all(file.exists(primary_files)) && primary_snp_files_match_config(config)) {
    message("Reusing saved primary SNP results.")
    primary <- load_primary_snp_results(config)
  } else {
    # Create primary SNP set results
    primary <- create_primary_snp_results(
      config = config,
      ld_results = load_ld_block_results(config),
      phenotypes = config$env$vars
    )
  }
  save_pipeline_stage(
    list(set = primary$set, output_files = primary_files),
    config, "primary_snps"
  )
})

# Adaptive scoring module stage
run_stage("adaptive_scoring", {
  # Define output files
  score_files <- file.path(config$adaptive_scoring$output_dir, c(
    "adaptive_snp_direction_table.csv",
    "adaptive_snp_direction_summary.csv",
    "accession_directional_scores.csv",
    "top_50_directionally_extreme_accessions_by_variable.csv",
    "top_50_higher_direction_accessions_by_variable.csv",
    "top_50_lower_direction_accessions_by_variable.csv",
    "directional_score_summary_by_variable.csv"
  ))
  if (!overwrite && all(file.exists(score_files))) {
    message("Reusing saved adaptive-scoring results.")
  } else {
    # Load QC + primary set stages
    qc_result <- load_qc_stage(config)
    primary <- load_primary_snp_results(config)
    # Run adaptive scoring
    adaptive_results <- run_adaptive_germplasm_scoring(
      primary_lead_snps = primary$lead_snps,
      gemma_results = load_gemma_results(config),
      rda_results = load_rda_results(config),
      lfmm_results = load_lfmm_results(config),
      qc_prefix = qc_result$qc_prefix,
      metadata = qc_result$metadata,
      output_dir = config$adaptive_scoring$output_dir,
      sample_col = config$metadata$sample_col,
      overwrite = overwrite
    )
  }
  save_pipeline_stage(
    list(primary_set = load_primary_snp_results(config)$set, output_files = score_files),
    config, "adaptive_scoring"
  )
})

# Biological interpretation module stage
run_stage("bio_interp", {
  bio_file <- file.path(
    config$biological_interpretation$output_dir,
    "mapped_candidate_genes_annotated.rds"
  )
  if (!overwrite && file.exists(bio_file)) {
    message("Reusing saved biological-interpretation results.")
  } else {
    # Load primary set stage
    primary <- load_primary_snp_results(config)
    # Run biological interpretation
    biointerp_df <- run_biointerpretation_workflow(
      config = config,
      primary_lead_snps = primary$lead_snps,
      overwrite_annotation = overwrite
    )
  }
  save_pipeline_stage(
    list(primary_set = load_primary_snp_results(config)$set, output_file = bio_file),
    config, "bio_interp"
  )
})

message("\nPipeline stages complete.\n
        Render the report with:")
message("  Rscript -e \"rmarkdown::render('Pipeline.Rmd')\"")
