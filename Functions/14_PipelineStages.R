# ---- Pipeline Stages ----

# pipeline_stage_dir()
# Defines where the stage .rds files are stored
pipeline_stage_dir <- function(config) {
  if (!is.null(config$pipeline$stage_dir)) {
    config$pipeline$stage_dir
  } else {
    "Output/PipelineStages"
  }
}

# pipeline_stage_file()
# Construts .rds file name for each stage
pipeline_stage_file <- function(config, stage) {
  file.path(pipeline_stage_dir(config), paste0(stage, ".rds"))
}

# save_pipeline_stage()
# Saves stage result as .rds output
save_pipeline_stage <- function(object, config, stage) {
  path <- pipeline_stage_file(config, stage)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  message("Saved stage result: ", path)
  invisible(object)
}

# load_pipeline_stage()
# Load a preexisting stage
load_pipeline_stage <- function(config, stage, required = TRUE) {
  path <- pipeline_stage_file(config, stage)
  
  if (!file.exists(path)) {
    if (required) {
      stop("Stage result not found: ", path,
           ". Run: Rscript Scripts/run_pipeline.R --stages ", stage)
    }
    return(NULL)
  }
  
  readRDS(path)
}

# pipeline_stage_enabled()
# Determines if a stage should be run according to config or user defined stages
pipeline_stage_enabled <- function(config, key, selected_stages = NULL) {
  if (!is.null(selected_stages)) {
    return(key %in% selected_stages)
  }
  
  isTRUE(config$analysis[[paste0("run_", key)]])
}

# pipeline_require_files()
# Checks all required files exist for a stage
pipeline_require_files <- function(paths, stage) {
  missing <- paths[!file.exists(paths)]
  
  if (length(missing) > 0) {
    stop("Stage '", stage, "' requires missing file(s):\n- ",
         paste(missing, collapse = "\n- "))
  }
  
  invisible(paths)
}

# load_environment_stage()
# Reconstructs environmental stage object using saved analysis output files
load_environment_stage <- function(config) {
  # Check for existing stage .rds file
  path <- pipeline_stage_file(config, "env")
  if (file.exists(path)) {
    return(readRDS(path))
  }
  
  climate_pca_file <- "Output/Structure/climate_pca.rds"
  
  pipeline_require_files(c(config$env$output_file, climate_pca_file), "env")
  
  metadata <- if (isTRUE(config$subsetting$use)) {
    read.csv(config$subsetting$sample_metadata, stringsAsFactors = FALSE)
  } else {
    read.csv(config$metadata$file, stringsAsFactors = FALSE)
  }
  
  list(
    metadata = metadata,
    climate_data = read.csv(config$env$output_file, check.names = FALSE),
    climate_pca = readRDS(climate_pca_file)
  )
}

# load_qc_stage()
# Reconstructs quality control stage object using saved analysis output files
load_qc_stage <- function(config) {
  # Check for existing stage .rds file
  path <- pipeline_stage_file(config, "qc")
  if (file.exists(path)) {
    return(readRDS(path))
  }
  
  marker_prefix <- config$marker_subsetting$output_prefix
  marker_files <- paste0(marker_prefix, c(".bed", ".bim", ".fam"))
  
  qc_prefix <- if (isTRUE(config$marker_subsetting$use_extract) && all(file.exists(marker_files))) {
    marker_prefix
  } else {
    config$qc_outputs$qc_prefix
  }
  pipeline_require_files(
    c(paste0(qc_prefix, c(".bed", ".bim", ".fam")), config$metadata$filtered_file),
    "qc"
  )
  
  climate_data <- load_environment_stage(config)$climate_data
  fam <- read.table(paste0(qc_prefix, ".fam"), stringsAsFactors = FALSE)
  sample_col <- config$metadata$sample_col
  climate_match <- match(as.character(fam[, 2]), as.character(climate_data[[sample_col]]))
  
  if (anyNA(climate_match)) {
    stop("Climate data not available for all saved QC samples")
  }
  
  list(
    qc_prefix = qc_prefix,
    metadata = readRDS(config$metadata$filtered_file),
    climate_data = climate_data[climate_match, , drop = FALSE],
    dimensions = count_plink_dataset(qc_prefix),
    marker_summary = data.frame()
  )
}

# load_kinship_stage()
# Reconstructs kinship estimation stage object using saved analysis output files
load_kinship_stage <- function(config) {
  # Check for existing stage .rds file
  path <- pipeline_stage_file(config, "kinship")
  if (file.exists(path)) {
    return(readRDS(path))
  }
  
  kinship_file <- file.path(config$gemma$output_dir, paste0(config$gemma$kinship_prefix, ".cXX.txt"))
  pipeline_require_files(kinship_file, "kinship")
  
  list(gemma_kinship_file = kinship_file)
}

# load_ld_pruning_stage()
# Reconstructs ld pruning stage object using saved analysis output files
load_ld_pruning_stage <- function(config) {
  # Check for existing stage .rds file
  path <- pipeline_stage_file(config, "ld_pruning")
  if (file.exists(path)) {
    return(readRDS(path))
  }
  
  ld_prefix <- config$qc_outputs$ld_prefix
  pipeline_require_files(paste0(ld_prefix, c(".bed", ".bim", ".fam")), "ld_pruning")
  
  list(ld_prefix = ld_prefix)
}

# load_pca_stage()
# Reconstructs PCA stage object using saved analysis output files
load_pca_stage <- function(config) {
  # Check for existing stage .rds file
  if (file.exists(config$pca$output_file)) {
    return(readRDS(config$pca$output_file))
  }
  
  load_pipeline_stage(config, "pca")
}

# load_overlap_stage()
# Reconstructs SNP overlap stage object using saved analysis output files
load_overlap_stage <- function(config) {
  load_pipeline_stage(config, "snp_overlap")
}

# load_biointerpretation_results()
# Reconstructs biological interpretation stage object using saved analysis output files
load_biointerpretation_results <- function(config) {
  path <- file.path(config$biological_interpretation$output_dir, "mapped_candidate_genes_annotated.rds")
  
  if (!file.exists(path)) {
    stop("Saved biological interpretation results were not found")
  }
  
  readRDS(path)
}

# load_candidate_gene_summary()
# Reconstructs candidate gene summary stage object using saved analysis output files
load_candidate_gene_summary <- function(config) {
  path <- file.path(config$biological_interpretation$output_dir, "candidate_gene_summary.rds")
  if (!file.exists(path)) {
    stop("Saved candidate gene summary was not found")
  }
  readRDS(path)
}

# update_pipeline_manifest()
# Tracks status of the latest execution of the pipeline (completed stages, execution time)
update_pipeline_manifest <- function(config, stage, status, started, finished) {
  dir.create(pipeline_stage_dir(config), recursive = TRUE, showWarnings = FALSE)
  path <- file.path(pipeline_stage_dir(config), "pipeline_manifest.csv")
  
  current <- if (file.exists(path)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    data.frame(
      stage = character(), status = character(), started = character(),
      finished = character(), elapsed_minutes = numeric()
    )
  }

  row <- data.frame(
    stage = stage,
    status = status,
    started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    finished = format(finished, "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_minutes = round(as.numeric(difftime(finished, started, units = "mins")), 3)
  )
  
  current <- current[current$stage != stage, , drop = FALSE]
  write.csv(rbind(current, row), path, row.names = FALSE)
  invisible(row)
}
