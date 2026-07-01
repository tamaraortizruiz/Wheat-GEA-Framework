# ---- Genotype Quality Control ----

# run_plink()
# Runs PLINK command from R using system2()
# plink = Path to the PLINK executable
# args = Character vector of PLINK command-line arguments
# Output: Runs PLINK and prints the command used
run_plink <- function(plink = "plink", args) {
  
  message("\nRunning PLINK:")
  message(plink, " ", paste(args, collapse = " "))
  
  status <- system2(command = plink, args = args)
  
  if (status != 0) {
    stop("PLINK failed. Check the .log file for details.")
  }
  
  invisible(status)
}

# create_keep_file()
# Creates a PLINK keep file using sample IDs in both metadata and PLINK .fam file
# metadata = Metadata data frame with sample information
# fam_file = Path to PLINK .fam file
# output_file = Path where the keep file will be saved
# sample_col = Column in metadata with sample IDs
# Output: PLINK keep file with family ID and sample ID
# Returns: Path to the generated keep file
create_keep_file <- function(
    metadata,
    fam_file,
    output_file,
    sample_col = "SeedID"
) {
  
  fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  # Set .fam column names
  colnames(fam) <- c(
    "family.ID",
    "sample.ID",
    "paternal.ID",
    "maternal.ID",
    "sex",
    "phenotype"
  )
  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  fam$sample.ID <- as.character(fam$sample.ID)
  
  keep_df <- fam %>%
    filter(sample.ID %in% metadata[[sample_col]]) %>%
    select(family.ID, sample.ID)
  
  if (nrow(keep_df) == 0) {
    stop("No matching samples between metadata subset and FAM file")
  }
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.table(keep_df, file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("Keep file saved to: ", output_file)
  message("Samples in keep file: ", nrow(keep_df))
  
  return(output_file)
}

# plink_filter()
# Performs genotype quality control using PLINK
# input = PLINK input prefix (without .bed/.bim/.fam)
# output = Output prefix for filtered PLINK files
# plink = Path to PLINK executable
# keep = Optional keep file for sample subsetting
# geno_na = Maximum NAs allowed per marker
# maf = Minimum minor allele frequency
# ind_na = Maximum NAs allowed per individual
# Output: Filtered PLINK binary files (.bed, .bim, .fam.)
# Returns: Output prefix of filtered PLINK dataset
plink_filter <- function(
    input,
    output,
    plink = "plink",
    keep = NULL,
    geno_na = 0.10,
    maf = 0.01,
    ind_na = NULL
) {
  
  # Check input files
  input_files <- paste0(input, c(".bed", ".bim", ".fam"))
  missing_input <- input_files[!file.exists(input_files)]
  
  if (length(missing_input) > 0) {
    stop("Missing PLINK input files:\n",
         paste(missing_input, collapse = "\n"))
  }
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  
  # PLINK arguments
  args <- c(
    "--bfile", input,
    "--allow-extra-chr"
  )
  
  # if sample filtering
  if (!is.null(keep)) {
    if (!file.exists(keep)) {
      stop("Sample keep file does not exist: ", keep)
    }
    args <- c(
      args,
      "--keep", keep
      )
  }
  
  args <- c(
    args,
    "--geno", as.character(geno_na),
    "--maf", as.character(maf),
    "--make-bed",
    "--out", output
  )
  
  # if filtering by missing individual data 
  if (!is.null(ind_na)) {
    args <- c(args, "--mind", as.character(ind_na))
  }
  
  # Run PLINK
  run_plink(plink = plink, args = args)
  
  # Check output
  output_files <- paste0(output, c(".bed", ".bim", ".fam"))
  missing_output <- output_files[!file.exists(output_files)]
  if (length(missing_output) > 0) {
    stop("Output files are missing:\n",
         paste(missing_output, collapse = "\n"))
  }
  
  message("PLINK filtering complete, saved to: ", output)
  
  return(output)
}

# plink_ld_prune(input, output, plink, window, step, r2)
# Performs LD pruning using PLINK --indep-pairwise
# input = PLINK input prefix
# output = Output prefix for LD-pruned PLINK files
# plink = Path to PLINK executable
# window = SNP window size for LD calculation
# step = Number of SNPs to shift the window each step
# r2 = LD threshold (SNPs with higher LD are pruned)
# Output: LD-pruned PLINK binary files and prune.in/prune.out files
# Returns: Output prefix of the LD-pruned dataset
plink_ld_prune <- function(
    input,
    output,
    plink = "plink",
    window = 50,
    step = 10,
    r2 = 0.2
) {
  
  # Check input
  input_files <- paste0(input, c(".bed", ".bim", ".fam"))
  missing_input <- input_files[!file.exists(input_files)]
  
  if (length(missing_input) > 0) {
    stop("Missing PLINK input files:\n",
         paste(missing_input, collapse = "\n"))
  }
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  prune_prefix <- paste0(output, "_prune")
  
  # Arguments for LD-pruned marker list
  args1 <- c(
    "--bfile", input,
    "--allow-extra-chr",
    "--indep-pairwise",
    as.character(window),
    as.character(step),
    as.character(r2),
    "--out", prune_prefix
  )
  
  run_plink(plink = plink, args = args1)
  
  # Extract LD-pruned markers
  args2 <- c(
    "--bfile", input,
    "--allow-extra-chr",
    "--extract", paste0(prune_prefix, ".prune.in"),
    "--make-bed",
    "--out", output
  )
  
  run_plink(plink = plink, args = args2)
  
  # Check output
  output_files <- paste0(output, c(".bed", ".bim", ".fam"))
  missing_output <- output_files[!file.exists(output_files)]
  if (length(missing_output) > 0) {
    stop("Output files are missing:\n",
         paste(missing_output, collapse = "\n"))
  }
  
  message("LD-pruned PLINK complete, saved to: ", output)
  
  return(output)
}

# filter_metadata()
# Filters and reorders metadata to match filtered PLINK .fam file
# metadata_file = Path to metadata .csv
# fam_file = Path to PLINK .fam file
# sample_col = Metadata column containing sample IDs
# output_file = Optional path to save filtered metadata as .rds
# Output: Metadata filtered and ordered to samples in .fam file
# Returns: Filtered metadata data frame
filter_metadata <- function(
    metadata_file,
    fam_file,
    sample_col = "SeedID",
    output_file = NULL
) {
  
  metadata <- read.csv(metadata_file)
  
  fam <- read.table(fam_file, header = FALSE)
  colnames(fam) <- c(
    "family.ID",
    "sample.ID",
    "paternal.ID",
    "maternal.ID",
    "sex",
    "phenotype"
  )
  
  # Check sample column
  if (!sample_col %in% colnames(metadata)) {
    stop("sample_col not found in metadata: ", sample_col)
  }
  
  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  fam$sample.ID <- as.character(fam$sample.ID)
  metadata_filtered <- metadata %>%
    filter(.data[[sample_col]] %in% fam$sample.ID)
  
  # Match genotype order
  metadata_filtered <- metadata_filtered[
    match(fam$sample.ID, metadata_filtered[[sample_col]]),
  ]
  
  if (any(is.na(metadata_filtered[[sample_col]]))) {
    stop("Some samples in the filtered fam file were not found in metadata.")
  }
  if (nrow(metadata_filtered) != nrow(fam)) {
    stop("Filtered metadata and fam file have different numbers of samples.")
  }
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(metadata_filtered, output_file)
    message("Filtered metadata saved to: ", output_file)
  }
  
  return(metadata_filtered)
}

# plink_extract_markers()
# Optionally extracts a selected list of markers from a PLINK dataset
# input = PLINK input prefix
# output = Output prefix for marker-filtered PLINK files
# plink = Path to PLINK executable
# extract_file = Text file containing marker IDs to keep
# overwrite = If FALSE, reuse existing output files
# Output: Marker-filtered PLINK .bed, .bim, and .fam files
# Returns: Output prefix of the marker-filtered dataset
plink_extract_markers <- function(
    input,
    output,
    plink = "plink",
    extract_file,
    overwrite = FALSE
) {
  
  output_files <- paste0(output, c(".bed", ".bim", ".fam"))
  
  if (!overwrite && all(file.exists(output_files))) {
    message("Reusing marker-filtered PLINK files: ", output)
    return(output)
  }
  
  if (!file.exists(extract_file)) {
    stop("Marker extract file does not exist: ", extract_file)
  }
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  
  # PLINK arguments
  args <- c(
    "--bfile", input,
    "--allow-extra-chr",
    "--extract", extract_file,
    "--make-bed",
    "--out", output
  )
  
  run_plink(plink = plink, args = args)
  
  message("Marker-filtered PLINK files saved to: ", output)
  
  return(output)
}

# plink_to_bigSNP()
# Converts a PLINK .bed file into a bigSNP object
# bed_file = Path to the PLINK .bed file
# overwrite = If TRUE, removes existing .rds and .bk files and remakes them
# Output: bigsnpr .rds and .bk files.
# Returns: Attached bigSNP object containing genotypes, map, and fam data
plink_to_bigSNP <- function(
    bed_file,
    overwrite = FALSE
) {
  
  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file)
  }
  
  rds_file <- sub("\\.bed$", ".rds", bed_file)
  bk_file  <- sub("\\.bed$", ".bk", bed_file)
  
  # if overwrite = TRUE, remove existing bigSNP .rds/.bk files
  if (overwrite) {
    if (file.exists(rds_file)) {
      file.remove(rds_file)
    }
    if (file.exists(bk_file)) {
      file.remove(bk_file)
    }
  }
  
  # if bigSNP files do not exist
  if (!file.exists(rds_file) || !file.exists(bk_file)) {
    snp_readBed(bed_file)
  }
  obj <- snp_attach(rds_file)
  
  return(obj)
}

# bigSNP_to_nummat()
# Converts a bigSNP genotype object into standard numeric matrix
# obj = bigSNP object created by bigsnpr
# Output: Numeric genotype matrix with samples as rows and markers as columns
# Returns: Genotype matrix
bigSNP_to_nummat <- function(
    obj
) {
  
  G <- obj$genotypes
  map <- obj$map
  fam <- obj$fam
  
  G_matrix <- as.matrix(G[])
  rownames(G_matrix) <- fam$sample.ID
  colnames(G_matrix) <- map$marker.ID
  
  return(G_matrix)
}



# plink_to_nummat()
# Converts PLINK files to bigSNP object to numeric genotype matrix
# plink_prefix = PLINK file prefix
# output_file = Optional path to save resulting object as .rds.
# overwrite_bigSNP = Whether to remake bigSNP files
# save_obj = Whether to keep the bigSNP object inside return list
# Output: Optional saved RDS containing genotype matrix, map, fam, and prefix
# Returns: List containing genotype matrix, map, fam, prefix (optionally obj)
plink_to_nummat <- function(
    plink_prefix,
    output_file = NULL,
    overwrite_bigSNP = TRUE,
    save_obj = FALSE
) {
  
  bed_file <- paste0(plink_prefix, ".bed")
  
  # PLINK to bigSNP
  obj <- plink_to_bigSNP(
    bed_file = bed_file,
    overwrite = overwrite_bigSNP
  )
  
  # bigSNP to numeric matrix
  G_matrix <- bigSNP_to_nummat(
    obj = obj
  )
  
  result <- list(
    obj = obj,
    G = G_matrix,
    map = obj$map,
    fam = obj$fam,
    prefix = plink_prefix
  )
  
  if (save_obj) {
    result$obj <- obj
  }
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, output_file)
    message("Numeric matrix object saved to: ", output_file)
  }
  
  message("Samples: ", nrow(G_matrix))
  message("Markers: ", ncol(G_matrix))
  
  return(result)
}

# count_plink_dataset()
# Returns number of samples and SNPs in PLINK dataset
# prefix = PLINK file prefix
# Returns: data frame of samples and SNPs
count_plink_dataset <- function(prefix) {
  
  fam <- read.table(paste0(prefix, ".fam"), header = FALSE, stringsAsFactors = FALSE)
  bim <- read.table(paste0(prefix, ".bim"), header = FALSE, stringsAsFactors = FALSE)
  
  data.frame(
    samples = nrow(fam),
    snps = nrow(bim)
  )
}

# read_plink_bim()
# Read PLINK BIM file for marker information
# prefix = PLINK QC genotype data prefix
# Output: Table of marker information from PLINK BIM file
read_plink_bim <- function(prefix) {
  
  bim_file <- paste0(prefix, ".bim")
  
  if (!file.exists(bim_file)) {
    stop("BIM file not found: ", bim_file)
  }
  
  bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bim) <- c(
    "chromosome",
    "marker.ID",
    "genetic_distance",
    "position_bp",
    "allele1",
    "allele2"
  )
  
  bim <- bim %>%
    mutate(
      chromosome = gsub("^chr", "", chromosome, ignore.case = TRUE),
      chromosome = as.character(chromosome),
      position_bp = as.numeric(position_bp)
    )
  
  return(bim)
}

# plot_marker_count_by_chr()
# Plot marker count by chromosome
# prefix = PLINK QC genotype data prefix
# output_file = Path to save plot
# Output: Returns plot
plot_marker_count_by_chr <- function(prefix, output_file = NULL) {
  
  bim <- read_plink_bim(prefix)
  chr_order <- sort(unique(bim$chromosome))
  
  # Marker counts
  marker_counts <- bim %>%
    count(chromosome, name = "n_markers") %>%
    mutate(chromosome = factor(chromosome, levels = chr_order))
  
  p <- ggplot(marker_counts, aes(x = chromosome, y = n_markers)) +
    geom_col() +
    theme_bw() +
    labs(title = "Marker count per chromosome",
         x = "Chromosome",
         y = "Number of markers") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_file, p, width = 9, height = 5, dpi = 300)
  }
  
  return(p)
}

# plot_marker_position_by_chr()
# Marker physical position by chromosome
# prefix = PLINK QC genotype data prefix
# output_file = Path to save plot
# Output: Returns plot
plot_marker_position_by_chr <- function(prefix, output_file = NULL) {
  
  bim <- read_plink_bim(prefix)
  chr_order <- sort(unique(bim$chromosome))
  
  bim_plot <- bim %>%
    filter(!is.na(position_bp), position_bp > 0) %>%
    mutate(
      chromosome = factor(chromosome, levels = chr_order),
      position_mb = position_bp / 1e6
    )
  
  p <- ggplot(bim_plot, aes(x = position_mb, y = chromosome)) +
    geom_point(alpha = 0.25, size = 0.4) +
    theme_bw() +
    labs(title = "Marker positions across chromosomes",
         x = "Physical position (Mb)",
         y = "Chromosome")
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  }
  
  return(p)
}

# plot_genotype_qc_overview()
# Produce all genotype QC plots and tables
# prefix = PLINK QC genotype data prefix
# output_dir = Directory path to write outputs to
# dataset_name = Name for output files
# Output: List object with QC plots and summary table
plot_genotype_qc_overview <- function(
    prefix,
    output_dir = "Output/QC",
    dataset_name = "QC_genotype"
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Marker count by chromosome
  p1 <- plot_marker_count_by_chr(
    prefix = prefix,
    output_file = file.path(output_dir, paste0(dataset_name, "_marker_count_by_chr.png"))
  )
  
  # Marker position by chromosome
  p2 <- plot_marker_position_by_chr(
    prefix = prefix,
    output_file = file.path(output_dir, paste0(dataset_name, "_marker_positions_by_chr.png"))
  )
  
  bim <- read_plink_bim(prefix)
  
  # Summary data frame grouped by chromosome
  summary_df <- bim %>%
    group_by(chromosome) %>%
    summarise(
      n_markers = n(),
      min_position_bp = min(position_bp, na.rm = TRUE),
      max_position_bp = max(position_bp, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(summary_df, file.path(output_dir, paste0(dataset_name, "_marker_summary.csv")),
            row.names = FALSE)
  
  return(list(
    marker_count_plot = p1,
    marker_position_plot = p2,
    marker_summary = summary_df
  ))
}