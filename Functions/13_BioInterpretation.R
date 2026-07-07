# ---- Biological Interpretation ----

# download_wheat_annotation()
# Downloads wheat annotation RefSeq v2.1 GFF3 annotation file
# annotation_file = File name for annotation file
# overwrite = defaults to FALSE, reuses existing annotation file
# Returns: Annotation file
download_wheat_annotation <- function(
    annotation_file = "RawData/IWGSC_RefSeq_v2.1_annotation.gff3.gz",
    overwrite = FALSE
) {
  
  dir.create(dirname(annotation_file), recursive = TRUE, showWarnings = FALSE)
  
  # Reuse existing file
  if (file.exists(annotation_file) && !overwrite) {
    message("Using existing annotation file: ", annotation_file)
    return(annotation_file)
  }
  
  # Download from url
  annotation_url <- paste0(
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/gff3/",
    "triticum_aestivum_refseqv2/",
    "Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.63.gff3.gz"
  )
  message("Downloading wheat gene annotation")
  download.file(url = annotation_url, destfile = annotation_file, mode = "wb")
  
  return(annotation_file)
}

# map_candidate_snps_to_genes()
# Maps selected primary lead SNPs to nearby genes using chromosome and position
# robustness_results = Robustness results object
# gene_annotation_file = Gene annotation file path
# output_dir = Biological interpretation output directory
# flank_bp = Flanking window, number of base pairs upstream and downstream of each gene to include
# overwrite_annotation = defaults to FALSE, uses existing annotation file
# Returns: A data frame of mapped SNP-gene relationships
map_candidate_snps_to_genes <- function(
    robustness_results,
    gene_annotation_file = "RawData/IWGSC_RefSeq_v2.1_annotation.gff3.gz",
    output_dir = "Output/BioInterpretation",
    flank_bp = 10000,
    overwrite_annotation = FALSE
) {
  
  message("\nMapping candidate SNPs to genes")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Download annotation if missing
  gene_annotation_file <- download_wheat_annotation(
    annotation_file = gene_annotation_file,
    overwrite = overwrite_annotation
  )
  
  # Extract candidate SNPs
  candidate_snps <- robustness_results$primary_lead_snps
  
  if (is.null(candidate_snps) || nrow(candidate_snps) == 0) {
    stop("No candidate SNPs found in robustness_results$primary_lead_snps")
  }
  
  # Format candidate SNP data
  candidate_snps <- candidate_snps %>%
    filter(
      selected_primary == TRUE,
      is_lead == TRUE
    ) %>%
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
      !is.na(position)
    ) %>%
    distinct()
  
  if (nrow(candidate_snps) == 0) {
    stop("No selected primary lead SNPs with chr and position were found.")
  }
  
  message("Candidate SNPs with genomic coordinates: ", nrow(candidate_snps))
  
  # Convert SNPs to GRanges
  snp_gr <- GenomicRanges::GRanges(
    seqnames = candidate_snps$chr,
    ranges = IRanges::IRanges(
      start = candidate_snps$position,
      end = candidate_snps$position
    )
  )
  GenomicRanges::mcols(snp_gr) <- S4Vectors::DataFrame(candidate_snps)
  
  # Read gene annotation
  message("Reading gene annotation")
  
  genes <- rtracklayer::import(gene_annotation_file)
  
  # Keep only gene features
  if ("type" %in% names(GenomicRanges::mcols(genes))) {
    genes <- genes[GenomicRanges::mcols(genes)$type == "gene"]
  }
  if (length(genes) == 0) {
    stop("No gene features found in the annotation file")
  }
  
  # Clean gene chromosome names
  GenomeInfoDb::seqlevels(genes) <- gsub("^chr", "",
                                         GenomeInfoDb::seqlevels(genes),
                                         ignore.case = TRUE)
  gene_info <- as.data.frame(GenomicRanges::mcols(genes))
  
  # Extract gene ID
  if ("ID" %in% names(gene_info)) {
    gene_id <- as.character(gene_info$ID)
  } else if ("gene_id" %in% names(gene_info)) {
    gene_id <- as.character(gene_info$gene_id)
  } else {
    gene_id <- as.character(names(genes))
  }
  
  # Extract gene name, if available
  if ("Name" %in% names(gene_info)) {
    gene_name <- as.character(gene_info$Name)
  } else if ("gene_name" %in% names(gene_info)) {
    gene_name <- as.character(gene_info$gene_name)
  } else {
    gene_name <- NA_character_
  }
  
  # Extract gene description, if available
  if ("description" %in% names(gene_info)) {
    gene_description <- as.character(gene_info$description)
  } else if ("Note" %in% names(gene_info)) {
    gene_description <- as.character(gene_info$Note)
  } else {
    gene_description <- NA_character_
  }
  
  # Clean gene IDs for joining with other resources
  gene_id_clean <- gene_id
  gene_id_clean <- gsub("^gene:", "", gene_id_clean)
  gene_id_clean <- gsub("\\.\\d+$", "", gene_id_clean)
  GenomicRanges::mcols(genes)$gene_id <- gene_id
  GenomicRanges::mcols(genes)$gene_id_clean <- gene_id_clean
  GenomicRanges::mcols(genes)$gene_name <- gene_name
  GenomicRanges::mcols(genes)$gene_description <- gene_description
  
  # Expand gene windows
  gene_windows <- genes
  original_gene_start <- GenomicRanges::start(genes)
  original_gene_end <- GenomicRanges::end(genes)
  GenomicRanges::start(gene_windows) <- pmax(1, original_gene_start - flank_bp)
  GenomicRanges::end(gene_windows) <- original_gene_end + flank_bp
  GenomicRanges::mcols(gene_windows)$gene_start <- original_gene_start
  GenomicRanges::mcols(gene_windows)$gene_end <- original_gene_end
  
  # Check chromosome overlap
  snp_chr <- unique(as.character(GenomicRanges::seqnames(snp_gr)))
  gene_chr <- unique(as.character(GenomicRanges::seqnames(gene_windows)))
  common_chr <- intersect(snp_chr, gene_chr)
  
  message("Candidate SNP chromosomes: ", paste(head(snp_chr, 15), collapse = ", "))
  message("Gene chromosomes: ", paste(head(gene_chr, 15), collapse = ", "))
  message("Common chromosomes: ", length(common_chr))
  
  if (length(common_chr) == 0) {
    warning("No chromosome names match between candidate SNPs and genes")
  }
  
  # Map SNPs to genes
  message("Mapping SNPs to genes with +/- ", flank_bp, " bp window")
  
  hits <- GenomicRanges::findOverlaps(
    query = snp_gr,
    subject = gene_windows,
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0) {
    warning("No SNPs mapped to genes using the selected flanking window.")
    mapped_candidate_genes <- tibble()
  } else {
    # Extract hits
    snp_hits <- snp_gr[S4Vectors::queryHits(hits)]
    gene_hits <- gene_windows[S4Vectors::subjectHits(hits)]
    snp_df <- as.data.frame(GenomicRanges::mcols(snp_hits))
    gene_df <- as.data.frame(GenomicRanges::mcols(gene_hits))
    
    # Create data frame of mapped SNPs with gene information
    mapped_candidate_genes <- bind_cols(
      snp_df,
      tibble(
        gene_chr = as.character(GenomicRanges::seqnames(gene_hits)),
        gene_start = gene_df$gene_start,
        gene_end = gene_df$gene_end,
        gene_window_start = GenomicRanges::start(gene_hits),
        gene_window_end = GenomicRanges::end(gene_hits),
        gene_id = gene_df$gene_id,
        gene_id_clean = gene_df$gene_id_clean,
        gene_name = gene_df$gene_name,
        gene_description = gene_df$gene_description
      )
    ) %>%
      mutate(
        distance_to_gene_bp = case_when(
          position >= gene_start & position <= gene_end ~ 0,
          position < gene_start ~ gene_start - position,
          position > gene_end ~ position - gene_end,
          TRUE ~ NA_real_
        ),
        snp_gene_position = case_when(
          position >= gene_start & position <= gene_end ~ "gene_body",
          position < gene_start ~ "upstream_or_before_gene",
          position > gene_end ~ "downstream_or_after_gene",
          TRUE ~ NA_character_
        ),
        flank_bp = flank_bp
      ) %>%
      dplyr::select(
        phenotype,
        marker,
        chr,
        position,
        n_methods,
        methods,
        min_p,
        min_q,
        consensus_set,
        gene_chr,
        gene_start,
        gene_end,
        gene_window_start,
        gene_window_end,
        gene_id,
        gene_id_clean,
        gene_name,
        gene_description,
        distance_to_gene_bp,
        snp_gene_position,
        flank_bp
      ) %>%
      arrange(
        phenotype,
        marker,
        distance_to_gene_bp
      )
  }
  
  write_csv(mapped_candidate_genes, file.path(output_dir, "mapped_candidate_genes.csv"))
  saveRDS(mapped_candidate_genes, file.path(output_dir, "mapped_candidate_genes.rds"))
  
  message("Mapped SNP-gene rows: ", nrow(mapped_candidate_genes))
  
  return(mapped_candidate_genes)
}


# get_biomart_gene_annotation()
# Retrieves gene annotations from Ensembl Plants BioMart for wheat RefSeq v2.1 genes
# gene_ids = Character vector of wheat gene IDs
# Returns: A table with gene ID, BioMart gene name, BioMart description, and gene biotype
get_biomart_gene_annotation <- function(gene_ids) {
  
  gene_ids <- unique(gene_ids)
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  gene_ids <- gsub("^gene:", "", gene_ids)
  gene_ids <- gsub("\\.\\d+$", "", gene_ids)
  
  # No genes
  if (length(gene_ids) == 0) {
    return(tibble())
  }
  
  message("Retrieving gene annotations from BioMart RefSeq v2.1")
  
  mart <- useEnsemblGenomes(
    biomart = "plants_mart",
    dataset = "tarefseqv2_eg_gene"
  )
  
  # Retrieve annotations
  annotations <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "description",
      "gene_biotype"
    ),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart
  )
  
  if (nrow(annotations) == 0) {
    warning("BioMart returned 0 annotations.")
    return(tibble())
  }
  
  # Format annotations tibble
  annotations %>%
    as_tibble() %>%
    dplyr::rename(
      gene_id_clean = ensembl_gene_id,
      biomart_gene_name = external_gene_name,
      biomart_description = description
    ) %>%
    mutate(
      biomart_gene_name = as.character(biomart_gene_name),
      biomart_description = as.character(biomart_description),
      gene_biotype = as.character(gene_biotype)
    ) %>%
    distinct(gene_id_clean, .keep_all = TRUE)
}

# add_biomart_annotation()
# Adds BioMart gene annotation to the mapped SNP-gene table.
# mapped_candidate_genes = Data frame returned by map_candidate_snps_to_genes().
# Returns: Mapped SNP-gene data frame with BioMart annotation columns added
add_biomart_annotation <- function(mapped_candidate_genes) {
  
  # Retrieve BioMart annotations
  biomart_annotation <- get_biomart_gene_annotation(
    mapped_candidate_genes$gene_id_clean
  )
  
  if (nrow(biomart_annotation) > 0) {
    
    mapped_candidate_genes <- mapped_candidate_genes %>%
      left_join(biomart_annotation, by = "gene_id_clean") %>%
      mutate(
        gene_name = coalesce(na_if(gene_name, ""), na_if(biomart_gene_name, "")),
        gene_description = coalesce(na_if(gene_description, ""),
                                    na_if(biomart_description, ""))
      ) %>%
      dplyr::select(
        -any_of(c(
          "biomart_gene_name",
          "biomart_description"
        ))
      )
  }
  
  # Include column even if no BioMart annotations
  if (!"gene_biotype" %in% names(mapped_candidate_genes)) {
    mapped_candidate_genes$gene_biotype <- NA_character_
  }
  
  mapped_candidate_genes
}

# read_tf_annotation()
# Reads a local wheat transcription factor annotation file
# tf_annotation_file = Local TF annotation CSV file path
# Returns: Tibble with gene ID, TF gene ID, TF gene name, and TF family
read_tf_annotation <- function(tf_annotation_file) {
  
  # Import annotation file
  tf_annotation <- read_csv(tf_annotation_file, show_col_types = FALSE)
  
  tf_annotation <- tf_annotation[, 1:2]
  colnames(tf_annotation) <- c("tf_gene_id", "tf_gene_name")
  
  # Format TF annotation
  tf_annotation %>%
    mutate(
      tf_gene_id = as.character(tf_gene_id),
      tf_gene_name = as.character(tf_gene_name),
      gene_id_clean = gsub("^gene:", "", tf_gene_id),
      gene_id_clean = gsub("\\.\\d+$", "", gene_id_clean),
      tf_family = sub("_.*", "", tf_gene_name)
    ) %>%
    dplyr::select(
      gene_id_clean,
      tf_gene_id,
      tf_gene_name,
      tf_family
    ) %>%
    distinct(gene_id_clean, .keep_all = TRUE)
}

# add_tf_annotation()
# Adds transcription factor annotation information to mapped candidate genes.
# mapped_candidate_genes = Mapped SNP-gene data frame
# tf_annotation_file = Local TF annotation CSV file path
# Returns: Mapped SNP-gene data frame with TF annotation
add_tf_annotation <- function(mapped_candidate_genes,
                              tf_annotation_file = NULL) {
  
  if (!"gene_biotype" %in% names(mapped_candidate_genes)) {
    mapped_candidate_genes$gene_biotype <- NA_character_
  }
  
  if (!"gene_description" %in% names(mapped_candidate_genes)) {
    mapped_candidate_genes$gene_description <- NA_character_
  }
  
  # if no TF annotation available
  if (is.null(tf_annotation_file) || is.na(tf_annotation_file) || !file.exists(tf_annotation_file)) {
    message("No TF annotation file provided. Skipping TF annotation.")
    return(
      mapped_candidate_genes %>%
        mutate(
          tf_gene_id = NA_character_,
          tf_gene_name = NA_character_,
          tf_family = NA_character_,
          annotation_level = case_when(
            !is.na(gene_description) & gene_description != "" ~ "gene_description",
            !is.na(gene_biotype) & gene_biotype != "" ~ "biotype_only",
            TRUE ~ "candidate_region_only"
          )
        )
    )
  }
  
  message("Adding local TF annotation")
  tf_annotation <- read_tf_annotation(tf_annotation_file)
  mapped_candidate_genes %>%
    left_join(tf_annotation, by = "gene_id_clean") %>%
    mutate(
      annotation_level = case_when(
        !is.na(tf_gene_name) & tf_gene_name != "" ~ "tf_family_annotation",
        !is.na(gene_description) & gene_description != "" ~ "gene_description",
        !is.na(gene_biotype) & gene_biotype != "" ~ "biotype_only",
        TRUE ~ "candidate_region_only"
      )
    )
}

# annotate_candidate_genes()
# Adds optional annotation layers to mapped SNP-gene table (BioMart abd TF)
# mapped_candidate_genes = Mapped SNP-gene data frame
# tf_annotation_file = Optional path to local TF annotation file
# use_biomart = if TRUE, adds BioMart gene annotation
# Returns: Annotated SNP-gene mapping data frame.
annotate_candidate_genes <- function(
    mapped_candidate_genes,
    tf_annotation_file = NULL,
    use_biomart = TRUE
) {
  
  annotated_genes <- mapped_candidate_genes
  
  # BioMart annotation
  if (isTRUE(use_biomart)) {
    annotated_genes <- add_biomart_annotation(annotated_genes)
  }
  
  if (!"gene_biotype" %in% names(annotated_genes)) {
    annotated_genes$gene_biotype <- NA_character_
  }
  
  # TF annotation
  annotated_genes <- add_tf_annotation(
    mapped_candidate_genes = annotated_genes,
    tf_annotation_file = tf_annotation_file
  )
  
  annotated_genes
}

# summarize_candidate_genes()
# Creates a gene-level summary from annotated SNP-gene mapping table
# mapped_candidate_genes = Annotated SNP-gene mapping data frame
# Returns: Gene-level summary table
summarize_candidate_genes <- function(mapped_candidate_genes) {
  
  if (nrow(mapped_candidate_genes) == 0) {
    return(tibble())
  }
  
  # Optional annotation columns
  optional_cols <- c("gene_biotype", "tf_gene_name", "tf_family", "annotation_level")
  
  for (col in optional_cols) {
    if (!col %in% names(mapped_candidate_genes)) {
      mapped_candidate_genes[[col]] <- NA
    }
  }
  
  # Summarize mapped candidate genes
  mapped_candidate_genes %>%
    group_by(
      phenotype,
      gene_id_clean,
      gene_name,
      gene_description,
      gene_biotype,
      tf_gene_name,
      tf_family,
      annotation_level
    ) %>%
    arrange(distance_to_gene_bp, min_q, min_p, .by_group = TRUE) %>%
    summarise(
      closest_lead_snp = dplyr::first(marker),
      closest_lead_snp_position = dplyr::first(paste0(chr, ":", position)),
      closest_snp_gene_position = dplyr::first(snp_gene_position),
      n_lead_snps = n_distinct(marker),
      all_lead_snps = paste(unique(marker), collapse = "; "),
      max_method_support = max(n_methods, na.rm = TRUE),
      methods_supported = paste(unique(methods), collapse = "; "),
      min_p = min(min_p, na.rm = TRUE),
      min_q = min(min_q, na.rm = TRUE),
      min_distance_to_gene_bp = min(distance_to_gene_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      phenotype,
      min_distance_to_gene_bp,
      desc(n_lead_snps),
      desc(max_method_support)
    )
}

# run_biointerpretation_workflow()
# Full biological interpretation workflow, SNP-to-genen mapping + annotations
# config = Pipeline configuration object
# robustness_results = Results object from robustness module
# overwrite_annotation = defaults to FALSE, reuses existing annotation file
# Returns: Annotated SNP-gene mapping data frame
run_biointerpretation_workflow <- function(
    config,
    robustness_results,
    overwrite_annotation = FALSE
) {
  
  # Configuration settings
  gene_annotation_file <- config$biological_interpretation$gene_annotation_file
  output_dir <- config$biological_interpretation$output_dir
  flank_bp <- config$biological_interpretation$flank_bp
  use_biomart <- isTRUE(config$biological_interpretation$use_biomart)
  tf_annotation_file <- config$biological_interpretation$tf_annotation_file
  
  if (is.null(tf_annotation_file) || length(tf_annotation_file) == 0 || is.na(tf_annotation_file)) {
    tf_annotation_file <- NULL
  }
  
  # SNP-to-gene mapping
  biointerp_df <- map_candidate_snps_to_genes(
    robustness_results = robustness_results,
    gene_annotation_file = gene_annotation_file,
    output_dir = output_dir,
    flank_bp = flank_bp,
    overwrite_annotation = overwrite_annotation
  )
  
  # Optional annotation layers
  biointerp_df <- annotate_candidate_genes(
    mapped_candidate_genes = biointerp_df,
    tf_annotation_file = tf_annotation_file,
    use_biomart = use_biomart
  )
  
  # Optional gene-level summary
  candidate_gene_summary <- summarize_candidate_genes(biointerp_df)
  
  # Save final outputs
  write_csv(biointerp_df, file.path(output_dir, "mapped_candidate_genes_annotated.csv"))
  write_csv(candidate_gene_summary, file.path(output_dir, "candidate_gene_summary.csv"))
  
  saveRDS(biointerp_df, file.path(output_dir, "mapped_candidate_genes_annotated.rds"))
  saveRDS(candidate_gene_summary, file.path(output_dir, "candidate_gene_summary.rds"))
  
  message("Mapped SNP-gene rows: ", nrow(biointerp_df))
  message("Candidate gene summary rows: ", nrow(candidate_gene_summary))
  
  return(biointerp_df)
}