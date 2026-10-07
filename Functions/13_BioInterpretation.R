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

# download_wheat_functional_annotation()
# Downloads the URGI functional annotation ZIP and extracts the CSV file
# annotation_file = Destination functional annotation CSV path
# overwrite = if TRUE replaces the CSV, FALSE reuses an existing file
# Output: Writes the CSV and removes temporary download/extraction files
# Returns: Path to the functional annotation CSV
download_wheat_functional_annotation <- function(
    annotation_file = "RawData/iwgsc_refseqv2.1_functional_annotation.csv",
    overwrite = FALSE
) {
  
  if (file.exists(annotation_file) && !overwrite) {
    message("Using existing functional annotation: ", annotation_file)
    return(annotation_file)
  }
  
  dir.create(dirname(annotation_file), recursive = TRUE, showWarnings = FALSE)
  
  annotation_url <- paste0(
    "https://urgi.versailles.inrae.fr/download/iwgsc/",
    "IWGSC_RefSeq_Annotations/v2.1/",
    "iwgsc_refseqv2.1_functional_annotation.zip"
  )
  
  temporary_dir <- tempfile()
  dir.create(temporary_dir)
  on.exit(unlink(temporary_dir, recursive = TRUE), add = TRUE)
  
  zip_file <- file.path(temporary_dir, "annotation.zip")
  
  message("Downloading wheat functional annotation")
  download.file(annotation_url, zip_file, mode = "wb")
  
  utils::unzip(zip_file, exdir = temporary_dir)
  
  csv_file <- list.files(
    temporary_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(csv_file) != 1L) {
    stop("Expected one functional annotation CSV in the ZIP.")
  }
  
  if (!file.copy(csv_file, annotation_file, overwrite = overwrite)) {
    stop("Could not save functional annotation: ", annotation_file)
  }
  
  annotation_file
}

# clean_annotation_gene_id()
# Standardizes gene IDs for joins with local functional annotations
# x = Character vector of gene IDs
# Removes unnecessary gene prefix or suffic
# Returns: Character vector of cleaned gene IDs
clean_annotation_gene_id <- function(x) {
  sub("\\.\\d+$", "", sub("^gene:", "", trimws(as.character(x))))
}

# map_candidate_snps_to_genes()
# Maps selected primary lead SNPs to nearby genes using chromosome and position
# primary_lead_snps = Final LD-pruned SNP table from the configured primary set
# gene_annotation_file = Gene annotation file path
# output_dir = Biological interpretation output directory
# flank_bp = Flanking window, number of base pairs upstream and downstream of each gene to include
# overwrite_annotation = defaults to FALSE, uses existing annotation file
# Returns: A data frame of mapped SNP-gene relationships
map_candidate_snps_to_genes <- function(
    primary_lead_snps,
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
  candidate_snps <- primary_lead_snps
  
  if (is.null(candidate_snps) || nrow(candidate_snps) == 0) {
    stop("No candidate SNPs found in primary_lead_snps")
  }
  
  # Format candidate SNP data
  candidate_snps <- candidate_snps %>%
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

# read_local_functional_annotation()
# Reads full local URGI records
# functional_annotation_file = local CSV functional annotation file, 
# downloads the file if missing and the functional downloader is available
# Returns: Gene ID plus annotation type, domain, name and label columns
read_local_functional_annotation <- function(functional_annotation_file) {
  if (!file.exists(functional_annotation_file)) {
    if (!exists("download_wheat_functional_annotation", mode = "function")) {
      stop("Functional annotation file not found: ", functional_annotation_file)
    }
    functional_annotation_file <- download_wheat_functional_annotation(
      annotation_file = functional_annotation_file, overwrite = FALSE
    )
  }
  
  raw <- read_csv(
    functional_annotation_file,
    col_types = cols(.default = col_character()),
    name_repair = "check_unique", show_col_types = FALSE
  )
  required <- c("g2.identifier", "f.type", "f.domain", "f.name", "f.label")
  if (!all(required %in% names(raw))) {
    stop("Missing annotation columns: ", paste(setdiff(required, names(raw)), collapse = ", "))
  }
  out <- tibble::tibble(
    gene_id_clean = clean_annotation_gene_id(raw$g2.identifier),
    annotation_type = raw$f.type,
    annotation_domain = raw$f.domain,
    annotation_name = raw$f.name,
    annotation_label = raw$f.label
  )
  out <- dplyr::mutate(out, dplyr::across(dplyr::everything(), ~ dplyr::na_if(trimws(.x), "")))
  out <- dplyr::filter(out, !is.na(.data$gene_id_clean))
  dplyr::distinct(out)
}

# add_local_functional_annotation()
# Adds local descriptions to existing SNP-gene table
# mapped_candidate_genes = SNP-gene mapping table with gene_id_clean
# functional_annotation_file = Path to the local annotation CSV
# Returns: Mapping table with updated gene_description
add_local_functional_annotation <- function(mapped_candidate_genes, functional_annotation_file) {
  annotations <- read_local_functional_annotation(functional_annotation_file)
  descriptions <- annotations %>%
    dplyr::filter(.data$annotation_type == "Human readable description", !is.na(.data$annotation_name)) %>%
    dplyr::group_by(.data$gene_id_clean) %>%
    dplyr::summarise(local_description = paste(unique(.data$annotation_name), collapse = "; "), .groups = "drop")
  out <- dplyr::left_join(mapped_candidate_genes, descriptions, by = "gene_id_clean")
  if (!"gene_description" %in% names(out)) out$gene_description <- NA_character_
  out$gene_description <- dplyr::coalesce(out$local_description, dplyr::na_if(out$gene_description, ""))
  dplyr::select(out, -local_description)
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
# Adds BioMart gene annotation to the mapped SNP-gene table
# mapped_candidate_genes = Data frame returned by map_candidate_snps_to_genes()
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
# Adds optional annotation layers to mapped SNP-gene table
# mapped_candidate_genes = Mapped SNP-gene data frame
# tf_annotation_file = Optional path to local TF annotation file
# use_biomart = if TRUE, adds BioMart gene annotation
# functional_annotation_file = Optional local functional annotation CSV, takes priority over BioMart
# Returns: Annotated SNP-gene mapping data frame
annotate_candidate_genes <- function(
    mapped_candidate_genes,
    tf_annotation_file = NULL,
    use_biomart = TRUE,
    functional_annotation_file = NULL
) {
  
  annotated_genes <- mapped_candidate_genes
  
  # Skip annotation lookups if there are no mapped genes.
  if (nrow(annotated_genes) == 0) {
    return(annotated_genes)
  }
  
  has_local_file <- !is.null(functional_annotation_file) &&
    length(functional_annotation_file) == 1L &&
    !is.na(functional_annotation_file) &&
    nzchar(trimws(functional_annotation_file))
  
  if (has_local_file) {
    message("Adding local functional annotation")
    annotated_genes <- add_local_functional_annotation(
      mapped_candidate_genes = annotated_genes,
      functional_annotation_file = functional_annotation_file
    )
    
  } else if (isTRUE(use_biomart)) {
    annotated_genes <- tryCatch(
      {
        add_biomart_annotation(annotated_genes)
      },
      error = function(e) {
        warning(
          "BioMart annotation failed. Continuing with existing annotations. ",
          "Error: ", conditionMessage(e)
        )
        
        annotated_genes$biomart_status <- "biomart_failed"
        annotated_genes
      }
    )
    
  }
  
  # Maintain column expected by TF annotation and gene summaries.
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

# prepare_candidate_annotations()
# Collects every local candidate-gene record and adds TF annotation records
# candidate_gene_summary = Gene summary created after TF annotation.
# functional_annotation_file = Path to the local URGI annotation CSV.
# Returns: Distinct phenotype-gene-annotation records with their source
prepare_candidate_annotations <- function(candidate_gene_summary, functional_annotation_file) {
  if (nrow(candidate_gene_summary) == 0) {
    return(tibble::tibble(
      phenotype = character(), gene_id_clean = character(),
      annotation_type = character(), annotation_domain = character(),
      annotation_name = character(), annotation_label = character(),
      annotation_source = character()
    ))
  }
  genes <- dplyr::distinct(candidate_gene_summary, .data$phenotype, .data$gene_id_clean)
  if (anyNA(genes$phenotype) || anyNA(genes$gene_id_clean)) stop("Candidate gene keys contain missing values.")
  annotations <- read_local_functional_annotation(functional_annotation_file)
  candidate_records <- dplyr::inner_join(
    genes, annotations, by = "gene_id_clean", relationship = "many-to-many"
  ) %>%
    dplyr::mutate(annotation_source = "Local functional annotation")
  
  # TF names/families come from your existing TF annotation step.
  tf <- candidate_gene_summary
  for (column in c("tf_gene_name", "tf_family")) {
    if (!column %in% names(tf)) tf[[column]] <- NA_character_
    tf[[column]] <- dplyr::na_if(trimws(as.character(tf[[column]])), "")
  }
  tf <- tf %>%
    dplyr::filter(!is.na(.data$tf_gene_name) | !is.na(.data$tf_family)) %>%
    dplyr::transmute(
      phenotype = .data$phenotype, gene_id_clean = .data$gene_id_clean,
      annotation_type = "Transcription factor", annotation_domain = NA_character_,
      annotation_name = .data$tf_gene_name, annotation_label = .data$tf_family,
      annotation_source = "Local TF annotation"
    )
  records <- dplyr::distinct(dplyr::bind_rows(candidate_records, tf))
  records <- dplyr::arrange(records, .data$phenotype, .data$gene_id_clean,
                            .data$annotation_type, .data$annotation_domain, .data$annotation_name)
  n_matched <- dplyr::n_distinct(candidate_records$gene_id_clean)
  message("Candidate gene IDs with local records: ", n_matched, " / ", dplyr::n_distinct(genes$gene_id_clean))
  if (nrow(genes) > 0 && n_matched == 0) warning("No local annotations matched. Check gene IDs and releases.")
  records
}

# show_candidate_annotation_table()
# Displays a searchable gene summary with expandable annotation tables
# candidate_gene_summary = One summary row per phenotype-gene combination
# candidate_annotations = Full records from prepare_candidate_annotations()
# Main rows show description, TF name/family, SNP support and record counts
# Expanded rows show Type, Domain, Name, Label and Source for all records
# Returns: reactable HTML widget
show_candidate_annotation_table <- function(candidate_gene_summary, candidate_annotations) {
  if (nrow(candidate_gene_summary) == 0) {
    return(htmltools::div("No candidate genes were mapped."))
  }

  required <- c("phenotype", "gene_id_clean")
  if (!all(required %in% names(candidate_gene_summary))){
    stop("Missing phenotype/gene columns")
  }
  
  if (anyDuplicated(candidate_gene_summary[required])){
    stop("Expected one summary row per phenotype-gene.")
  }
  
  counts <- candidate_annotations %>%
    dplyr::count(.data$phenotype, .data$gene_id_clean, name = "n_annotation_records")
  display <- candidate_gene_summary %>%
    dplyr::select(dplyr::any_of(c(
      "phenotype",
      "gene_id_clean",
      "gene_description",
      "tf_gene_name",
      "tf_family",
      "n_lead_snps",
      "max_method_support",
      "min_distance_to_gene_bp",
      "closest_lead_snp",
      "closest_lead_snp_position"
    ))) %>%
    dplyr::left_join(counts, by = required)
  display$n_annotation_records <- dplyr::coalesce(display$n_annotation_records, 0L)
  definitions <- list(
    phenotype = reactable::colDef(name = "Variable"),
    gene_id_clean = reactable::colDef(name = "Gene ID", minWidth = 200),
    gene_description = reactable::colDef(name = "Description", minWidth = 250),
    tf_gene_name = reactable::colDef(name = "TF name"),
    tf_family = reactable::colDef(name = "TF family"),
    n_lead_snps = reactable::colDef(name = "Lead SNPs"),
    max_method_support = reactable::colDef(name = "Max. GEA support"),
    min_distance_to_gene_bp = reactable::colDef(name = "Closest SNP distance (bp)"),
    n_annotation_records = reactable::colDef(name = "Annotation records"),
    closest_lead_snp = reactable::colDef(name = "Closest lead SNP", minWidth = 220),
    closest_lead_snp_position = reactable::colDef(name = "SNP position", minWidth = 150)
  )
  reactable::reactable(
    display, searchable = TRUE, filterable = TRUE, defaultPageSize = 10,
    columns = definitions[intersect(names(definitions), names(display))],
    defaultColDef = reactable::colDef(na = "—", style = list(whiteSpace = "normal")),
    details = function(index) {
      records <- candidate_annotations %>%
        dplyr::filter(.data$phenotype == display$phenotype[index],
                      .data$gene_id_clean == display$gene_id_clean[index]) %>%
        dplyr::select(.data$annotation_type, .data$annotation_domain,
                      .data$annotation_name, .data$annotation_label, .data$annotation_source)
      if (nrow(records) == 0) return(htmltools::div("No local functional or TF annotations available."))
      htmltools::div(style = "padding: 12px;",
                     reactable::reactable(
                       records, searchable = TRUE, filterable = TRUE, defaultPageSize = 10,
                       defaultColDef = reactable::colDef(na = "—", style = list(whiteSpace = "normal")),
                       columns = list(
                         annotation_type = reactable::colDef(name = "Type"),
                         annotation_domain = reactable::colDef(name = "Domain"),
                         annotation_name = reactable::colDef(name = "Name", minWidth = 220),
                         annotation_label = reactable::colDef(name = "Label", minWidth = 300),
                         annotation_source = reactable::colDef(name = "Source")
                       )
                     )
      )
    }
  )
}

# plot_annotation_coverage()
# Plots candidate-gene annotation coverage by phenotype and annotation type.
# candidate_gene_summary = Gene summary used to count all candidate genes.
# candidate_annotations = Full local and TF records.
# Counts each gene once per type; labels show annotated/all candidate genes.
# Types overlap, so bars are not additive. This is not an enrichment test.
# Returns: ggplot object, or NULL when there are no annotation types.
plot_annotation_coverage <- function(candidate_gene_summary, candidate_annotations) {
  counts <- candidate_annotations %>%
    dplyr::filter(!is.na(.data$annotation_type)) %>%
    dplyr::distinct(.data$phenotype, .data$gene_id_clean, .data$annotation_type) %>%
    dplyr::count(.data$phenotype, .data$annotation_type, name = "n_genes")
  totals <- candidate_gene_summary %>%
    dplyr::distinct(.data$phenotype, .data$gene_id_clean) %>%
    dplyr::count(.data$phenotype, name = "total_genes")
  types <- sort(unique(counts$annotation_type))
  
  if (!length(types)){
    return(NULL)
  }
  
  grid <- expand.grid(phenotype = totals$phenotype, annotation_type = types, stringsAsFactors = FALSE)
  counts <- dplyr::left_join(grid, counts, by = c("phenotype", "annotation_type")) %>%
    dplyr::left_join(totals, by = "phenotype") %>%
    dplyr::mutate(n_genes = dplyr::coalesce(.data$n_genes, 0L),
                  count_label = paste0(.data$n_genes, "/", .data$total_genes))
  
  ggplot(counts, aes(x = n_genes, y = annotation_type)) +
    geom_col(fill = "#0072B2") +
    geom_text(aes(label = count_label), hjust = -0.1, size = 3) +
    facet_wrap(~ phenotype) +
    scale_x_continuous(breaks = scales::breaks_width(1), expand = ggplot2::expansion(mult = c(0, 0.25))) +
    labs(title = "Candidate gene annotation coverage", x = "Distinct candidate genes", y = NULL,
                  caption = "Labels: annotated genes / all candidate genes per variable. Annotation types overlap.") +
    theme_minimal(base_size = 11)
}

# plot_annotation_terms()
# Plots frequencies of GO terms, Pfam domains or other annotation records.
# candidate_annotations = Full local and TF records.
# phenotype_name = Environmental variable to display.
# annotation_type = Exact type to display, e.g. Gene Ontology, Pfam or InterPro.
# top_n = Maximum terms shown, ordered by count; NULL displays all terms.
# Counts distinct genes per term. GO BP/MF/CC remain distinguished.
# Returns: ggplot object, or NULL when the requested annotations are absent.
plot_annotation_terms <- function(candidate_annotations, phenotype_name,
                                  annotation_type = "Gene Ontology", top_n = 15) {
  type_name <- annotation_type
  counts <- candidate_annotations %>%
    dplyr::filter(.data$phenotype == phenotype_name, .data$annotation_type == type_name) %>%
    dplyr::distinct(.data$gene_id_clean, .data$annotation_domain,
                    .data$annotation_name, .data$annotation_label) %>%
    dplyr::group_by(.data$annotation_domain, .data$annotation_name, .data$annotation_label) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(.data$gene_id_clean), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$n_genes), .data$annotation_name, .data$annotation_domain)
  
  if (!nrow(counts)){
    return(NULL)
  }
  if (!is.null(top_n)){
    counts <- dplyr::slice_head(counts, n = top_n)
  }
  
  counts <- counts %>%
    dplyr::mutate(
      domain = dplyr::coalesce(.data$annotation_domain, "Not specified"),
      term = paste0(dplyr::coalesce(.data$annotation_name, "Unnamed"),
                    ifelse(is.na(.data$annotation_label), "", paste0(" — ", .data$annotation_label)),
                    ifelse(is.na(.data$annotation_domain), "", paste0(" [", .data$annotation_domain, "]")))
    )
  counts$term <- factor(counts$term, levels = rev(unique(counts$term)))
  
  ggplot(counts, aes(x = n_genes, y = term, fill = domain)) +
    geom_col() +
    scale_y_discrete(labels = function(x) str_wrap(x, 65)) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +
    labs(title = paste(phenotype_name, type_name, sep = ": "),
                  subtitle = if (is.null(top_n)) "All annotation terms" else paste("Up to", top_n, "most frequent terms"),
                  x = "Distinct candidate genes", y = NULL, fill = "Domain",
                  caption = "Descriptive annotation frequencies, not an enrichment test") +
    theme_minimal(base_size = 11) +
    theme(legend.position = if (type_name == "Gene Ontology") "bottom" else "none")
}

# plot_tf_families()
# Plots candidate transcription factor family counts for each phenotype.
# candidate_annotations = Full annotation records including Transcription factor.
# Counts each gene once within a TF family and phenotype.
# Returns: ggplot object, or NULL when no TF families are available.
plot_tf_families <- function(candidate_annotations) {
  counts <- candidate_annotations %>%
    dplyr::filter(.data$annotation_type == "Transcription factor", !is.na(.data$annotation_label)) %>%
    dplyr::distinct(.data$phenotype, .data$gene_id_clean, .data$annotation_label) %>%
    dplyr::count(.data$phenotype, .data$annotation_label, name = "n_genes")
  if (!nrow(counts)){
    return(NULL)
  }
  
  ggplot(counts, aes(x = n_genes, y = annotation_label)) +
    geom_col(fill = "#009E73") + ggplot2::facet_wrap(~ phenotype) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +
    labs(title = "Candidate transcription factor families", x = "Distinct candidate genes", y = "TF family") +
    theme_minimal(base_size = 11)
}


# run_biointerpretation_workflow()
# Full biological interpretation workflow, SNP-to-gene mapping + annotations
# config = Pipeline configuration object
# primary_lead_snps = Final LD-pruned SNP table from the configured primary set
# overwrite_annotation = defaults to FALSE, reuses existing annotation file
# Returns: Annotated SNP-gene mapping data frame
# Output: Also saves full candidate annotation records and annotation plots
run_biointerpretation_workflow <- function(
    config,
    primary_lead_snps,
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
    primary_lead_snps = primary_lead_snps,
    gene_annotation_file = gene_annotation_file,
    output_dir = output_dir,
    flank_bp = flank_bp,
    overwrite_annotation = overwrite_annotation
  )
  
  # Optional annotation layers
  biointerp_df <- annotate_candidate_genes(
    mapped_candidate_genes = biointerp_df,
    tf_annotation_file = tf_annotation_file,
    use_biomart = use_biomart,
    functional_annotation_file = config$biological_interpretation$functional_annotation_file
  )
  
  # Optional gene-level summary
  candidate_gene_summary <- summarize_candidate_genes(biointerp_df)
  
  # Preserve every local annotation and TF record for display and plotting.
  candidate_annotations <- prepare_candidate_annotations(
    candidate_gene_summary = candidate_gene_summary,
    functional_annotation_file =
      config$biological_interpretation$functional_annotation_file
  )
  
  write_csv(candidate_annotations,
            file.path(output_dir, "candidate_gene_annotations_long.csv"))
  saveRDS(candidate_annotations,
          file.path(output_dir, "candidate_gene_annotations_long.rds"))
  
  # Prepare plots here; the report only loads and prints them.
  annotation_plots <- list(
    coverage = plot_annotation_coverage(candidate_gene_summary, candidate_annotations),
    tf_families = plot_tf_families(candidate_annotations),
    by_variable = list()
  )
  
  for (variable in unique(candidate_gene_summary$phenotype)) {
    annotation_plots$by_variable[[variable]] <- list()
    for (type in c("Gene Ontology", "Pfam", "InterPro")) {
      annotation_plots$by_variable[[variable]][type] <- list(
        plot_annotation_terms(candidate_annotations, variable, type, top_n = 15)
      )
    }
  }
  
  saveRDS(annotation_plots, file.path(output_dir, "annotation_plots.rds"))
  
  # Save final outputs
  write_csv(biointerp_df, file.path(output_dir, "mapped_candidate_genes_annotated.csv"))
  write_csv(candidate_gene_summary, file.path(output_dir, "candidate_gene_summary.csv"))
  
  saveRDS(biointerp_df, file.path(output_dir, "mapped_candidate_genes_annotated.rds"))
  saveRDS(candidate_gene_summary, file.path(output_dir, "candidate_gene_summary.rds"))
  
  message("Mapped SNP-gene rows: ", nrow(biointerp_df))
  message("Candidate gene summary rows: ", nrow(candidate_gene_summary))
  
  return(biointerp_df)
}


