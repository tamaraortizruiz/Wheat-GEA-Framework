# ---- Environmental variables ----
# extract_climate_variables()
# Extracts WorldClim climate variables for each sample using its longitude and latitude
# metadata = Metadata data frame containing coordinates
# lon_col = Column name for longitude
# lat_col = Column name for latitude
# sample_col = Column name for sample ID
# var = WorldClim variable type, defaults to "bio" for bioclimatic variables
# res = Spatial resolution of WorldClim data
# climate_dir = Folder where climate raster files are stored
# output_file = Path to save extracted climate variables
# overwrite = If FALSE, reuse existing climate output file
# Output: .csv file containing sample IDs and extracted climate variables
# Returns: Data frame with one row per sample and climate variables as columns
extract_climate_variables <- function(
    metadata,
    lon_col = "Longitude",
    lat_col = "Latitude",
    sample_col = "SeedID",
    var = "bio",
    res = 5,
    climate_dir = "Data/Climate",
    output_file = "Data/climate_variables.csv",
    overwrite = FALSE
) {
  
  # if file exists and overwrite is FALSE, load existing files
  if (file.exists(output_file) && !overwrite) {
    message("Loading existing climate data: ", output_file)
    return(read.csv(output_file, check.names = FALSE))
  }
  
  needed_cols <- c(lon_col, lat_col, sample_col)
  missing_cols <- needed_cols[!needed_cols %in% colnames(metadata)]
  
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }
  
  # Samples with no longitude/latitude information
  n_before <- nrow(metadata)
  # Filter samples with no longitude/latitude information
  metadata <- metadata %>%
    filter(!is.na(.data[[lon_col]]),
           !is.na(.data[[lat_col]]),
           !( .data[[lon_col]] == 0 & .data[[lat_col]] == 0 ))
  n_after <- nrow(metadata)
  
  # Create directory for extracted climate data
  dir.create(climate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Extract raster data
  clim_rasters <- worldclim_global(
    var = var,
    res = res,
    path = climate_dir,
    version = "2.1"
  )
  
  # Bioclimatic variables
  if (var == "bio") {
    names(clim_rasters) <- paste0("BIO", 1:19)
  }
  
  # Create spatial points
  points <- vect(
    metadata,
    geom = c(lon_col, lat_col),
    crs = "EPSG:4326"
  )
  
  # Extract climatic variables
  climate_df <- extract(clim_rasters, points)
  
  message("Climate data extracted for ", nrow(climate_df), " samples.")
  
  climate_df <- climate_df[, -1, drop = FALSE]
  climate_df[[sample_col]] <- metadata[[sample_col]]
  climate_df <- climate_df[, c(sample_col, setdiff(colnames(climate_df), sample_col))]
  
  # Keep only complete cases (no NAs)
  climate_df <- climate_df[complete.cases(climate_df),]
  
  write.csv(climate_df, output_file, row.names = FALSE)

  message("\nClimate extraction summary\n",
          "Samples in original dataset: ", n_before, "\n",
          "Samples with coordinate data: ", n_after, "\n",
          "BIO variables extracted: ", ncol(climate_df) - 1, "\n",
          "Samples with climate data: ", nrow(climate_df))
  message("Climate variables saved to: ", output_file)
  
  return(climate_df)
}

# plot_climate_map()
# Creates a map of geographic locations of accessions with extracted climate data
# Points are colored by selected climatic variable
# climate_data = Data frame of extracted climate variables
# metadata = Metadata data frame with longitude and latitude columns
# sample_col = Column used to match climate_data and metadata
# lon_col = Metadata column containing longitude
# lat_col = Metadata column containing latitude
# color_col = Climate variable used to color map points
# Output: A ggplot map with accession locations over a world map
# Returns: ggplot object.
plot_climate_map <- function(
    climate_data,
    metadata,
    sample_col = "SeedID",
    lon_col = "Longitude",
    lat_col = "Latitude",
    color_col = "BIO1"
) {
  
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  plot_data <- climate_data %>%
    left_join(metadata[, c(sample_col, lon_col, lat_col)], by = sample_col)
  
  p <- ggplot() +
    geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
    geom_point(data = plot_data,
               aes(x = .data[[lon_col]],
                   y = .data[[lat_col]],
                   color = .data[[color_col]]),
               alpha = 0.8,
               size = 2) +
    coord_sf() +
    labs(title = paste("Geographic distribution colored by", color_col),
         x = "Longitude",
         y = "Latitude",
         color = color_col) +
    theme_classic()
  
  return(p)
}

# plot_climate_cor_heatmap()
# Creates a correlation matrix and heatmap among bioclimatic variables to identify highly correlated environmental predictors
# climate_data = Data frame containing extracted climate variables.
# bio_prefix = Prefix to identify bioclimatic variables, defaults to "BIO".
# Output: A heatmap showing pairwise Pearson correlations among BIO variables
# Returns: ggplot object
plot_climate_cor_heatmap <- function(
    climate_data,
    bio_prefix = "BIO"
) {
  
  # BIO variable names
  bio_cols <- grep(paste0("^", bio_prefix), colnames(climate_data), value = TRUE)
  
  # Generate correlation matrix
  cor_mat <- cor(climate_data[, bio_cols], use = "pairwise.complete.obs")
  
  # As data frame
  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Var1", "Var2", "Correlation")
  
  # Plot
  p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = 0,
                         limits = c(-1, 1)) +
    coord_fixed() +
    labs(title = "Correlation among bioclimatic variables",
         x = NULL,
         y = NULL,
         fill = "r") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8))
  
  return(p)
}

# run_climate_pca()
# Performs PCA on extracted bioclimatic variables to summarize main climatic gradients
# climate_data = Data frame of extracted climate variables
# bio_prefix = Prefix used to identify bioclimatic variables, defaults to "BIO"
# output_file = Optional path to save the PCA result as .rds file
# Output: PCA object, PCA scores, variance explained, and list of variables used
# Returns: List containing pca, scores, variance and variables
run_climate_pca <- function(
    climate_data,
    bio_prefix = "BIO",
    output_file = NULL
) {
  
  bio_cols <- grep(paste0("^", bio_prefix), colnames(climate_data), value = TRUE)
  
  env_data <- climate_data[, bio_cols]
  env_data <- env_data[complete.cases(env_data), ]
  
  # PCA
  env_pca <- prcomp(env_data, center = TRUE, scale. = TRUE)
  
  # Variance explained
  variance <- data.frame(
    PC = paste0("EnvPC", seq_along(env_pca$sdev)),
    Variance = (env_pca$sdev^2 / sum(env_pca$sdev^2)) * 100,
    CumulativeVariance = cumsum((env_pca$sdev^2 / sum(env_pca$sdev^2)) * 100)
  )
  
  # PC scores
  scores <- as.data.frame(env_pca$x)
  colnames(scores) <- paste0("EnvPC", seq_len(ncol(scores)))
  
  # Variable loadings
  loadings <- as.data.frame(env_pca$rotation)
  colnames(loadings) <- paste0("EnvPC", seq_len(ncol(loadings)))
  loadings$Variable <- rownames(loadings)
  loadings <- loadings[, c("Variable", setdiff(colnames(loadings), "Variable"))]
  
  result <- list(
    pca = env_pca,
    scores = scores,
    loadings = loadings,
    variance = variance,
    variables = bio_cols
  )
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, output_file)
    message("Climate PCA saved to: ", output_file)
  }
  
  return(result)
}

# plot_climate_pca()
# Plots selected environmental principal components from climate PCA result
# climate_pca_result = Output from run_climate_pca()
# pc_x = Environmental PC for x-axis
# pc_y = Environmental PC for y-axis
# Output: PCA scatterplot of bioclimatic variation among samples
# Returns: ggplot object
plot_climate_pca <- function(
    climate_pca_result,
    pc_x = 1,
    pc_y = 2
) {

  scores <- climate_pca_result$scores
  variance <- climate_pca_result$variance
  
  # PCs to plot
  x_col <- paste0("EnvPC", pc_x)
  y_col <- paste0("EnvPC", pc_y)
  
  x_lab <- paste0(x_col, " (", round(variance$Variance[pc_x], 2), "%)")
  y_lab <- paste0(y_col, " (", round(variance$Variance[pc_y], 2), "%)")
  
  p <- ggplot(scores, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.8) +
    labs(title = "PCA of bioclimatic variables",
         x = x_lab,
         y = y_lab) +
    theme_classic()
  
  return(p)
}