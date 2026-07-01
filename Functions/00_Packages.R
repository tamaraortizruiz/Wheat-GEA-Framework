# ---- Package installation and loading ----

cran_packages <- c(
  "bigsnpr",
  "bigstatsr",
  "dplyr",
  "readr",
  "yaml",
  "ggplot2",
  "adegenet",
  "statgenGWAS",
  "ade4",
  "terra",
  "geodata",
  "rnaturalearth",
  "rnaturalearthdata",
  "reactable",
  "vegan",
  "pcadapt",
  "igraph",
  "htmltools"
)

bio_packages <- c(
  "LEA"
)

# Install missing CRAN packages
missing_cran <- cran_packages[!cran_packages %in% rownames(installed.packages())]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

# Install BiocManager if needed
if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager")
}

# Install missing Bioconductor packages
missing_bio <- bio_packages[!bio_packages %in% rownames(installed.packages())]
if (length(missing_bio) > 0) {
  BiocManager::install(missing_bio, ask = FALSE, update = FALSE)
}

# Load all packages
all_packages <- c(cran_packages, bio_packages)
invisible(
  lapply(all_packages, function(pkg) {
    library(pkg, character.only = TRUE)
  })
)