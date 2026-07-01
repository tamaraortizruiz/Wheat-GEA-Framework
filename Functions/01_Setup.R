# ---- Setup ----

# create_project_dirs()
# Creates the pipeline's project folder structure
# overwrite = If TRUE, removes and remakes pipeline folders
# Output: Creates required project directories
create_project_dirs <- function(overwrite = FALSE) {
  dirs <- c(
    "Data",
    "Output",
    "Output/Structure",
    "Output/GEA",
    "Output/GEA/GEMMA",
    "Output/GEA/LFMM",
    "Output/GEA/RDA",
    "Output/GEA/pcadapt",
    "Output/ConsensusSNP",
    "Output/ConsensusLD",
    "Output/ConsensusRobustness"
  )
  
  # Remove existing directories if overwrite = TRUE
  if (overwrite) {
    message("overwrite = TRUE")
    message("Removing existing Data and Output directories")
    for (d in rev(dirs)) {
      if (dir.exists(d)) {
        unlink(d, recursive = TRUE, force = TRUE)
      }
    }
  }
  
  # Create directories
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}