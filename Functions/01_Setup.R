# ---- Setup ----

# create_project_dirs()
# Creates the pipeline's folder structure
# overwrite = If TRUE, removes and remakes existing folders
# Output: Creates directories if they do not already exist:
# Data, Output, Output/Kinship, Output/Structure, Output/GEA and Output/Plots
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
    "Output/Plots",
    "Logs"
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
  
  # Remake directories
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}