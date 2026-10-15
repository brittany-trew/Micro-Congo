#' Main run script.

scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))


#' Part 1: Process LiDAR data.
source(paste0(scripts.path, "ProcessLiDAR.R"))

#' Part 2: Derive drivers for microclimate modelling.
source(paste0(scripts.path, "land_cover/process_landcover_",sample.name,".R"))
source(paste0(scripts.path, "GroundMetrics.R"))
