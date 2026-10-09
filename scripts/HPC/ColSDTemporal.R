#' Calc. Temporal Heterogeneity by Stratum.
scripts.path <- "scripts/"
source(paste0(scripts.path, "imbalanga/microParameters.R"))

# -------------------------------
# Cluster batch command
# -------------------------------
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
ft <- array_id + offset

# Extract height from parent folder name
get_height <- function(f) {
  as.numeric(gsub("m$", "", basename(dirname(f))))
}

# Load raster and mask heights at or above the canopy
load_masked <- function(f, chm) {
  hh <- get_height(f)
  r <- resample(rast(f), chm)
  ifel(chm <= hh, NA, r)
}

# Calculate temporal SD across all days and years at one height
collapse_temporal <- function(files, chm) {
  r.lst <- lapply(files, load_masked, chm = chm)
  r.stk <- rast(r.lst)
  r.stk[r.stk < 0] <- NA
  
  app(r.stk, sd, na.rm = TRUE)
}

# Pathways
data.path <- file.path(head.path, "dailyTemps")
strata.path <- file.path(head.path, "strata")
outpath <- file.path(head.path, "strataSummary")
dir.create(outpath, showWarnings = FALSE)

strata.file <- file.path(
  strata.path,
  paste0("StrataZones_ft_", ft, ".RDS")
)
strata <- readRDS(strata.file)
strata <- strata[
  order(strata$zone_start, strata$zone_end),
  ,
  drop = FALSE
]
strata$zone <- seq_len(nrow(strata))

chm <- rast(file.path(head.path, "chm.tif"))

# Find daily maximum temperature files for all years
yr.seq <- 2004:2024

files_max <- list.files(
  data.path,
  pattern = paste0(
    "^DailyMax_(",
    paste(yr.seq, collapse = "|"),
    ")\\.tif$"
  ),
  recursive = TRUE,
  full.names = TRUE
)

heights <- vapply(files_max, get_height, numeric(1))
file.order <- order(heights, files_max)
files_max <- files_max[file.order]
heights <- heights[file.order]

# Process each stratum
for (z in seq_len(nrow(strata))) {
  zone.start <- strata$zone_start[z]
  zone.end <- strata$zone_end[z]
  
  zone_heights <- sort(unique(heights[
    !is.na(heights) &
      heights >= zone.start &
      heights < zone.end
  ]))
  
  # Temporal SD at each height, using all years together
  temporal_sd <- lapply(zone_heights, function(hh) {
    height_files <- files_max[which(heights == hh)]
    collapse_temporal(height_files, chm)
  })
  
  # Mean temporal SD across heights within the stratum
  mean_sd <- app(
    rast(temporal_sd),
    mean,
    na.rm = TRUE
  )
  names(mean_sd) <- "mean_temporal_sd"
  
  writeRaster(
    mean_sd,
    file.path(
      outpath,
      paste0("meanTemporalSD_f", ft, "_s", z, ".tif")
    ),
    overwrite = TRUE
  )
  
  rm(temporal_sd, mean_sd)
  gc()
}