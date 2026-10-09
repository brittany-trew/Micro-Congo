#' Calc. Vertical Heterogeneity.
scripts.path <- "scripts/"
source(paste0(scripts.path,"imbalanga/microParameters.R")) # loads worker functions

# -------------------------------
# Cluster batch command
# -------------------------------
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
ft <- array_id + offset

#' Helpers:
collapse_vertical <- function(files, chm) {
  r.lst <- lapply(files, load_masked, chm = chm)
  r.stk <- rast(r.lst)
  r.stk[r.stk < 0] <- NA
  
  n_days <- nlyr(r.lst[[1]])
  index  <- rep(seq_len(n_days), times = length(r.lst))
  
  colmax <- tapp(
    r.stk, index,
    fun = max,
    na.rm = TRUE
  )
  
  colmin <- tapp(
    r.stk, index,
    fun = min,
    na.rm = TRUE
  )
  
  colsd <- tapp(
    r.stk, index,
    fun = sd,
    na.rm = TRUE
  )
  
  col95 <- tapp(
    r.stk,
    index,
    fun = function(x, ...) {
      quantile(
        x,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
      )
    }
  )
  
  gc()
  
  return(list(
    colmax = colmax,
    colmin = colmin,
    colsd  = colsd,
    col95  = col95
  ))
}
# Extract height (metres) from parent folder name
# e.g. ".../dailyTemps/5.0m/DailyMax_2009.tif" -> 5
get_height <- function(f) {
  as.numeric(gsub("m$", "", basename(dirname(f))))
}
# Load a single height-band raster, resample to CHM grid, then mask pixels
# where canopy height is below this height band (i.e. above canopy)
load_masked <- function(f, chm) {
  hh <- get_height(f)
  r  <- resample(rast(f), chm)
  ifel(chm <= hh, NA, r)
}

#' Pathways.
data.path <- paste0(head.path,"dailyTemps/")
strata.path <- paste0(head.path,"/strata/")
outpath <- paste0(head.path,"/strataSummary/")
dir.create(outpath, showWarnings = F)

strata.file <- paste0(strata.path, "StrataZones_ft_", ft, ".RDS")
strata <- readRDS(strata.file)
strata <- strata[order(strata$zone_start, strata$zone_end), , drop = FALSE]
strata$zone <- seq_len(nrow(strata))

chm <- rast(paste0(head.path, "chm.tif"))

# ── Process each year and stratum ---------------------------------------------
yr.seq <- 2004:2024

for (yr in yr.seq) {
  message("Processing forest type ", ft, ", year ", yr, " ...")
  
  files_max <- list.files(
    data.path,
    pattern    = paste0("DailyMax_", yr, "\\.tif$"),
    recursive  = TRUE,
    full.names = TRUE
  )
  
  heights <- vapply(files_max, get_height, numeric(1))
  file.order <- order(heights)
  files_max  <- files_max[file.order]
  heights    <- heights[file.order]
  
  for (z in seq_len(nrow(strata))) {
    zone.start <- strata$zone_start[z]
    zone.end   <- strata$zone_end[z]

    suffix <- paste0(
      "_f", ft,
      "_s", z,
      "_", yr, ".tif"
    )
    
    zone_files <- files_max[
      !is.na(heights) &
        heights >= zone.start &
        heights < zone.end
    ]
    
    # Max, Min and SD of vertical columns.
    minmax <- collapse_vertical(zone_files, chm)
    gc()
    # Daily vertical spread in peak temperatures
    r_range <- minmax$colmax - minmax$colmin
    
    writeRaster(minmax$colmax, paste0(outpath, "colmax_", suffix), overwrite = TRUE)
    writeRaster(r_range,  paste0(outpath, "verticalRange_", suffix), overwrite = TRUE)
    writeRaster(minmax$colsd,  paste0(outpath, "verticalSD_", suffix), overwrite = TRUE)
    writeRaster(minmax$col95, paste0(outpath, "colmax95_", suffix), overwrite = TRUE)
    
    message("Output written for ", yr)
    rm(minmax, r_range)
    gc()
  }
  
  message("Outputs written for forest type ", ft, ", year ", yr)
}