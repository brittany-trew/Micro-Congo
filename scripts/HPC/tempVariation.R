#' Calc. Vertical Heterogeneity.
scripts.path <- "scripts/"
source(paste0(scripts.path,"imbalanga/microParameters.R")) # loads worker functions

# -------------------------------
# Cluster batch command
# -------------------------------
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
d <- array_id + offset

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


chm <- rast(paste0(head.path,"chm.tif"))

inpath <- paste0(head.path,"dailyTemps/")
outpath <- paste0(head.path,"verticalSummary/")
dir.create(outpath, showWarnings = F)

#' Loop for each year.
#' Loads daily max for each year. 
#' Calculates summaries of the vertical columns for each day. 
#' and saves as yearly files: eg., verticalSD_2018.tif.
yr.seq <- seq(2004, 2024, 1)

for (i in 1:length(yr.seq)) {
  yr <- yr.seq[[i]]
  message("Processing ", yr, " ...")
  
  files_max <- list.files(
    inpath,
    pattern    = paste0("DailyMax_", yr, ".tif$"),
    recursive = T,
    full.names = TRUE
  )
  
  heights <- as.numeric(sub(".*/dailyTemps/([0-9.]+)m/.*", "\\1", files_max))
  files_max <- files_max[order(heights)]
  
  # Max, Min and SD of vertical columns.
  minmax <- collapse_vertical(files_max, chm)
  gc()
  # Daily vertical spread in peak temperatures
  r_range <- minmax$colmax - minmax$colmin
  
  writeRaster(minmax$colmax, paste0(outpath, "colmax_", yr, ".tif"), overwrite = TRUE)
  writeRaster(r_range,  paste0(outpath, "verticalRange_",  yr, ".tif"), overwrite = TRUE)
  writeRaster(minmax$colsd,  paste0(outpath, "verticalSD_",  yr, ".tif"), overwrite = TRUE)
  writeRaster(minmax$col95, paste0(outpath, "colmax95_", yr, ".tif"), overwrite = TRUE)
  
  message("Output written for ", yr)
  rm(minmax, r_range)
  gc()
}