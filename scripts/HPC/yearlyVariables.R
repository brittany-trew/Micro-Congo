#' Calculate climate diversity in the vertical column
#' Climate variables (for each  height)
#' Mean daily tmin
#' q90 tmax
#' Mean diurnal range
#' SD of mean daily temp

#' ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path,"imbalanga/microParameters.R")) # loads worker functions
#' ----------------------------------------------------------------------------

array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
d <- array_id + offset

#' ----------------------------------------------------------------------------
#' Helper functions
q90_fun <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(quantile(x, probs = 0.90, na.rm = TRUE))
}

#' ----------------------------------------------------------------------------

inpath  <- paste0(head.path, "/dailyTemps/")
outpath <- paste0(head.path, "/yearlyTemps/")
dir.create(outpath, showWarnings = F)

yr.seq <- seq(2004, 2024, 1)
dates <- unlist(lapply(yr.seq, function(yr) {
  seq(
    as.Date(paste0(yr, "-01-01")),
    as.Date(paste0(yr, "-12-31")),
    by = "day"
  )
}))
dates <- as.Date(dates, origin = "1970-01-01")
idx <- format(dates, "%Y")

height_folders <- list.dirs(
  inpath,
  full.names = FALSE,
  recursive = FALSE
)

# Keep only folders named as a numeric height followed by "m"
height_folders <- height_folders[
  grepl("^[0-9]+(\\.[0-9]+)?m$", height_folders)
]

model.heights <- sort(as.numeric(sub("m$", "", height_folders)))

h <- sprintf("%.1f", model.heights[d])
print(h)

#' Load base files.
files_min <- list.files(
  paste0(inpath,h,"m/"),
  pattern    = paste0("DailyMin_.*tif$"),
  full.names = TRUE
)
min.stk  <- rast(files_min)
mon.tmin <- tapp(
  min.stk,
  index = idx,
  fun = mean,
  na.rm = TRUE
)
writeRaster(mon.tmin, paste0(outpath,"yearlyMean_tmin_",h,".tif"))

files_max <- list.files(
  paste0(inpath,h,"m/"),
  pattern    = paste0("DailyMax_.*tif$"),
  full.names = TRUE
)
max.stk  <- rast(files_max)
mon.tmax <- tapp(
  max.stk,
  index = idx,
  fun = mean,
  na.rm = TRUE
)
writeRaster(mon.tmax, paste0(outpath,"yearlyMean_tmax_",h,".tif"))
gc()
rm(min.stk)

#' Diurnal Range (Monthly mean)
mon.range <- mon.tmax - mon.tmin
writeRaster(mon.range, paste0(outpath,"yearlyDiurnalRange_",h,".tif"))
gc()

#' Q90 of Tmax (Monthly)
mon.q90 <- tapp(
  max.stk,
  index = idx,
  fun = q90_fun
)
writeRaster(mon.q90, paste0(outpath,"yearlyq90_tmax_",h,".tif"))
rm(max.stk)
gc()

#' Monthly SD of mean.
files_mean <- list.files(
  paste0(inpath,h,"m/"),
  pattern    = paste0("DailyMean_.*tif$"),
  full.names = TRUE
)
mean.stk  <- rast(files_mean)
mon.sdmean <- tapp(
  mean.stk,
  index = idx,
  fun = sd,
  na.rm = TRUE
)
writeRaster(mon.sdmean, paste0(outpath,"yearlySD_tmean_",h,".tif"))
gc()
rm(mean.stk)
