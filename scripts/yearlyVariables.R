#' Calculate climate diversity in the vertical column
#' (A) Whole column
#' (B) By zone in different forest types
#' (C) Does climate diversity correlate with warming in tmax over time?

#' Climate variables (for each  height)
#' Mean daily tmin
#' q90 tmax
#' Mean diurnal range
#' SD of mean daily temp

# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

terraOptions(
  tempdir = "/Volumes/MicroMaze2/tmp",
  memfrac = 0.7
)

sample.name <- all.samples[[4]]
head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")

inpath  <- paste0(volumes3, sample.name, "/")
chm           <- rast(paste0(head.path, "/chm.tif"))
pai.layers    <- list.files(paste0(head.path, "/pai"), 
                            pattern = "_to_canopy.tif$")
n.h           <- length(pai.layers)
dzd = 1
model.heights <- c(0.2,seq(dzd,(n.h-1)*dzd,by = dzd))
yr.seq        <- seq(2004, 2024, 1)
if(sample.name == "MobaCluster"){yr.seq <- yr.seq[yr.seq != 2007]}
dates <- unlist(lapply(yr.seq, function(yr) {
  seq(
    as.Date(paste0(yr, "-01-01")),
    as.Date(paste0(yr, "-12-31")),
    by = "day"
  )
}))

dates <- as.Date(dates, origin = "1970-01-01")
idx <- format(dates, "%Y")

outpath2 <- paste0(inpath, "monthlyTemps/")
dir.create(outpath2, showWarnings = F)

q90_fun <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(quantile(x, probs = 0.90, na.rm = TRUE))
}

#' Variable creation
for(i in 1:length(model.heights)){
  h <- sprintf("%.1f", model.heights[i])
  print(h)
  
  #' Load base files.
  files_min <- list.files(
    paste0(inpath,"dailyTemps/",h,"m/"),
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
  writeRaster(mon.tmin, paste0(outpath2,"yearlyMean_tmin_",h,".tif"))
  
  files_max <- list.files(
    paste0(inpath,"dailyTemps/",h,"m/"),
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
  writeRaster(mon.tmax, paste0(outpath2,"yearlyMean_tmax_",h,".tif"))
  gc()
  rm(min.stk)
  
  #' Diurnal Range (Monthly mean)
  mon.range <- mon.tmax - mon.tmin
  writeRaster(mon.range, paste0(outpath2,"yearlyDiurnalRange_",h,".tif"))
  gc()
  
  #' Q90 of Tmax (Monthly)
  mon.q90 <- tapp(
    max.stk,
    index = idx,
    fun = q90_fun
  )
  writeRaster(mon.q90, paste0(outpath2,"yearlyq90_tmax_",h,".tif"))
  rm(max.stk)
  gc()
  
  #' Monthly SD of mean.
  files_mean <- list.files(
    paste0(inpath,"dailyTemps/",h,"m/"),
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
  writeRaster(mon.sdmean, paste0(outpath2,"yearlySD_tmean_",h,".tif"))
  gc()
  rm(mean.stk)
}
