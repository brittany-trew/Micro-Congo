scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R"))

inpath <- paste0(volumes3,sample.name,"/")
outpath <- paste0(out.data,sample.name,"/")
dir.create(outpath)
outpath <- paste0(outpath,"Vertical_Summaries/")
dir.create(outpath)

pai.layers <- list.files(paste0(head.path,"/pai"))
n.h <- length(pai.layers)
model.heights <- seq(0,(n.h-1)*dzd,by = dzd)

chm <- rast(paste0(head.path,"/chm.tif"))
plot(chm)

yr.seq <- seq(2004, 2024, 1)

for(i in 1:length(yr.seq)){
  yr <- yr.seq[i]
  
  files_tmax <- list.files(inpath, pattern = paste0("DailyMax_",yr,".tif$"), 
                           recursive = TRUE, full.names = TRUE)
  
  tmax_full <- rast(files_tmax)
  #' remove any values above canopy
  r.lst <- list()
  for(b in 1:length(files_tmax)){
    f <- files_tmax[[b]]
    hh <- as.numeric(unlist(strsplit(unlist(strsplit(f,"/"))[[7]],"m")))
    r <- rast(f)
    r <- resample(r,chm)
    r_masked <- ifel(chm < hh, NA, r)
    r.lst[[b]] <- r_masked
  }
  r.stk.max <- rast(r.lst)
  n_days <- nlyr(rast(r.lst[[1]]))  # 366
  index <- rep(1:n_days, times = length(r.lst))
  r_max <- tapp(r.stk.max, index, fun = max, na.rm = TRUE)
  rm(r.stk.max)
  plot(r_max[[1]])
  writeRaster(r_max, paste0(outpath,"tmax_",yr,".tif"))
  rm(r_max)
  gc()
}


for(i in 1:length(yr.seq)){
  yr <- yr.seq[i]
  
  files_tmin <- list.files(inpath, pattern = paste0("DailyMin_",yr,".tif$"), 
                           recursive = TRUE, full.names = TRUE)
  
  tmin_full <- rast(files_tmin)
  #' remove any values above canopy
  r.lst <- list()
  for(b in 1:length(files_tmin)){
    f <- files_tmin[[b]]
    hh <- as.numeric(unlist(strsplit(unlist(strsplit(f,"/"))[[7]],"m")))
    r <- rast(f)
    r <- resample(r,chm)
    r_masked <- ifel(chm < hh, NA, r)
    r.lst[[b]] <- r_masked
  }
  r.stk.min <- rast(r.lst)
  n_days <- nlyr(rast(r.lst[[1]]))  # 366
  index <- rep(1:n_days, times = length(r.lst))
  r_min <- tapp(r.stk.min, index, fun = min, na.rm = TRUE)
  rm(r.stk.min)
  plot(r_min[[1]])
  writeRaster(r_min, paste0(outpath,"tmin_",yr,".tif"))
  rm(r_min)
  gc()
}

#' Compute Range (Tmax - Tmin) of Vertical Profiles
for(i in 1:length(yr.seq)){
  
  yr <- yr.seq[i]
  rmax <- rast(paste0(outpath,"tmax_",yr,".tif"))
  rmin <- rast(paste0(outpath,"tmin_",yr,".tif"))
  
  er <- rmin > rmax
  rmin[er] <- NA
  rmax[er] <- NA
  
  r.range <- rmax-rmin
  plot(r.range[[1]])
  writeRaster(r.range, paste0(outpath,"range_",yr,".tif"), overwrite = T)
}

all.files <- list.files(outpath, pattern = paste0("range_.*.tif$"), 
                        full.names = TRUE)
range.allyrs <- rast(all.files)

#' Flatten Vertical Profiles to Summaries
mean_range <- app(range.allyrs, mean, na.rm = TRUE)
sd_range   <- app(range.allyrs, sd,   na.rm = TRUE)
cv_range <- sd_range / mean_range

writeRaster(mean_range, paste0(outpath,"Mean_range_2004to2024.tif"), overwrite= T)
writeRaster(cv_range, paste0(outpath,"Cvariance_range_2004to2024.tif"), overwrite= T)
