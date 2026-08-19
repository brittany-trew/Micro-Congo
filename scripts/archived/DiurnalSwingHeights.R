#' Diurnal swing by height
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

for(h in 1:length(model.heights)){
  height <- sprintf("%04.1fm", model.heights[h])
  
  file.path <- paste0(inpath,"dailyTemps/",height,"/")
  files_tmax <- list.files(file.path, pattern = paste0("DailyMax_.*.tif$"), 
                           recursive = T, full.names = TRUE)
  tmax_full <- rast(files_tmax)
  
  r <- resample(tmax_full,chm)
  r_masked <- ifel(chm < h-1, NA, r)
  
  r_max <- max(r_masked, na.rm = TRUE)
  
  n_days <- nlyr(rast(r.lst[[1]]))  # 366
  index <- rep(1:n_days, times = length(r.lst))
  r_max <- tapp(r.stk.max, index, fun = max, na.rm = TRUE)
  rm(r.stk.max)
  plot(r_max[[1]])
  
}


for(i in 1:length(yr.seq)){
  yr <- yr.seq[i]
  
  files_tmax <- list.files(inpath, pattern = paste0("DailyMax_",yr,".tif$"), 
                           recursive = TRUE, full.names = TRUE)
  
  tmax_full <- rast(files_tmax)
  #' remove any values above canopy
  r.lst <- list()
  for(b in 1:length(files_tmax)){
    f <- files_tmax[[b]]
    hh <- as.numeric(unlist(strsplit(unlist(strsplit(f,"/"))[[9]],"m")))
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
