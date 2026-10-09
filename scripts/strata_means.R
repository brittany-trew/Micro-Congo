# Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

#' Take the mean overall.
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tempdir())

for(i in 1:length(all.samples)){
  
  sample.name <- all.samples[[i]]
  rpath <- paste0(out.data,sample.name,"/rasters")
  dir.create(rpath)
  
  pca.path <- paste0(out.data,"pca_clusters/")
  pca <- rast(paste0(pca.path, sample.name, "_pca_C4.tif"))
  
  for(f in 1:4){
    ft <- f
    
    spath <- paste0(out.data, "strata/")
    strata <- readRDS(paste0(spath,"StrataZones_ft_",ft,".RDS"))
    
    inpath <- paste0(out.data,sample.name,"/strataSummary/")
    
    for(z in 1:nrow(strata)){
      
      z1 <- strata[z,]
      
      if(!file.exists(paste0(rpath, "/VerticalSDinDailyTmax_2004to2024_f",f,"_",z1$zone_start,"to",z1$zone_end,".tif"))
      ){
        sd_all  <- rast(list.files(inpath, pattern = paste0("^verticalSD__f",f,"_s",z,".*.tif$"),  full.names = TRUE))
        mean_sd <- app(sd_all, mean, na.rm = TRUE)
        
        sd_q98 <- terra::global(
          mean_sd,
          fun = stats::quantile,
          probs = 0.75,
          na.rm = TRUE
        )[1, 1]
        mean_sd_fixed <- ifel(mean_sd > sd_q98, NA, mean_sd)
        
        mean_sd_fixed[pca != ft] <- NA
        plot(mean_sd_fixed)
        
        writeRaster(
          mean_sd_fixed,
          paste0(paste0(rpath, "/VerticalSDinDailyTmax_2004to2024_f",f,"_",z1$zone_start,"to",z1$zone_end,".tif")),
          overwrite = TRUE
        )
      }
      
    }
  }
  
  file.remove(list.files(rpath, pattern = paste0("^sdVertical.*",z1$zone_start,"to",z1$zone_end,".tif$"),  full.names = TRUE))
  
  
  #' Temporal means.
  for(f in 1:4){
    ft <- f
    
    spath <- paste0(out.data, "strata/")
    strata <- readRDS(paste0(spath,"StrataZones_ft_",ft,".RDS"))
    
    inpath <- paste0(out.data,sample.name,"/verticalSummary/")
    
    for(z in 1:nrow(strata)){
      
      z1 <- strata[z,]
      
      if(!file.exists(paste0(rpath, "/TemporalSDinDailyTmax_2004to2024_f",f,"_",z1$zone_start,"to",z1$zone_end,".tif"))
      ){
        
        files <- list.files(
          inpath,
          pattern = "^temporalSD_[0-9.]+m\\.tif$",
          full.names = TRUE
        )
        
        heights <- as.numeric(
          sub("^temporalSD_([0-9.]+)m\\.tif$", "\\1", basename(files))
        )
        
        keep <- heights >= z1$zone_start & heights < z1$zone_end
        
        zone_files <- files[keep]
        zone_heights <- heights[keep]
        
        # Sort numerically by height
        zone_files <- zone_files[order(zone_heights)]
        
        r <- terra::rast(zone_files)
        
        mean_sd <- app(r, mean, na.rm = TRUE)
        
        writeRaster(
          mean_sd,
          paste0(paste0(rpath, "/TemporalSDinDailyTmax_2004to2024_f",f,"_",z1$zone_start,"to",z1$zone_end,".tif")),
          overwrite = TRUE
        )
      }
      
    }
  }
  
}



