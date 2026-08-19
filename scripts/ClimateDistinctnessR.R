#' Distinctness Calcs. 
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

#' ----------------------------------------------------------------------------
#' Helper functions
clim_dist_fun <- function(x) {
  #' Convert the vector of climate values into a matrix where:
  #' rows = heights within the forest column
  #' columns = climate variables
  mat <- matrix(
    x,
    nrow = n.h,
    ncol = 4,
    byrow = FALSE)
  #' Remove heights containing missing values
  mat <- mat[complete.cases(mat), , drop = FALSE]
  #' Require at least two valid heights to calculate
  #' pairwise climatic distinctness
  if(nrow(mat) < 2) {
    return(c(mean_dist = NA_real_, sd_dist = NA_real_))}
  #' Calculate Euclidean pairwise climatic distance
  #' among all heights within the vertical column
  d <- as.numeric(dist(mat))
  #' mean_dist:
  #' Average climatic distinctness among heights
  #' sd_dist:
  #' Variability in climatic distinctness within the column
  c(mean_dist = mean(d),
    sd_dist   = sd(d))
}


load_scaled <- function(files, scale_pars) {
  r <- rast(files)
  names(r) <- sub("^X", "", names(r))
  r[!is.finite(r)] <- NA
  r <- (r - scale_pars$mean) / scale_pars$sd
  return(r)
}

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
inpath_all <- "/Volumes/MicroMaze2/"
all.scales <- readRDS(paste0(out.data,"climate_distinctness/scaling.RDS"))


#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' By site:
all.samples <- all.samples[-1]

for(a in 1:length(all.samples)){
  sample.name <- all.samples[[a]]
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  inpath  <- paste0(volumes3, sample.name, "/")
  
  chm           <- rast(paste0(head.path, "/chm.tif"))
  pai.layers    <- list.files(paste0(head.path, "/pai"), 
                              pattern = "_to_canopy.tif$")
  n.h           <- length(pai.layers)
  dzd = 1
  model.heights <- c(0.2,seq(dzd,(n.h-1)*dzd,by = dzd))
  yr.seq        <- seq(2004, 2024, 1)
  
  dates <- seq(
    as.Date(paste0(yr.seq[1],"-01-01")),
    as.Date(paste0(yr.seq[length(yr.seq)],"-12-31")),
    by = "day"
  )
  idx <- format(dates, "%Y")
  outpath2 <- paste0(inpath, "monthlyTemps/")
  
  tmin.mean <- list.files(outpath2, pattern = paste0("_tmin_"), full.names = T)
  tmin.mean <- load_scaled(tmin.mean, all.scales$tmin)
  
  tmean.sd <- list.files(outpath2, pattern = paste0("SD_tmean_"), full.names = T)
  tmean.sd <- load_scaled(tmean.sd, all.scales$tmeanSD)
  
  q90.tmax <- list.files(outpath2, pattern = paste0("q90_tmax_"), full.names = T)
  q90.tmax <- load_scaled(q90.tmax, all.scales$q90tmax)
  
  DR.mean <- list.files(outpath2, pattern = paste0("DiurnalRange_"), full.names = T)
  DR.mean <- load_scaled(DR.mean, all.scales$meanDR)
  
  distinct_mean <- rast()
  distinct_sd <- rast()
  for(i in 1:length(yr.seq)){
    
    yr <- yr.seq[i]
    if(sample.name=="MobaCluster" & yr == 2007) next()
    idx <- which(names(tmin.mean) == yr)
    tmin.r <- tmin.mean[[idx]]  
    
    idx <- which(names(tmean.sd) == yr)
    tmeanSD.r <- tmean.sd[[idx]]  
    
    idx <- which(names(q90.tmax) == yr)
    q90max.r <- q90.tmax[[idx]]  
    
    idx <- which(names(DR.mean) == yr)
    DR.r <- DR.mean[[idx]]  
    
    n.h <- nlyr(tmin.r)  # number of heights
    
    clim_stack <- c(
      tmin.r,
      q90max.r,
      DR.r,
      tmeanSD.r
    )
    
    clim_dist <- app(clim_stack, clim_dist_fun)
    
    names(clim_dist) <- c("mean_climate_distinctness", "sd_climate_distinctness")
    
    distinct_mean <- c(distinct_mean,clim_dist[[1]])
    distinct_sd <- c(distinct_sd,clim_dist[[2]])
  }
  
  writeRaster(distinct_mean, 
              paste0(out.data,"climate_distinctness/",sample.name,"_Mean_allyears.tif"), 
              overwrite = T)
  
  writeRaster(distinct_sd, 
              paste0(out.data,"climate_distinctness/",sample.name,"_SD_allyears.tif"), 
              overwrite = T)
  
  dmean <- mean(distinct_mean, na.rm = T)
  writeRaster(dmean, 
              paste0(out.data,"climate_distinctness/",sample.name,"_Mean.tif"), 
              overwrite = T)
  
  #' By forest layers:
  pca.path <- paste0(out.data,"pca_clusters/")
  pca <- rast(paste0(pca.path, sample.name, "_pca_Clust4.tif"))
  in.path <- paste0(out.data,sample.name,"/veg_zoning/")
  
  forest.types <- unique(values(pca, na.rm = T))
  
  distinct_mean <- rast(paste0(out.data,"climate_distinctness/",sample.name,"_Mean_allyears.tif"))
  names(distinct_mean) <- yr.seq
  
  for(f in 1:length(forest.types)){
    #' For forest type f:
    ft <- forest.types[f]
    print(paste0("Forest type: ",ft))
    zones <- readRDS(paste0(in.path,"VegZoning_FT_",ft,".RDS"))
    for(z in 1:nrow(zones)){
      z1 <- zones[z,]
      
      #' Check for layers:
      idx <- which(names(tmin.mean) == 2004)
      tmin.r <- tmin.mean[[idx]]  
      names(tmin.r) <- model.heights[1:nlyr(tmin.r)]
      h <- as.numeric(names(tmin.r))
      check <- which(h >= z1$zone_start & h < z1$zone_end)
      if (length(check) == 0) next
      
      print(paste0("Zone: ",z1$zone_start," to ", z1$zone_end))
      distinct_mean <- rast()
      distinct_sd <- rast()
      for(i in 1:length(yr.seq)){
        
        yr <- yr.seq[i]
        if(sample.name=="MobaCluster"&yr==2007){next()}
        idx <- which(names(tmin.mean) == yr)
        tmin.r <- tmin.mean[[idx]]  
        names(tmin.r) <- model.heights[1:nlyr(tmin.r)]
        h <- as.numeric(names(tmin.r))
        
        tmin.sub <- tmin.r[[
          h >= z1$zone_start &
            h <  z1$zone_end
        ]]
        
        idx <- which(names(tmean.sd) == yr)
        tmeanSD.r <- tmean.sd[[idx]]  
        names(tmeanSD.r) <- model.heights[1:nlyr(tmeanSD.r)]
        tmeanSD.sub <- tmeanSD.r[[
          h >= z1$zone_start &
            h <  z1$zone_end
        ]]
        
        idx <- which(names(q90.tmax) == yr)
        q90max.r <- q90.tmax[[idx]]  
        names(q90max.r) <- model.heights[1:nlyr(q90max.r)]
        q90max.sub <- q90max.r[[
          h >= z1$zone_start &
            h <  z1$zone_end
        ]]
        
        idx <- which(names(DR.mean) == yr)
        DR.r <- DR.mean[[idx]]  
        names(DR.r) <- model.heights[1:nlyr(DR.r)]
        DR.sub <- DR.r[[
          h >= z1$zone_start &
            h <  z1$zone_end
        ]]
        
        n.h <- nlyr(tmin.r)  # number of heights
        clim_stack <- c(
          tmin.sub,
          q90max.sub,
          DR.sub,
          tmeanSD.sub
        )
        
        clim_dist <- app(clim_stack, clim_dist_fun)
        
        #' Mask to forest type.
        clim_dist[pca != ft] <- NA
        
        names(clim_dist) <- c("mean_climate_distinctness", "sd_climate_distinctness")
        
        distinct_mean <- c(distinct_mean,clim_dist[[1]])
        distinct_sd <- c(distinct_sd,clim_dist[[2]])
      }
      
      writeRaster(distinct_mean, 
                  paste0(out.data,"climate_distinctness/",sample.name,"/FT_",
                         ft,"_Mean_allyears_",
                         z1$zone_start,"_",z1$zone_end,".tif"),
                  overwrite = T)
      writeRaster(distinct_sd, 
                  paste0(out.data,"climate_distinctness/",sample.name,"/FT_",
                         ft,"_SD_allyears_",
                         z1$zone_start,"_",z1$zone_end,".tif"),
                  overwrite = T)
    }
  }
}








