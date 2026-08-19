#' Final workflow: Driving microclimate models with LiDAR-derived vegetation metrics. 
#' LiDAR Pre-processing steps. 
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R"))
#source(paste0(scripts.path,"LandCovertoExtCo.R")) # loads extinction coefficient map.

#' Checks for .las files and derives DEM & CHM...
.LASprocess(all.tiles, defined_proj, outdir = head.path, res)
#       las <- decimate_points(las, homogenize(2))
pad.path <- paste0(head.path,"pad/")
dir.create(pad.path, showWarnings = F)

#' Extinction Coefficient:
k.forest <- 0.5

chm <- list.files(paste0(head.path,"/chm/"), full.names = T)
chm_list <- lapply(chm, rast)
chm_mosaic <- do.call(mosaic, chm_list)
plot(chm_mosaic)
#unlink(chm)
writeRaster(chm_mosaic, paste0(head.path,"chm.tif"), overwrite = T)

dem <- list.files(paste0(head.path,"/dem/"), full.names = T)
dem_list <- lapply(dem, rast)
dem_mosaic <- do.call(mosaic, dem_list)
plot(dem_mosaic)
#unlink(dem)
writeRaster(dem_mosaic, paste0(head.path,"dem.tif"), overwrite = T)

#' Estimate PAD by voxel.
for (a in 1:length(all.tiles)) {
  tile <- sub(".*_(\\d+)\\.(?:las|laz)$", "\\1", all.tiles[a])
  print(paste0("Processing tile: ", tile, "... ", a, " of ", length(all.tiles), "."))
  
  #' Convert las file to a voxelised lidar array. 
  require(canopyLazR)
  tryCatch({
    
    message("Organising voxels...")
    
    las.voxel <- laz.to.array(
      all.tiles[a],
      voxel.resolution = res,
      z.resolution = dzd,
      use.classified.returns = TRUE)
      
    #' Level the voxelized array to mimic a canopy height model
    level.canopy <- canopy.height.levelr.new(lidar.array = las.voxel)
    
    print("Estimating LAD...")
    # Estimate LAD for each voxel in leveled array
    lad.estimates <- machorn.lad(leveld.lidar.array = level.canopy, 
                                 voxel.height = dzd, 
                                 beer.lambert.constant = k.forest)
  
    # Convert the LAD array into a single raster stack
    lad.raster <- lad.array.to.raster.stack(lad.array = lad.estimates, 
                                            laz.array = las.voxel, 
                                            epsg.code = defined_proj)
    
    #' convert to Spat-raster and save.
    lad.r <- rast(lad.raster)
    plot(lad.r[[c(1:2)]])
    crs(lad.r) <- defined_proj
    ext(lad.r) <- ext(lad.raster)
    writeRaster(lad.r, paste0(pad.path,"pad_voxels_",tile,".tif"), overwrite = T)
    saveRDS(las.voxel, paste0(pad.path,"pad_voxels_",tile,".RDS"))
  
  }, error = function(e) {
    
    message(paste0("⚠️ Tile ", tile, " skipped due to error:"))
    message(e$message)
    
  })
}

chm_mosaic <- rast(paste0(head.path,"chm.tif"))
#' Estimate PAI (cumulative LAD) from Leaf Area Density (m2m3)
all.pad <- list.files(paste0(head.path,"pad"), pattern = ".tif", full.names = T)
pad.r <- lapply(all.pad, rast)
pad.r <- do.call(mosaic,pad.r)
pad.r <- terra::resample(pad.r, chm_mosaic)

n <- nlyr(pad.r)
z0 <- (0:(n-1)) * dzd # bottom of each height (Z 0)
keep <- which(z0 >= dzd)
pad.aboveG <- pad.r[[keep]] # remove ground layer 0-0.5m
plot(pad.aboveG[[1]])

nG  <- nlyr(pad.aboveG)
z0G <- z0[keep]

chm <- clamp(chm_mosaic, lower=0, upper=max(z0G) + dzd, values=TRUE)
chm_max <- as.numeric(global(chm_mosaic, "max", na.rm=TRUE)[1,1])
max_start <- floor(chm_max / dzd) * dzd
start_heights <- z0G[z0G <= max_start]

if(site.location == "Gorongosa"){start_heights <- start_heights[start_heights <= 15]}
if(site.location == "Selati"){start_heights <- start_heights[start_heights <= 15]}
if(site.location == "Kruger"){start_heights <- start_heights[start_heights <= 15]}
if(site.location == "Cederberg"){start_heights <- start_heights[start_heights <= 15]}

#' Calc PAI at every height to canopy.
PAI.list <- vector("list", length(start_heights))
for (i in seq_along(start_heights)) {
  
  h0 <- start_heights[i]
  PAI.list[[i]] <- app(c(pad.aboveG, chm), fun = function(v) {
    ch  <- v[nG + 1]
    lad <- v[1:nG]
    if (is.na(ch)) return(NA_real_)
    if (ch <= h0) return(0)
    sum(lad[z0G >= h0 & z0G < ch] * dzd, na.rm = TRUE)
  })
}
PAI.stk <- rast(PAI.list)

mx <- global(PAI.stk, "max", na.rm = TRUE)[,1]
keep_out <- which(mx > 0 & !is.na(mx))
PAI.stk <- PAI.stk[[keep_out]]
start_heights <- start_heights[keep_out]
names(PAI.stk) <- sprintf("PAI_%.1fm_to_canopy", start_heights)
plot(PAI.stk[[1:6]])

#' Adjust for seasonality:
pai.path <- paste0(head.path,"pai/")
dir.create(pai.path, showWarnings = F)

writeRaster(PAI.stk, paste0(pai.path,"pai_noSeasonalAdjust.tif"))

#' Adjust PAI for seasonality to produce 12 months.
if(site.location == "Cederberg"| site.location == "Gorongosa"){
  MODISpath <- paste0(in.data,"MODIS_LAI/",site.location,"/Full_extent/")
}else{
  MODISpath <- paste0(in.data,"MODIS_LAI/",site.location,"/",sample.name,"/")
}

pai_list <- .MODISAdjust(MODISpath, PAI.stk, dzd)

for(i in 1:length(pai_list)){
  hm <- names(pai_list)[[i]]
  print(hm)
  r <- pai_list[[i]]
  plot(r[[1]])
  writeRaster(pai_list[[i]], paste0(pai.path,hm,".tif"), overwrite = T)
}

#' Create fine-scale habitat layers
file1 <- paste0(scripts.path, "LandCover/process_landcover_", sample.name, ".R")
file2 <- paste0(scripts.path, "LandCover/process_landcover_", site.location, ".R")

if (file.exists(file1)) {
  source(file1)
} else if (file.exists(file2)) {
  source(file2)
} else {
  warning("Neither landcover script exists.")
}

source(paste0(scripts.path,"/process_groundmetrics.R")) # processes albedo and soil

#' Define modelling spaces based on habitat.

  
# Calculate max LAD and height of max LAD
epsg_num <- as.numeric(sub("EPSG:", "", defined_proj))

max.lad <- lad.ht.max(lad.array = lad.estimates, 
                      laz.array = las.voxel, 
                      ht.cut = 0, 
                      epsg.code = epsg_num)

# Calculate the ratio of filled and empty voxels in a given column of the canopy
# (1) ratio of voxels in a given column that contain a LAD estimate 
# (2) ratio of voxels in a given column that are empty (i.e. no LAD estimates)
empty.filled.ratio <- canopy.porosity.filled.ratio(lad.array = lad.estimates,
                                                   laz.array = las.voxel,
                                                   ht.cut = 0,
                                                   epsg.code = epsg_num)

# Calculate the volume of filled and empty voxles in a given column of the canopy
empty.filled.volume <- canopy.porosity.filled.volume(lad.array = lad.estimates,
                                                     laz.array = las.voxel,
                                                     ht.cut = 0,
                                                     xy.res = res,
                                                     z.res = dzd,
                                                     epsg.code = epsg_num)

# Calculate the within canopy rugosity
within.can.rugosity <- rugosity.within.canopy(lad.array = lad.estimates,
                                              laz.array = las.voxel,
                                              ht.cut = 0,
                                              epsg.code = epsg_num)

# Calculate the heights of various LAD quantiles
ht.quantiles <- lad.quantiles(lad.array = lad.estimates,
                              laz.array = las.voxel,
                              ht.cut = 0,
                              epsg.code = epsg_num)




