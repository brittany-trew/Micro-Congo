#' Manuscript: Vertical microclimate diversity shapes thermal refugia in Congolese rainforests.
#' Date: August 2026
#' Author: Brittany Trew (*brittany.trew@fas.harvard.edu*)
#' Description: Converts LiDAR point cloud tiles into a Canopy Height Model, Digital Elevation Model
#' and calculates Plant Area Density (PAD) by voxels.
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R"))

#' Step 1: Convert tiled point clouds into DEM and CHM tiles.
#' Checks for .las files and derives DEM & CHM...
.LASprocess(all.tiles, defined_proj, outdir = head.path, res)

#' Step 2: Merge CHM tiles.
chm <- list.files(paste0(head.path,"/chm/"), full.names = T)
chm_list <- lapply(chm, rast)
chm_mosaic <- do.call(mosaic, chm_list)
plot(chm_mosaic)
writeRaster(chm_mosaic, paste0(head.path,"chm.tif"), overwrite = T)

#' Step 3: Merge DEM tiles.
dem <- list.files(paste0(head.path,"/dem/"), full.names = T)
dem_list <- lapply(dem, rast)
dem_mosaic <- do.call(mosaic, dem_list)
plot(dem_mosaic)
writeRaster(dem_mosaic, paste0(head.path,"dem.tif"), overwrite = T)


#' Estimate PAD by voxel.
epsg_code <- as.numeric(sub("EPSG:", "", defined_proj))
for (a in 1:length(all.tiles)) {
  tile <- sub(
    ".*_([0-9]+)$",
    "\\1",
    tools::file_path_sans_ext(basename(all.tiles[a]))
  )  
  print(paste0("Processing tile: ", tile, "... ", a, " of ", length(all.tiles), "."))
  las <- readLAS(all.tiles[a])
  chk <- las_check(las, print = FALSE)
  #' Convert las file to a voxelised lidar array. 
  tryCatch({
    if ("The point cloud is ground classified" %in% chk$messages) {
      print("The point cloud is ground classified. Organising voxels...")
      las.voxel <- laz.to.array(
        all.tiles[a],
        voxel.resolution = res,
        z.resolution = dzd,
        use.classified.returns = TRUE)
    } else {
      print("The point cloud isn't ground classified. Classifying ground points...")
      las <- classify_ground(
        las,
        algorithm = pmf()
      )
      # Write classified LAS temporarily because laz.to.array() requires a file
      tmp.las <- tempfile(fileext = ".las")
      writeLAS(las, tmp.las)
      
      print("Organising voxels...")
      las.voxel <- laz.to.array(
        tmp.las,
        voxel.resolution = res,
        z.resolution = dzd,
        use.classified.returns = TRUE
      )
      # Remove temporary file
      unlink(tmp.las)
    }
    
    
    #' Level the voxelized array to mimic a canopy height model
    level.canopy <- canopy.height.levelr.fixed(lidar.array = las.voxel)
    
    print("Estimating LAD...")
    # Estimate LAD for each voxel in leveled array
    lad.estimates <- machorn.lad(leveld.lidar.array = level.canopy, 
                                 voxel.height = dzd, 
                                 beer.lambert.constant = k.forest)
    
    
    # Convert the LAD array into a single raster stack
    lad.raster <- lad.array.to.raster.stack(lad.array = lad.estimates, 
                                            laz.array = las.voxel, 
                                            epsg.code = epsg_code)
    
    #' convert to Spat-raster and save.
    lad.r <- rast(lad.raster)
    terra::crs(lad.r) <- defined_proj
    
    heights <- seq(
      from = dzd,
      by = dzd,
      length.out = terra::nlyr(lad.r)
    )
    names(lad.r) <- paste0("PAD_", heights, "m")
    
    writeRaster(lad.r, paste0(pad.path,"pad_voxels_",tile,".tif"), overwrite = T)
    saveRDS(lad.estimates, paste0(pad.path,"pad_voxels_",tile,".RDS"))
    
  }, error = function(e) {
    
    message(paste0("⚠️ Tile ", tile, " skipped due to error:"))
    message(e$message)
    
  })
}
