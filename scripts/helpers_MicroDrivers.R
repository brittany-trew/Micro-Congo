#' Manuscript: Vertical microclimate diversity shapes thermal refugia in Congolese rainforests.
#' Date: August 2026
#' Author: Brittany Trew (*brittany.trew@fas.harvard.edu*)
#' Description: This script contains all the helper functions required for processing the LiDAR point clouds.
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' @title Check if Input is a LAS Object
#' @description
#' Internal utility function that throws an error if the input is not a `LAS` object (as defined by the `{lidR}` package). Primarily used for early validation in point cloud processing workflows.
#' @param las Object to be tested. Expected to be a `LAS` object.
#' @return Invisibly returns the input if valid; otherwise, throws an informative error.
#' @noRd
.stopifnotlas <- function(las) {
  name <- deparse(substitute(las))
  if (!inherits(las, "LAS")) {
    stop(
      sprintf("`%s` must be a LAS object, but has class: %s",
              name, paste(class(las), collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(las)
}

#' @title Generate Digital Elevation Model (DEM) from a Point Cloud.
#' @description
#' This internal helper function creates a Digital Elevation Model (DEM) from a LAS point cloud by extracting ground-classified points and applying a specified interpolation algorithm. The DEM is returned as a raster.
#' @param las A `LAS` or `LAScluster` object from the `{lidR}` package. Must contain ground-classified points (classification code 2).
#' @param res Numeric. Desired output raster resolution (in the same units as the LAS projection).
#' @param dem.algorithm Character. Interpolation method to use for rasterizing terrain. Options are `"tin"` (default), `"knnidw"`, or `"kriging"`.
#' @return A `RasterLayer` or `SpatRaster` (depending on your environment) representing the interpolated DEM.
#' @import lidR
#' @noRd
.makeDEM <- function(las, res, dem.algorithm) {
  if(dem.algorithm == "tin") dtm <- rasterize_terrain(
    las, res = res, algorithm = tin(), use_class = 2L)
  if(dem.algorithm == "knnidw") dtm <- rasterize_terrain(
    las, res = res, algorithm = knnidw(), use_class = 2L)
  if(dem.algorithm == "kriging") dtm <- rasterize_terrain(
    las, res = res, algorithm = kriging(), use_class = 2L)
  return(dtm)
}

#' @title Generate Canopy Height Model (CHM) from Normalized Point Cloud
#' @description
#' Internal helper function that computes a Canopy Height Model (CHM) from a height-normalized LAS object. Users can specify one of three algorithms: `"pitfree"` (default), `"p2r"`, or `"dsmtin"`. Optional post-processing (e.g., smoothing) is applied when appropriate.
#' @param las_norm A height-normalized `LAS` object from the `{lidR}` package.
#' @param res Numeric. Desired resolution of the output CHM raster.
#' @param chmA Character. Algorithm to use for CHM generation. Options are `"pitfree"` (default), `"p2r"`, or `"dsmtin"`. If `NULL`, `"pitfree"` is used.
#' @return A `RasterLayer` or `SpatRaster` object representing the Canopy Height Model.
#' @import lidR
#' @importFrom terra focal
#' @noRd
.CanopyModel <- function(las_norm, res, chmA = NULL){
  if(chmA == "pitfree" || is.null(chmA)){
    chm <- rasterize_canopy(las_norm, res, algorithm = pitfree())
  } else if (chmA == "p2r"){
    chm <- rasterize_canopy(las_norm, res, algorithm = p2r(na.fill = tin()))
    w <- matrix(1, 3, 3)
    chm <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
  } else if (chmA == "dsmtin"){
    chm <- rasterize_canopy(las_norm, res, algorithm = dsmtin())
  } else {
    print("CHM Algorithm not recognised.")
  }
  return(chm)
}

#' @title Derive CHMs and DEMs.
#' @description
#' This function splits large LiDAR point cloud files (`.las` or `.laz`) into smaller tiles when their file size exceeds a user-defined threshold (`max.size`). 
#' Each output tile will be approximately `cut.off` GB in size, split along the X-axis with an optional spatial buffer to avoid edge artifacts. 
#' This is useful for managing large files that are too big to process efficiently or that exceed memory constraints. Note, any original files which have been split will be deleted if delete = T (default).
#' @param all.tiles Character vector. Full file paths to `.las` or `.laz` files to be evaluated and potentially split. Use `list.files(..., full.names = TRUE)` to generate.
#' @param defined_proj Numeric. Threshold file size in gigabytes (GB) above which a file will be split. Files smaller than this are ignored.
#' @param outdir Numeric. Approximate target size in gigabytes (GB) for each split tile. Determines the number of output chunks.
#' @param res Numeric. Spatial resolution of outputs (in metres). 
#' @param canopyHeight Numeric. Width of buffer to include on either side of each tile (in the same units as the point cloud's coordinate system). Helps prevent artifacts along tile edges.
#' @param dem.algorithm Character. Interpolation method to use for rasterizing terrain. Options are `"tin"` (default), `"knnidw"`, or `"kriging"`.
#' @import lidR
#' @export
.LASprocess <- function(all.tiles, defined_proj, outdir,res,
                        canopyHeight = NULL, dem.algorithm = "tin"){
  
  #' Loops through multiple tiles if needed.
  for (a in 1:length(all.tiles)) {
    tryCatch({
      #' **Read in .las/.laz tile**
      las <- readLAS(all.tiles[a])
      tile <- sub(
        ".*_([0-9]+)$",
        "\\1",
        tools::file_path_sans_ext(basename(all.tiles[a]))
      )
      print(paste0("Processing tile: ", tile, "... ", a, " of ", length(all.tiles), "."))
      
      #' **Initial Checks**
      .stopifnotlas(las)
      st_crs(las) <- defined_proj
      chk <- las_check(las, print = FALSE)
      
      #' **Ground Classification**
      if ("The point cloud is ground classified" %in% chk$messages) {
        print("The point cloud is ground classified.")
      } else {
        print("The point cloud isn't ground classified. Classifying ground points...")
        las <- classify_ground(las, algorithm = pmf())
      }
      
      #' **Normalise Heights**
      if ("The point cloud is height normalized" %in% chk$messages) {
        print("The point cloud is height normalized.")
        las_norm <- las
      } else {
        print("The point cloud isn't height normalised. Normalising heights...")
        las_norm <- normalize_height(las, tin())
      }
      
      #' **Digital Elevation Model**
      print("Calculating a digital elevation model...")
      dem <- .makeDEM(las, res, dem.algorithm)
      
      path <- paste0(outdir, "dem/")
      dir.create(path, showWarnings = FALSE)
      writeRaster(dem, paste0(path, "dem_", tile, ".tif"), overwrite = TRUE)
      rm(las)  # save memory
      
      #' **Canopy Height Model**
      if (is.null(canopyHeight)) {
        chm <- .CanopyModel(las_norm, res)
        path <- paste0(outdir, "chm/")
        dir.create(path, showWarnings = FALSE)
        writeRaster(chm, paste0(path, "chm_", tile, ".tif"), overwrite = TRUE)
      }
      rm(las_norm)
      gc()
    }, error = function(e) {
      warning(paste("Error processing tile", all.tiles[a], ":", e$message))
    })
    gc() }
}


canopy.height.levelr.fixed <- function(lidar.array)
{
  l.array <- lidar.array$array
  chm.pulses <- array(data = NA, dim = dim(l.array))
  
  for (q in 1:dim(l.array)[2]) {
    for (z in 1:dim(l.array)[3]) {
      
      lidar.col <- l.array[, q, z]
      ground <- lidar.col[1]
      canopy <- lidar.col[2]
      
      if (is.finite(ground)) {
        
        canopy.ht <- canopy - ground
        ch.col <- lidar.col
        ch.col[1] <- 0
        ch.col[2] <- canopy.ht
        
        g.pulse <- which(ch.col[3:length(ch.col)] > 0)
        
        # ** The Main Fix: Handle columns with no positive voxel pulse counts
        if (length(g.pulse) == 0) {
          chm.pulses[, q, z] <- NA
          next
        }
        # ---------------------------------------------
        
        g.ind <- g.pulse[1] + 2
        cht.col <- c(ch.col[1:2], ch.col[g.ind:length(ch.col)])
        
        new.index <- length(cht.col) + 1
        if (new.index > length(ch.col)) {
          chm.pulses[, q, z] <- cht.col
        } else {
          cht.col[new.index:length(ch.col)] <- 0
          chm.pulses[, q, z] <- cht.col
        }
        
      } else {
        chm.pulses[, q, z] <- NA
      }
    }
  }
  
  return.data <- base::list(
    array = chm.pulses,
    x.bin = lidar.array$x.bin,
    y.bin = lidar.array$y.bin,
    z.bin = lidar.array$z.bin
  )
  return(return.data)
}



machorn.lad.fixed <- function (leveld.lidar.array, voxel.height, beer.lambert.constant = NULL) 
{
  voxel.N.pulse <- leveld.lidar.array$array[3:dim(leveld.lidar.array$array)[1], 
                                            , ]
  pulse.accum <- array(0, dim = dim(voxel.N.pulse))
  voxel.N.pulse[is.na(voxel.N.pulse)] <- 0
  for (i in (dim(voxel.N.pulse)[1]):1) {
    if (i == (dim(voxel.N.pulse)[1])) {
      pulse.accum[i, , ] <- voxel.N.pulse[i, , ]
    }
    else {
      pulse.accum[i, , ] <- pulse.accum[i + 1, , ] + voxel.N.pulse[i, 
                                                                   , ]
    }
  }
  pulse.all <- array(rep(pulse.accum[1, , ], each = dim(voxel.N.pulse)[1]), 
                     dim = dim(voxel.N.pulse))
  shots.through <- pulse.all - pulse.accum
  shots.in <- array(dim = dim(voxel.N.pulse))
  shots.in[1:(dim(voxel.N.pulse)[1] - 1), , ] <- shots.through[2:(dim(voxel.N.pulse)[1]), 
                                                               , ]
  shots.in[(dim(voxel.N.pulse)[1]), , ] <- pulse.accum[1, , 
  ]
  k <- beer.lambert.constant
  dz <- voxel.height
  if (is.null(k)) {
    rLAD <- log(shots.in/shots.through) * (1/dz)
    rLAD[is.infinite(rLAD) | is.nan(rLAD)] <- NA
  }
  else {
    print("MacArthur-Horn constant is set! k = ")
    print(k)
    rLAD <- log(shots.in/shots.through) * (1/(k * dz))
    rLAD[is.infinite(rLAD) | is.nan(rLAD)] <- NA
  }
  for (r in 1:dim(rLAD)[2]) {
    for (c in 1:dim(rLAD)[3]) {
      na.cut <- ceiling(
        leveld.lidar.array$array[2,r,c] / voxel.height
      ) + 2 # Fix is here so that canopy isn't cut according to 1m voxels only.
      if (is.na(na.cut) == FALSE) {
        if (na.cut > dim(rLAD)[1]) {
        }
        else {
          rLAD[na.cut:dim(rLAD)[1], r, c] <- NA
        }
      }
      else {
      }
    }
  }
  out <- list()
  out$rLAD <- rLAD
  out$shots.in <- shots.in
  gc()
  remove(leveld.lidar.array)
  remove(pulse.accum)
  remove(voxel.N.pulse)
  gc()
  return(out)
}

