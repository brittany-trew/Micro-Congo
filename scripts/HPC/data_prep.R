#' @title Load and Install Required R Packages
#' @description
#' Loads a vector of R package names. If a package is not already installed, the function installs it from CRAN and then loads it. Useful for reproducible workflows or scripts with many dependencies.
#' @param packages Character vector. Names of the packages to load.
#' @return No return value. Side effect: packages are loaded into the session (and installed if missing).
#' @export
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

#' @title Split Large LAS/LAZ Files.
#' @description
#' This function splits large LiDAR point cloud files (`.las` or `.laz`) into smaller tiles when their file size exceeds a user-defined threshold (`max.size`). 
#' Each output tile will be approximately `cut.off` GB in size, split along the X-axis with an optional spatial buffer to avoid edge artifacts. 
#' This is useful for managing large files that are too big to process efficiently or that exceed memory constraints. Note, any original files which have been split will be deleted if delete = T (default).
#' @param las.files Character vector. Full file paths to `.las` or `.laz` files to be evaluated and potentially split. Use `list.files(..., full.names = TRUE)` to generate.
#' @param max.size Numeric. Threshold file size in gigabytes (GB) above which a file will be split. Files smaller than this are ignored.
#' @param cut.off Numeric. Approximate target size in gigabytes (GB) for each split tile. Determines the number of output chunks.
#' @param buffer Numeric. Width of buffer to include on either side of each tile (in the same units as the point cloud's coordinate system). Helps prevent artifacts along tile edges.
#' @import lidR
#' @export
.Splitlas <- function(las.files, max.size, cut.off, buffer, delete = T){
  
  file.sizes <- sapply(las.files, function(f) file.info(f)$size / (1024^3))  # size in GB
  names(file.sizes) <- basename(las.files)
  too.big <- file.sizes > max.size
  big.tiles <- las.files[too.big]
  
  if (length(big.tiles) == 0) {
    message("No tiles exceed the size threshold. Nothing to split.")
    return(invisible(NULL))
  }
  print("The following files are too big and will be split:")
  print(big.tiles)
  
  for (f in 1:length(big.tiles)) {
    
    file <- big.tiles[f]
    print(paste0("Splitting: ", file))
    size_gb <- file.info(file)$size / (1024^3)
    if (is.na(size_gb)) next
    
    las <- readLAS(file)
    
    n_chunks <- ceiling(size_gb / cut.off)
    bounds <- st_bbox(las)  # xmin, xmax, ymin, ymax
    xmin <- bounds["xmin"]
    xmax <- bounds["xmax"]
    xbreaks <- seq(xmin, xmax, length.out = n_chunks + 1)
    
    parts <- seq(1,n_chunks,1)
    
    for (i in 1:length(parts)) {
      x_min <- max(xbreaks[i]   - buffer, xmin)
      x_max <- min(xbreaks[i+1] + buffer, xmax)
      las_part <- filter_poi(las, X >= xbreaks[i] & X < xbreaks[i+1])

      outfile <- file.path(dirname(file),
                           paste0(tools::file_path_sans_ext(basename(file)),
                                  sprintf("%02d", i),".las"))
      
      writeLAS(las_part, outfile)
    }
    if(delete) file.remove(file)
  }
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
      tile <- sub(".*_(\\d+)\\.(?:las|laz)$", "\\1", all.tiles[a])
      print(paste0("Processing tile: ", tile, "... ", a, " of ", length(all.tiles), "."))
      
      #' **Initial Checks**
      .stopifnotlas(las)
      st_crs(las) <- defined_proj
      ext <- tryCatch({
        terra::ext(las)
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(ext)) next
      chk <- las_check(las, print = FALSE)
      
      #' **Ground Classification**
      if (chk$messages[4] == "The point cloud is ground classified") {
        print("The point cloud is ground classified.")
      } else {
        print("The point cloud isn't ground classified. Classifying ground points...")
        las <- classify_ground(las, algorithm = pmf()) 
      }
      
      #' **Digital Elevation Model**
      print("Calculating a digital elevation model...")
      dem <- .makeDEM(las, res, dem.algorithm)
      
      path <- paste0(outdir, "dem/")
      dir.create(path, showWarnings = FALSE)
      writeRaster(dem, paste0(path, "dem_", tile, ".tif"), overwrite = TRUE)
      
      #' **Normalise Height**
      if (chk$messages[5] == "The point cloud is not height normalized") {
        print("Normalising heights...")
        las_norm <- normalize_height(las, tin())
      } else {
        print("Point cloud already height normalised")
        las_norm <- las
      }
      rm(las)  # save memory
      
      #' **Canopy Height Model**
      if (is.null(canopyHeight)) {
        chm <- .CanopyModel(las_norm, res)
        path <- paste0(outdir, "chm/")
        dir.create(path, showWarnings = FALSE)
        writeRaster(chm, paste0(path, "chm_", tile, ".tif"), overwrite = TRUE)
      }
      
    }, error = function(e) {
      warning(paste("Error processing tile", all.tiles[a], ":", e$message))
    })
  }
}


vegpfromhab2 <- function (habitats, hgts = NA, pai = NA, lat = NA, long = NA, 
          tme = NA, clump0 = TRUE) 
{
  .poparray <- function(a, sel, v) {
    for (i in 1:length(v)) {
      m <- a[, , i]
      m[sel] <- v[i]
      a[, , i] <- m
    }
    a
  }
  if (class(habitats)[1] == "PackedSpatRaster") 
    habitats <- rast(habitats)
  if (class(lat) == "logical") {
    h_latlong <- project(habitats, "EPSG:4326")
    e <- ext(h_latlong)
    lat <- (e$ymin + e$ymax)/2
    long <- (e$ymin + e$ymax)/2
  }
  if (class(tme)[1] == "logical") {
    current_year <- as.numeric(format(Sys.time(), "%Y"))
    hiy <- ifelse(current_year%%4 == 0, 366 * 24, 365 * 24) - 
      1
    tme <- as.POSIXlt(c(0:hiy), origin = paste0(current_year, 
                                                "-01-01 00:00"), tz = "UTC")
  }
  m <- .is(habitats)
  uh <- unique(as.vector(m))
  uh <- uh[is.na(uh) == F]
  pte <- mean(pai, na.rm = TRUE)
  if (is.na(pte)) {
    paii <- .paifromhabitat(1, lat, long, tme)
    pai <- array(NA, dim = c(dim(m), length(paii)))
  }
  x <- m
  gsmax <- m
  leafr <- m
  leafd <- m
  hgt <- m
  for (i in uh) {
    sel <- which(m == i)
    if (is.na(pte)) {
      paii <- .paifromhabitat(i, lat, long, tme)
      pai <- .poparray(pai, sel, paii)
    }
    vegi <- microclimf:::.onehab(i)
    x[sel] <- vegi$x
    gsmax[sel] <- vegi$gsmax
    leafr[sel] <- vegi$leafr
    leafd[sel] <- vegi$leafd
    hgt[sel] <- vegi$hgt
  }
  clump <- pai * 0
  pai <- .rast(pai, habitats)
  leaft <- 0.5 * leafr
  if (clump0 == F) {
    for (i in 1:dim(pai)[3]) clump[, , i] <- clumpestimate(hgt, 
                                                           leafd, .is(pai)[, , i])
    clump <- .rast(clump, habitats)
  }
  if (class(hgts) == "logical") {
    hgt <- .rast(hgt, habitats)
  }
  else hgt <- hgts
  x <- .rast(x, habitats)
  gsmax <- .rast(gsmax, habitats)
  leafr <- .rast(leafr, habitats)
  leaft <- .rast(leaft, habitats)
  leafd <- .rast(leafd, habitats)
  vegp <- list(pai = wrap(pai), hgt = wrap(hgt), x = wrap(x), 
               gsmax = wrap(gsmax), leafr = wrap(leafr), clump = wrap(clump), 
               leafd = wrap(leafd), leaft = wrap(leaft))
  class(vegp) <- "vegparams"
  return(vegp)
}

#' Version of machorn.lad() from canopyLazR which takes varying extinction coefficient.
machorn.lad.kMap <- function (leveld.lidar.array, voxel.height, 
                              beer.lambert.constant = NULL, k.map = NULL) 
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
  dz <- voxel.height
  if (!is.null(k.map)) {
    # voxel-specific k
    rLAD <- array(NA, dim = dim(shots.in))
    for (r in 1:dim(rLAD)[2]) {
      for (c in 1:dim(rLAD)[3]) {
        kval <- k.map[r, c]
        if (!is.na(kval) && kval > 0) {
          rLAD[, r, c] <- log(shots.in[, r, c] / shots.through[, r, c]) * (1/(kval * dz))
        }
      }
    }
  } else if (is.null(beer.lambert.constant)) {
    rLAD <- log(shots.in/shots.through) * (1/dz)
  } else {
    message("MacArthur-Horn constant is set! k = ", beer.lambert.constant)
    rLAD <- log(shots.in/shots.through) * (1/(beer.lambert.constant * dz))
  }
  rLAD[is.infinite(rLAD) | is.nan(rLAD)] <- NA
  
  for (r in 1:dim(rLAD)[2]) {
    for (c in 1:dim(rLAD)[3]) {
      na.cut <- ceiling(leveld.lidar.array$array[2, r, 
                                                 c]) + 2
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


lad.array.to.raster.stack <- function (lad.array, laz.array, epsg.code) 
{
  lad.testr <- lad.array$rLAD
  lad.mat <- list()
  for (i in seq_along(1:dim(lad.testr)[1])) {
    m.lad <- matrix(data = NA, nrow = dim(lad.testr)[2], 
                    ncol = dim(lad.testr)[3])
    for (q in seq_along(1:dim(lad.testr)[3])) {
      m.lad[, q] <- lad.testr[i, , q]
    }
    lad.mat[[i]] <- m.lad
  }
  raw.lad.rasters <- list()
  for (f in seq_along(1:length(lad.mat))) {
    crs.proj <- base::paste0("+init=epsg:", epsg.code)
    lad.raster <- raster::raster(lad.mat[[f]], xmn = laz.array$x.bin[1], 
                                 xmx = laz.array$x.bin[length(laz.array$x.bin)], ymn = laz.array$y.bin[1], 
                                 ymx = laz.array$y.bin[length(laz.array$y.bin)], crs = crs.proj)
    lad.raster.flip <- flip(lad.raster, direction = "y")
    raw.lad.rasters[[f]] <- lad.raster.flip
  }
  lad.rasters <- do.call(raster::stack, raw.lad.rasters)
  lad.ras <- raster::dropLayer(lad.rasters, 1)
  names(lad.ras) <- paste0("m.", rep(1:nlayers(lad.ras)))
  return(lad.ras)
}


canopy.height.levelr <- function (lidar.array) 
{
  l.array <- lidar.array$array
  chm.pulses <- array(data = NA, dim = dim(l.array))
  for (q in 1:dim(l.array)[2]) {
    for (z in 1:dim(l.array)[3]) {
      lidar.col <- l.array[, q, z]
      ground <- lidar.col[1]
      canopy <- lidar.col[2]
      if (ground > 0 & !is.na(ground)) {
        canopy.ht <- canopy - ground
        returns <- lidar.col[3:length(lidar.col)]
        ch.col <- lidar.col
        ch.col[1] <- 0
        ch.col[2] <- canopy.ht
        g.pulse <- which(ch.col[3:length(ch.col)] > 0)
        g.ind <- g.pulse[1] + 2
        cht.col <- c(ch.col[1:2], ch.col[g.ind:length(ch.col)])
        new.index <- length(cht.col) + 1
        if (new.index > length(ch.col)) {
          chm.pulses[, q, z] <- cht.col
        }
        else {
          cht.col[new.index:length(ch.col)] <- 0
          chm.pulses[, q, z] <- cht.col
        }
      }
      else {
        chm.pulses[, q, z] <- NA
      }
    }
  }
  return.data <- base::list(array = chm.pulses, x.bin = lidar.array$x.bin, 
                            y.bin = lidar.array$y.bin, z.bin = lidar.array$z.bin)
  gc()
  remove(lidar.array)
  gc()
  return(return.data)
}

laz.to.array <- function (laz.file.path, voxel.resolution, z.resolution, use.classified.returns) 
{
  laz.data <- rlas::read.las(laz.file.path)
  laz.xyz.table <- as.data.frame(cbind(x = laz.data$X, y = laz.data$Y, 
                                       z = laz.data$Z, class = laz.data$Classification))
  gc()
  rm(laz.data)
  gc()
  stat.q <- summary(laz.xyz.table[, 3])
  t.lower <- stat.q[2] - 1.5 * (stat.q[5] - stat.q[2])
  t.upper <- stat.q[5] + 1.5 * (stat.q[5] - stat.q[2])
  laz.xyz <- laz.xyz.table[(laz.xyz.table[, 3] >= t.lower) & 
                             (laz.xyz.table[, 3] <= t.upper), ]
  gc()
  rm(laz.xyz.table)
  rm(stat.q)
  rm(t.lower)
  rm(t.upper)
  gc()
  x.range.raw <- range(laz.xyz[, 1], na.rm = T)
  y.range.raw <- range(laz.xyz[, 2], na.rm = T)
  z.range.raw <- range(laz.xyz[, 3], na.rm = T)
  x.y.grain <- voxel.resolution
  z.grain <- z.resolution
  x.range <- c(floor(x.range.raw[1]/x.y.grain) * x.y.grain, 
               ceiling(x.range.raw[2]/x.y.grain) * x.y.grain)
  y.range <- c(floor(y.range.raw[1]/x.y.grain) * x.y.grain, 
               ceiling(y.range.raw[2]/x.y.grain) * x.y.grain)
  z.range <- c(floor(z.range.raw[1]/z.grain) * z.grain, ceiling(z.range.raw[2]/z.grain) * 
                 z.grain)
  x.bin <- seq(x.range[1], x.range[2], x.y.grain)
  y.bin <- seq(y.range[1], y.range[2], x.y.grain)
  z.bin <- seq(z.range[1], z.range[2], z.grain)
  x.cuts <- round((laz.xyz[, 1] - x.bin[1] + x.y.grain/2)/x.y.grain)
  y.cuts <- round((laz.xyz[, 2] - y.bin[1] + x.y.grain/2)/x.y.grain)
  z.cuts <- round((laz.xyz[, 3] - z.bin[1] + z.grain/2)/z.grain)
  y.cuts.dec <- y.cuts/(length(y.bin))
  x.y.cuts <- x.cuts + y.cuts.dec
  gc()
  rm(x.cuts)
  rm(y.cuts)
  gc()
  x.y.levels <- as.data.frame(expand.grid(1:(length(x.bin) - 
                                               1), (1:(length(y.bin) - 1))/length(y.bin)))
  colnames(x.y.levels) <- c("x.level", "y.level")
  x.y.levels <- x.y.levels[order(x.y.levels[, "x.level"], x.y.levels[, 
                                                                     "y.level"]), ]
  x.y.levels.char <- as.character(x.y.levels[, "x.level"] + 
                                    x.y.levels[, "y.level"])
  x.y.cuts.factor <- factor(x.y.cuts, levels = x.y.levels.char)
  xyz.matrix <- as.data.frame(matrix(NA, nrow = length(x.y.levels.char), 
                                     ncol = 4, dimnames = list(NULL, c("x", "y", "z", "class"))))
  xyz.table <- rbind(laz.xyz, xyz.matrix)
  gc()
  rm(laz.xyz)
  rm(xyz.matrix)
  gc()
  z.index <- c(z.cuts, rep(NA, length(x.y.levels.char)))
  x.y.index <- c(x.y.cuts.factor, as.factor(x.y.levels.char))
  lidar.table <- data.frame(xyz.table, x.y.index = x.y.index, 
                            z.index = z.index)
  gc()
  rm(xyz.table)
  gc()
  lidar.array.populator.classified <- function(x) {
    ground.pts <- x[x[, 4] == 2, ]
    as.numeric(c(quantile(ground.pts[, 3], prob = 0, na.rm = TRUE), 
                 quantile(x[, 3], prob = 1, na.rm = TRUE), (table(c(x[, 
                                                                      6], 1:length(z.bin) - 1)) - 1)))
  }
  lidar.array.populator.lowest <- function(x) {
    ground.pts <- x[x[, 4] == 2, ]
    as.numeric(c(quantile(x[, 3], prob = 0, na.rm = TRUE), 
                 quantile(x[, 3], prob = 1, na.rm = TRUE), (table(c(x[, 
                                                                      6], 1:length(z.bin) - 1)) - 1)))
  }
  if (use.classified.returns == TRUE) {
    print("Using classified LiDAR returns!")
    z.vox.test <- length(dlply(lidar.table[1:2, ], "x.y.index", 
                               lidar.array.populator.classified)[[1]])
    lidar.array <- array(as.vector(unlist(dlply(lidar.table, 
                                                "x.y.index", lidar.array.populator.classified, .progress = "text"))), 
                         dim = c(z.vox.test, (length(y.bin) - 1), (length(x.bin) - 
                                                                     1)))
    return.data <- base::list(array = lidar.array, x.bin = x.bin, 
                              y.bin = y.bin, z.bin = z.bin)
  }
  if (use.classified.returns == FALSE) {
    print("Using unclassified LiDAR returns!")
    z.vox.test <- length(dlply(lidar.table[1:2, ], "x.y.index", 
                               lidar.array.populator.lowest)[[1]])
    lidar.array <- array(as.vector(unlist(dlply(lidar.table, 
                                                "x.y.index", lidar.array.populator.lowest, .progress = "text"))), 
                         dim = c(z.vox.test, (length(y.bin) - 1), (length(x.bin) - 
                                                                     1)))
    return.data <- base::list(array = lidar.array, x.bin = x.bin, 
                              y.bin = y.bin, z.bin = z.bin)
  }
  return(return.data)
  gc()
  remove(lidar.table)
  remove(lidar.array)
  gc()
}


#' @title Derive Plant Area Index (PAI) from Height-Normalized LiDAR Tiles
#' @description
#' Processes multiple `.las` or `.laz` tiles to generate Plant Area Index (PAI) rasters from a height-normalized point cloud. For each vertical window defined in `heights`, PAI is computed using vertical profiles, optionally adjusted using MODIS-derived correction factors.
#' @param all.tiles Character vector. Full paths to `.las` or `.laz` point cloud tiles.
#' @param defined_proj Character or numeric. Coordinate reference system (CRS) to assign to the point cloud (e.g., EPSG code).
#' @param heights Numeric vector. Lower bounds of vertical bins for calculating PAI (e.g., `c(0, 2, 5, 10)`).
#' @param maxHeight Numeric. Upper bound for the tallest vertical bin (e.g., 25 or 30).
#' @param outdir Character. Directory where output PAI rasters will be saved.
#' @param MODISpath Character. File path to MODIS-derived correction factor or adjustment surface to refine PAI estimates.
#' @details
#' For each tile, the LAS file is normalized for height if not already done. Then PAI is calculated in vertical bins defined by `heights` and `maxHeight`. Results are saved as multi-layer rasters, optionally corrected using MODIS.
#' @import lidR
#' @importFrom terra rast ext
#' @export
.LASNorm <- function(all.tiles, defined_proj, heights, maxHeight, k, res, outdir, MODISpath){
  #' **Create out directory**
  pai.path <- paste0(outdir,"/pai/")
  dir.create(pai.path, showWarnings = F)
  
  layers <- setNames(
    lapply(heights, function(minh) c(minh, maxHeight)),
    paste0("band_", heights)
  )
  print(layers)
  for(a in 1:length(all.tiles)){
    #' **Read in .las/.laz tile**
    las <- readLAS(paste0(all.tiles[a]))
    tile <- sub(".*_(\\d+)\\.(?:las|laz)$", "\\1", all.tiles[a])
    print(paste0("Processing tile: ",tile, "... ",a, " of ", length(all.tiles),"."))
    
    #' **Initial Checks**
    .stopifnotlas(las)
    projection(las) <- defined_proj
    ext <- tryCatch({
      terra::ext(las)
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(ext)) next
    chk <- las_check(las, print = FALSE)
    
    #' **Normalise Height**
    if(chk$messages[5] == "The point cloud is not height normalized"){
      print("Normalising heights...")
      las_norm <- normalize_height(las, tin())
    } else {
      print("Point cloud already height normalised")
    }
    
    return(las)
    #' **Deriving Plant Area Index from LAD**
    results <- list()
    for (i in 1:length(layers)) {    
      window <- layers[[i]]
      cat(paste0("Computing fixed PAI layer: ", i, " (", window[1], "-", window[2], "m)...\n"))
      zmin <- as.numeric(window[1])
      print(zmin)
      zxi <- as.numeric(window[2])
      dz <- 0.1
      
      pai_layer <- pixel_metrics(
        las_norm,
        ~ compute_pai_fixed(Z, zmin, dz, k),
        res = res
      )
      
      names(pai_layer) <- paste0(window[1],"m")
      results[[i]] <- pai_layer
      plot(results[[i]])
    }
    
    pai_r <- rast(results)
    .MODISAdjust(MODISpath, pai_r)
  }
}



albedo_process<-function(r, pathin)  {
  lst <- list.files(pathin)
  if (length(lst) > 0) {
    pb <- utils::txtProgressBar(min = 0, max = length(lst) + 1, style = 3)
    fi <- paste0(pathin, lst[1])
    modr <- rast(fi)[["Albedo_WSA_shortwave"]]
    print(summary(values(modr, na.rm = TRUE)))
    # get extent of r
    bbx <- rast(ext(r))
    terra::crs(bbx) <- terra::crs(r)
    bbx <- project(bbx, terra::crs(modr))
    e <- ext(bbx)
    e$xmin <- e$xmin - 1000
    e$xmax <- e$xmax + 1000
    e$ymin <- e$ymin - 1000
    e$ymax <- e$ymax + 1000
    modr <- extend(modr, e)
    modr <- crop(modr, e)
    utils::setTxtProgressBar(pb, 1)
    for (i in 2:length(lst)) {
      fi <- paste0(pathin, lst[i])
      mr <- rast(fi)[["Albedo_WSA_shortwave"]]
      mr <- extend(mr, e)
      mr <- crop(mr, e)
      modr <- c(modr, mr)
      utils::setTxtProgressBar(pb, i)
    }
  } else stop("No files to process!")
  m <- apply3D(as.array(modr))
  albedo <- .rast(m, modr)
  albedo <- project(albedo, terra::crs(r))
  albedo <-crop(albedo, ext(r))
  utils::setTxtProgressBar(pb, i + 1)
  return(albedo)
}


albedo_fromaerial <- function(RGBimage, CIRimage, RGBbandmins = c(620, 495, 450),
                              RGBbandmaxs = c(750, 570, 495),
                              CIRbandmins = c(750, 620, 495),
                              CIRbandmaxs = c(900, 750, 570)) {
  # Check resolutions and resample to coarser resolutions if not matching
  res1 <- res(RGBimage)[1]
  res2 <- res(CIRimage)[1]
  if (res1 > res2) CIRimage <- resample(CIRimage, RGBimage)
  if (res2 > res1) RGBimage <- resample(RGBimage, CIRimage)
  # Check extents and intersect
  e1 <- ext(RGBimage)
  e2 <- ext(CIRimage)
  e <- ext(max(e1$xmin, e2$xmin), min(e1$xmax, e2$xmax),
           max(e1$ymin, e2$ymin), min(e1$ymax, e2$ymax))
  RGBimage <- crop(RGBimage, e)
  CIRimage <- crop(CIRimage, e)
  
  RGBimage <- RGBimage / 10000
  CIRimage <- CIRimage / 10000
  # Create weights
  wgt1 <- sum(.Planck(seq(RGBbandmins[1], RGBbandmaxs[1], by = 1))/(10^13))
  wgt2 <- sum(.Planck(seq(RGBbandmins[2], RGBbandmaxs[2], by = 1))/(10^13))
  wgt3 <- sum(.Planck(seq(RGBbandmins[3], RGBbandmaxs[3], by = 1))/(10^13))
  wgt4 <- sum(.Planck(seq(CIRbandmins[1], CIRbandmaxs[1], by = 1))/(10^13))
  wgt5 <- sum(.Planck(seq(CIRbandmins[2], CIRbandmaxs[2], by = 1))/(10^13))
  wgt6 <- sum(.Planck(seq(CIRbandmins[3], CIRbandmaxs[3], by = 1))/(10^13))
  # calculate albedo
  albedo <- (RGBimage[[1]] * wgt1 + RGBimage[[2]] * wgt2 + RGBimage[[3]] * wgt3 +
               CIRimage[[1]] * wgt4 + CIRimage[[2]] * wgt5 + CIRimage[[3]] * wgt6)  /
    (wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6)
  # albedo <- albedo / 255 REMOVED TO REFLECT SENTINEL-2 DATA INPUT
  albedo <- mask(albedo, RGBimage[[1]])
  return(albedo)
}

albedo_adjust<-function(photoalbedo, modisalbedo) {
  # crop modis to match extent of aerial
  if (terra::crs(photoalbedo) != terra::crs(modisalbedo)) modisalbedo<-project(modisalbedo, terra::crs(photoalbedo))
  v1 <- as.vector(crop(modisalbedo, ext(photoalbedo)))
  v1 <-v1[is.na(v1) == FALSE]
  v2 <- as.vector(resample(photoalbedo, modisalbedo))
  v2 <-v2[is.na(v2) == FALSE]
  # logit transform
  lv1 <- log(v1/ (1 - v1))
  lv2 <- log(v2/ (1 - v2))
  # if length of either is 1
  n <- min(length(v1), length(v2))
  if (n == 1) {
    mu <- lv1 - lv2
    lphoto <- log(photoalbedo / (1 - photoalbedo)) + mu
    albedo <- 1 / (1 + exp(-lphoto))
  } else {
    mum <- mean(lv1) - mean(lv2)
    mus <- sd(lv1) / sd(lv2)
    lphoto <- log(photoalbedo / (1 - photoalbedo))
    me <- mean(as.vector(lphoto), na.rm = TRUE)
    lphoto <- ((lphoto - me) * mus) + mum + me
    albedo <- 1 / (1 + exp(-lphoto))
  }
  return(albedo)
}

reflectance_calc <- function(alb, lai, x, plotprogress = TRUE, maxiter = 50, tol = 0.001, bwgt = 0.5){
  e1 <- intersect(ext(lai), ext(alb))
  e <- intersect(e1, ext(x))
  lai <- crop(lai, e)
  alb <- crop(alb, e)
  x <- crop(x, e)
  # check dims
  all_same <- compareGeom(lai, alb, x)
  if (all_same) {
    lref <- x * 0 + 0.5
    gref <- x * 0 + 0.15
    mxdif <- tol * 10
    paim <- as.matrix(lai, wide = TRUE)
    xm <- as.matrix(x, wide = TRUE)
    albm <- as.matrix(alb, wide = TRUE)
    itr <- 1
    while (mxdif > tol) {
      gref2 <- .rast(find_gref(as.matrix(lref, wide = TRUE), paim, xm, albm), x)
      gref2 <- .fillna(gref2, x, zerotoNA = FALSE)
      lref2 <- .rast(find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
      lref2 <- .fillna(lref2, x, zerotoNA=FALSE)
      gref <- bwgt * gref + (1 - bwgt) * gref2
      lref <- bwgt * lref + (1 - bwgt) * lref2
      mxdif1 <- mean(abs(as.vector(gref) - as.vector(gref2)), na.rm = TRUE)
      mxdif2 <- mean(abs(as.vector(lref) - as.vector(lref2)), na.rm = TRUE)
      mxdif <- max(mxdif1, mxdif2)
      if (plotprogress & itr%%3 == 0) {
        tp1 <- paste0("Ground difference from previous: ", round(mxdif1, 4))
        tp2 <- paste0("Leaf difference from previous: ", round(mxdif2, 4))
        par(mfrow = c(1, 2))
        plot(gref, main = tp1, cex.main = 1)
        plot(lref, main = tp2, cex.main = 1)
      }
      itr <- itr+1
      if (itr > maxiter) mxdif <- 0
    }
  } else (stop("Geometries of input rasters do not match"))
  return(list(gref = gref, lref = lref))
}

#' edit: removed rectify()
soildata_download <- function(r, pathdir = getwd(), deletefiles = TRUE) {
  # get bounding box
  dir.create(pathdir,showWarnings=FALSE)
  e<-ext(r)
  subsetx<-paste0("&SUBSET=X(",e$xmin,",",e$xmax,")")
  subsety<-paste0("&SUBSET=Y(",e$ymin,",",e$ymax,")")
  # get subsetting crs
  crs_info <- terra::crs(r)
  epsg_code <- sub(".*ID\\[\"EPSG\",(\\d+)\\].*", "\\1", crs_info)
  subsetcrs<-paste0("&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/",epsg_code)
  outcrs<-paste0("&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/",epsg_code)
  # datasets
  # classes
  vars<-c("bdod","clay","sand","silt")
  dps_mn<-c(0,5,15,30,60,100)
  dps_mx<-c(5,15,30,60,100,200)
  ro<-r
  for (ii in 1:4) {
    for (jj in 1:6) {
      base_url<-paste0("https://maps.isric.org/mapserv?map=/map/",vars[ii],".map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage")
      idform<-paste0("&COVERAGEID=",vars[ii],"_",dps_mn[jj],"-",dps_mx[jj],"cm_Q0.5&FORMAT=GEOTIFF_INT16")
      request_url<-paste0(base_url,idform,subsetx,subsety,subsetcrs,outcrs)
      response <- GET(request_url)
      response
      if (response$status_code == 200) {
        fo<-paste0(pathdir,vars[ii],"_",dps_mn[jj],"_",dps_mx[jj],"cm.tif")
        writeBin(content(response, "raw"), fo)
      } else {
        stop(paste0("Bad query request for ",vars[ii]," ",dps_mn[jj],"-",dps_mx[jj],"cm"))
      }
    }
    dps<-c(5,10,15,30,40,100)
    wgts<-c(1,0.5,0.25,0.125,0.0625,0.03125)*dps
    wgts<-wgts/sum(wgts)
    fi<-paste0(pathdir,vars[ii],"_",dps_mn[1],"_",dps_mx[1],"cm.tif")
    ri<-rast(fi)*wgts[1]
    if(deletefiles) unlink(fi)
    for (jj in 2:6) {
      fi<-paste0(pathdir,vars[ii],"_",dps_mn[jj],"_",dps_mx[jj],"cm.tif")
      ri<-rast(fi)*wgts[jj]+ri
      if(deletefiles) unlink(fi)
    }
    if(ii == 1) ro<-resample(ro,ri)
    ro<-c(ro,ri)
  }
  ro<-ro[[-1]]
  names(ro)<-vars
  return(ro)
}

soildata_downscale <- function(soildata, landcover, water = 80) {
  # get coarse mask
  landcover[landcover == water] <- NA
  msk <- resample(landcover, soildata[[1]])
  bden <- soildata[[1]]
  bden[bden == 0] <- NA
  bden <- .fillna(bden, msk)
  # clay mass fraction
  clay <- soildata[[2]]
  clay[clay == 0] <- NA
  clay <- .fillna(clay, msk)
  # sand mass fraction
  sand <- soildata[[3]]
  sand[sand == 0] <- NA
  sand <- .fillna(sand, msk)
  # silt mass fraction
  silt <- soildata[[4]]
  silt[silt == 0] <- NA
  silt <- .fillna(silt, msk)
  soildata <- c(bden,clay,sand,silt)
  soildataf <- .sharpen(soildata[[1]], landcover)
  for (i in 2:4)  soildataf <- c(soildataf, .sharpen(soildata[[i]], landcover))
  return(soildataf)
}

era5_download<-function(r, tme, credentials, file_prefix, pathout) {
  # set credentials
  sel <- which(credentials$Site == "CDS")
  uid <- credentials$username[sel]
  cds_access_token <- credentials$password[sel]
  ecmwfr::wf_set_key(user = uid,
                     key = cds_access_token)
  # build request
  e<-ext(r)
  ro<-rast(e)
  terra::crs(ro)<-terra::crs(r)
  rll<-project(ro,"EPSG:4326")
  ell<-ext(rll)
  xmn<-floor(ell$xmin*4)/4
  xmx<-ceiling(ell$xmax*4)/4
  ymn<-floor(ell$ymin*4)/4
  ymx<-ceiling(ell$ymax*4)/4
  req <- mcera5::build_era5_request(xmin = xmn, xmax = xmx,
                                    ymin = ymn, ymax = ymx,
                                    start_time = tme[1],
                                    end_time = tme[length(tme)],
                                    outfile_name = file_prefix)
  # check which downloads already exist
  keep<-rep(TRUE,length(req))
  for (i in 1:length(req)) {
    fi<-paste0(pathout,req[[i]]$target)
    if (file.exists(fi)) keep[i]<-FALSE
  }
  s<-which(keep)
  req2<-req[s]
  # download data
  dir.create(pathout,showWarnings=FALSE)
  mcera5::request_era5(request = req2, uid = uid, out_path = pathout, combine = F)
  return(req)
}

era5_process <- function(req, pathin, r, tme, out = "grid", lat = NA, long = NA, resampleout = TRUE) {
  if (class(req) == "logical") {
    files<-list.files(pathin)
    n<-length(files)
  } else {
    n<-length(req)
    files <- ""
    for (i in 1:n) files[i]<-req[[i]]$target
  }
  # Get raster template for nc file
  nc<-paste0(pathin,files[1])
  # Grid process
  climdata<-.extract_clima(nc, r, resampleout)
  if (n > 1) {
    for (i in 2:n) {
      nc<-paste0(pathin,files[i])
      climone<-.extract_clima(nc, r, resampleout)
      for (j in 1:9)  climdata[[j]]<-c(climdata[[j]], climone[[j]])
    }
  }
  tc<-climdata[[1]]
  tmenc<-as.POSIXlt(time(tc),tz="UTC")
  if (length(tmenc) > length(tme)) {
    s<-which(tmenc >= tme[1] & tmenc <= tme[length(tme)])
    for (i in 1:9) {
      ro<-climdata[[i]]
      ro<-ro[[s]]
      climdata[[i]]<-ro
    }
  }
  if (out == "point") {
    if (class(r) != "logical") {
      ll<-.latlongfromrast(r)
      lat <- ll$lat
      long <- ll$long
    }
    xy<-data.frame(x=long,y=lat)
    m<-matrix(0,ncol=10,nrow=length(tme))
    for (i in 1:9) {
      v<-as.numeric(extract(climdata[[i]],xy))
      m[,i+1]<-v[-1]
    }
    climdata<-data.frame(m)
    names(climdata)<-c("obs_time","temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
    climdata$obs_time<-tme
  } else {
    for (i in 1:9)  climdata[[i]]<-wrap(climdata[[i]])
  }
  return(climdata)
}



soildata_gettype <- function(soildata) {
  bden <- soildata[[1]] / 100
  clay <- soildata[[2]] / 1000
  sand <- soildata[[3]] / 1000
  silt <- soildata[[4]] / 1000
  # soil type
  soiltype <- .rast(getsoiltypecpp(as.matrix(bden, wide = TRUE), as.matrix(clay, wide = TRUE),
                                   as.matrix(sand, wide = TRUE), as.matrix(silt, wide = TRUE)), bden)
  return(soiltype)
}


albedo_download<-function(r, tme, pathout, credentials) {
  tsed<-tme[length(tme)]
  tmod<-as.POSIXlt(0,origin="2000-02-18",tz="UTC")
  if (tsed < tmod) stop("No data available prior to 2000-02-18")
  # Download MODIS DATA
  e<-ext(r)
  r2<-rast(e)
  crs(r2)<-crs(r)
  r2<-project(r2, "EPSG:4326")
  e2<-ext(r2)
  st<-substr(as.character(tme[1]),1,10)
  ed<-substr(as.character(tme[length(tme)]),1,10)
  mf<-luna::getNASA("MCD43A3", st, ed, aoi = e2, version = "061", download=FALSE)
  username <- credentials$username[1]
  password <- credentials$password[1]
  if (length(mf) > 0) {
    luna::getNASA("MCD43A3", st, ed, aoi = e2, version = "061", download=TRUE,
                  path=pathout,username=username,password=password,server="LPDAAC_ECS")
  } else {
    stop("No data for specified location or time period")
  }
}


#' get the right UKCP decade
.find_ukcp_decade<-function(collection=c('land-gcm','land-rcm'),startdate,enddate){
  collection<-match.arg(collection)
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")
  rcm_decades<-c('19801201-19901130','19901201-20001130','20001201-20101130','20101201-20201130',
                 '20201201-20301130','20301201-20401130','20401201-20501130','20501201-20601130',
                 '20601201-20701130','20701201-20801130')
  gcm_decades<-c('19791201-19891130','19891201-19991130','19991201-20091130','20091201-20191130',
                 '20191201-20291130','20291201-20391130','20391201-20491130','20491201-20591130',
                 '20591201-20691130','20691201-20791130')
  if(collection=='land-gcm') ukcp_decades<-gcm_decades else ukcp_decades<-rcm_decades
  decade_start<-lubridate::ymd(sapply(strsplit(ukcp_decades,"-"), `[`, 1))
  decade_end<-lubridate::ymd(sapply(strsplit(ukcp_decades,"-"), `[`, 2))
  decades<-ukcp_decades[which(decade_end>startdate & decade_start<enddate)]
  return(decades)
}
#' extract dates from downloaded ukcp18 nc file
.get_ukcp18_dates<-function(ncfile){
  netcdf_data <-ncdf4::nc_open(ncfile)
  time_hours <- ncdf4::ncvar_get(netcdf_data,"time")
  ncdf4::nc_close(netcdf_data)
  years<-floor(time_hours/(360*24))+1970
  months<-ceiling((time_hours-(years-1970)*(360*24)) / (30*24) )
  days<-ceiling( (time_hours-((years-1970)*(360*24)) - ((months-1)*30*24) )/ 24 )
  hours<-time_hours - ((years-1970)*(360*24)) - ((months-1)*30*24) - ((days-1)*24)
  ukcp_dates<-paste(years,sprintf("%02d",months),sprintf("%02d",days),sep="-")
  return(ukcp_dates)
}
#' correct dates form downloaded ukcp18 nc file
.correct_ukcp_dates<-function(ukcp_dates){
  years<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 1))
  months<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 2))
  days<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 3))
  
  # Shift March day by plus one
  sel<-which(months==3)
  days[sel] <- days[sel] + 1
  
  # Shift Feb day by minus 1
  sel <- which(months==2)
  days[sel] <-  days[sel] - 1
  
  # 1st of Feb as 31 Jan
  sel <- which(days == 0)
  days[sel] <- 31
  months[sel] <- 1
  
  # 29 of Feb as 1 Mar
  sel <- which(months == 2 & days == 29)
  months[sel] <- 3
  days[sel] <- 1
  
  # Construct new date variable
  real_dates<-as.POSIXlt(ISOdate(years, months, days))
  real_dates<-trunc(real_dates,"day")
  return(real_dates)
}
#' convert downloaded ukcp18 nc from 360 day to 365 day
.fill_calendar_data<-function(ukcp_r, real_dates, testplot=FALSE, plotdays=c(89:91,242:244)){
  # Assign real dates to layer names and time values
  terra::time(ukcp_r)<-real_dates
  names(ukcp_r)<-real_dates
  # Insert and fill missing date layers of spatrast - could use zoo::na.approx or na.spline
  ukcp_r<-terra::fillTime(ukcp_r)
  if(testplot)  plot(ukcp_r[[plotdays]],main=paste(time(ukcp_r)[plotdays], 'before filling'))
  ukcp_r<-terra::approximate(ukcp_r,method="linear")
  if(testplot)  plot(ukcp_r[[plotdays]],main=paste(time(ukcp_r)[plotdays], 'after filling'))
  return(ukcp_r)
}
#' apply bias correction to UKCP18 data
.biascorrect <- function(hist_obs, hist_mod, fut_mod, rangelims = 1.05, silent = FALSE) {
  if (!inherits(hist_obs, "SpatRaster")) stop("hist_obs must be a SpatRaster")
  if (!inherits(hist_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  if (!inherits(fut_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  # reproject and crop if necessary
  if (crs(hist_mod) != crs(hist_obs)) hist_mod<-project(hist_mod,hist_obs)
  if (crs(fut_mod) != crs(hist_obs)) fut_mod<-project(fut_mod,hist_obs)
  hist_mod<-resample(hist_mod,hist_obs)
  fut_mod<-resample(fut_mod,hist_obs)
  # Convert to arrays
  a1<-.is(hist_obs)
  a2<-.is(hist_mod)
  a3<-.is(fut_mod)
  # Mask out any cells that are missing
  msk1<-apply(a1,c(1,2),mean,na.rm=T)
  msk2<-apply(a2,c(1,2),mean,na.rm=T)
  msk3<-apply(a3,c(1,2),mean,na.rm=T)
  msk<-msk1*msk2*msk3
  # Check whether dataset has more than 1000 entries per time-series and
  # subset if it does
  n<-dim(a1)[3]
  if (n > 1000) s<-sample(0:n,998,replace = FALSE)  # 998 so min and max can be tagged on
  ao<-array(NA,dim=dim(a3))
  # Create array for storing data
  counter<-0
  nn<-dim(a1)[1]*dim(a1)[2]
  if (silent == FALSE) pb<-utils::txtProgressBar(min = 0, max = nn, style = 3)
  for (i in 1:dim(a1)[1]) {
    for (j in 1:dim(a1)[2]) {
      if (silent == FALSE) utils::setTxtProgressBar(pb,counter)
      if (is.na(msk[i,j]) == FALSE) {
        v1<-a1[i,j,]
        v2<-a2[i,j,]
        v1 <- v1[order(v1)]
        v2 <- v2[order(v2)]
        if (n > 1000) {
          v1<-c(v1[1],v1[s],v1[length(v1)]) # tags on min and max value
          v2<-c(v2[1],v2[s],v2[length(v2)]) # tags on min and max value
        }
        # Apply gam
        m1 <- gam(v1~s(v2))
        v3<-a3[i,j,]
        xx <- as.numeric(predict.gam(m1, newdata = data.frame(v2 = v3)))
        if (is.na(rangelims) == FALSE) xx<-rangelimapply(v1, v2, v3, xx)
        ao[i,j,]<-xx
      }
      counter<-counter+1
    }
  }
  ao<-.rast(ao,hist_obs)
  return(ao)
}
#' @title Calculate moving average
#' @noRd
.mav <- function(x, n = 5) {
  y <- stats::filter(x, rep(1 / n, n), circular = TRUE, sides = 1)
  y
}
#' Correct precipitation
.precipcorrect <- function(hist_obs, hist_mod, fut_mod, rangelim = 1.05) {
  if (!inherits(hist_obs, "SpatRaster")) stop("hist_obs must be a SpatRaster")
  if (!inherits(hist_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  if (!inherits(fut_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  # reproject and crop if necessary
  if (crs(hist_mod) != crs(hist_obs)) hist_mod<-project(hist_mod,hist_obs)
  if (crs(fut_mod) != crs(hist_obs)) fut_mod<-project(fut_mod,hist_obs)
  hist_mod<-resample(hist_mod,hist_obs)
  fut_mod<-resample(fut_mod,hist_obs)
  # mask data
  hist_mod<-mask(hist_mod,hist_obs)
  # Calculate observed rainfall total and rain day frac
  rcount<-hist_obs
  rcount[rcount > 0] <-1
  rtot1<-apply(.is(hist_obs),c(1,2),sum)
  tfrac1<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate observed rainfall total and rain day frac
  rcount<-hist_mod
  rcount[rcount > 0] <-1
  rtot2<-apply(.is(hist_mod),c(1,2),sum)
  tfrac2<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate ratios
  mu_tot<-rtot1/rtot2
  mu_frac<-tfrac1/tfrac2
  # Apply correction to modelled future data
  rcount<-fut_mod
  rcount[rcount > 0] <-1
  rtot<-apply(.is(fut_mod),c(1,2),sum)*mu_tot
  tfrac<-(apply(.is(rcount),c(1,2),sum)/dim(rcount)[3])*mu_frac
  # Calculate regional precipitation
  rrain<-apply(.is(fut_mod),3,sum,na.rm=TRUE)
  rr2<-as.numeric(.mav(rrain,10))
  s<-which(rrain==0)
  rrain[s]<-rr2[s]
  s<-which(rrain==0)
  rrain[s]<-0.1
  rrain<-(rrain/max(rrain))*10
  # Adjust rainfall
  fut_mod<-mask(fut_mod,hist_obs)
  a<-as.array(fut_mod)
  mm<-matrix(as.vector(a),nrow=dim(a)[1]*dim(a)[2],ncol=dim(a)[3])
  rtot<-.rast(rtot,hist_obs)
  tfrac<-.rast(tfrac,hist_obs)
  rtot<-as.vector(t(rtot))
  rfrac<-as.vector(t(tfrac))
  mm<-rainadjustm(mm,rrain,rfrac,rtot)
  a2<-array(mm,dim=dim(a))
  a2<-prangelimapply(.is(hist_obs), .is(hist_mod), .is(fut_mod), a2, rangelims)
  out<-.rast(a2,fut_mod)
  return(out)
}


#' Extract daily min, mean, max.
compute.daily <- function(tme, model_sub){
  
  udays <- unique(as.Date(tme))
  ndays <- length(udays)
  
  daily.min <- array(NA, dim = c(dim(model_sub)[c(1:2)], ndays))
  daily.max <- array(NA, dim = dim(daily.min))
  daily.mean <- array(NA, dim = dim(daily.min))
  
  # Loop over days and compute summaries
  for (i in seq_along(udays)) {
    day_indices <- which(as.Date(tme) == udays[i])
    daily_slice <- model_sub[,,day_indices]
    daily_slice[daily_slice==-Inf] <- 0
    
    daily.min[,,i] <- apply(daily_slice, c(1, 2), min, na.rm = TRUE)
    daily.max[,,i] <- apply(daily_slice, c(1, 2), max, na.rm = TRUE)
    daily.mean[,,i] <- apply(daily_slice, c(1, 2), mean, na.rm = TRUE)
  }
  daily.list <- list(
    dailyMin = daily.min, 
    dailyMax = daily.max, 
    dailyMean = daily.mean,
    day = udays)
  return(daily.list)
}

compute.monthly <- function(tme, model_sub) {
  
  # Convert time to monthly format (year-month)
  umonths <- unique(format(tme, "%Y-%m"))
  nmonths <- length(umonths)
  
  monthly.min <- array(NA, dim = c(dim(model_sub)[1:2], nmonths))
  monthly.max <- array(NA, dim = dim(monthly.min))
  monthly.mean <- array(NA, dim = dim(monthly.min))
  
  for (i in seq_along(umonths)) {
    month_indices <- which(format(tme, "%Y-%m") == umonths[i])
    month_slice <- model_sub[,,month_indices]
    month_slice[month_slice == -Inf] <- 0
    
    monthly.min[,,i] <- apply(month_slice, c(1, 2), min, na.rm = TRUE)
    monthly.max[,,i] <- apply(month_slice, c(1, 2), max, na.rm = TRUE)
    monthly.mean[,,i] <- apply(month_slice, c(1, 2), mean, na.rm = TRUE)
  }
  
  monthly.list <- list(
    monthlyMin = monthly.min,
    monthlyMax = monthly.max,
    monthlyMean = monthly.mean,
    month = umonths)
  
  return(monthly.list)
}



bbox.create <- function(bbox.file, las.path){
  if(!file.exists(bbox.file)){
    # Create a LAScatalog (automatically loads all LAS/LAZ files in folder).
    ctg <- readLAScatalog(las.path, recursive = T)
    bbox_vals <- st_bbox(ctg)  # returns xmin, ymin, xmax, ymax
    bbox_vect <- vect(rbind(c(bbox_vals$xmin, bbox_vals$ymin),
                            c(bbox_vals$xmax, bbox_vals$ymin),
                            c(bbox_vals$xmax, bbox_vals$ymax),
                            c(bbox_vals$xmin, bbox_vals$ymax),
                            c(bbox_vals$xmin, bbox_vals$ymin)),
                      type = "polygons",
                      crs = terra::crs(defined_proj))
    writeVector(bbox_vect, bbox.file, overwrite = F)
  } else {bbox_vect <- vect(bbox.file)}
  
  return(bbox_vect)
}


pca_heights <- function(apl){
  #' Combine profiles:
  flat_profiles <- purrr::compact(unlist(apl, recursive = FALSE))
  flat_dfs <- purrr:::map(flat_profiles, ~ {
    df <- as.data.frame(.x)
    names(df) <- c("z","lad")
    df
  })
  #' Pick a z sequence that covers all profiles
  all_z <- seq(
    from = min(purrr::map_dbl(flat_dfs, ~min(.x$z))), 
    to   = max(purrr::map_dbl(flat_dfs, ~max(.x$z))),
    length.out = 100  # choose as fine as you like
  )
  #' For each profile, interpolate lad onto all_z
  profile_mat <- purrr::map_dfr(
    flat_dfs,
    ~{y <- approx(x = .x$z, y = .x$lad, xout = all_z, rule = 2)$y
    tibble(z = all_z, lad = y)},
    .id = "tile"
  ) %>%
    pivot_wider(names_from = z, values_from = lad, names_prefix = "z_") %>%
    tibble::column_to_rownames("tile")
  #' Normalise between 0-1 and compute PCA
  row_sums <- rowSums(profile_mat, na.rm = TRUE)
  profile_prop <- profile_mat / row_sums
  profile_prop[is.na(profile_prop)] <- 0
  pca_prop <- prcomp(profile_prop, center = TRUE, scale. = TRUE)
  
  scores_prop <- as.data.frame(pca_prop$x) %>% 
    tibble::rownames_to_column("tile")
  
  set.seed(123)
  #' Feature space:
  pc_feats <- scores_prop[, c("PC1","PC2","PC3")]
  km <- kmeans(pc_feats, centers = 3, nstart = 25)
  scores_prop$cluster <- factor(km$cluster)
  ggplot(scores_prop, aes(PC1, PC2, color = cluster)) +
    geom_point(size = 2, alpha = 0.8) +
    stat_ellipse(
      aes(fill = cluster),
      geom   = "polygon",
      type   = "norm",     # ← use the normal‐theory ellipse
      level  = 0.8,
      alpha  = 0.2,
      color  = NA,
      show.legend = FALSE
    ) +
    theme_minimal(base_size = 14)
  tileGroups <- cbind.data.frame(scores_prop$tile, scores_prop$cluster)
  colnames(tileGroups) <- c("tile","cluster")
  clustered_profiles <- split(
    flat_dfs,
    tileGroups$cluster[ match(names(flat_dfs), tileGroups$tile) ]
  )
  
  return(clustered_profiles)
  
}



daily.mean <- function(v.list){
  v.unwrapped <- lapply(v.list, unwrap)
  remove_feb29 <- function(r) {
    dts <- time(r)
    keep <- format(dts, "%m-%d") != "02-29"
    r[[which(keep)]]
  }
  v.trimmed <- lapply(v.unwrapped, remove_feb29)
  
  stopifnot(all(sapply(v.trimmed, nlyr) == 365))
  
  daily_means <- lapply(1:365, function(i) {
    layers <- lapply(v.trimmed, \(r) subset(r, i))
    # Defensive filter: remove NULLs just in case
    layers <- layers[!sapply(layers, is.null)]
    if (length(layers) == 0) {
      warning(paste("No layers for day", i))
      return(NULL)
    }
    mean(rast(layers))  # Use rast() instead of stack()
  })
  
  daily_mean_raster <- rast(daily_means)
  names(daily_mean_raster) <- format(seq.Date(as.Date("2001-01-01"), length.out = 365, by = "day"), "%m-%d")
  return(daily_mean_raster)
}



era5_process_rec <- function(req, pathin, r, tme, out = "grid", lat = NA, long = NA, resampleout = TRUE) {
  if (class(req) == "logical") {
    files<-list.files(pathin, recursive = T)
    n<-length(files)
  } else {
    n<-length(req)
    files <- ""
    for (i in 1:n) files[i]<-req[[i]]$target
  }
  # Get raster template for nc file
  nc<-paste0(pathin,files[1])
  # Grid process
  climdata<-.extract_clima(nc, r, resampleout)
  if (n > 1) {
    for (i in 2:n) {
      nc<-paste0(pathin,files[i])
      climone<-.extract_clima(nc, r, resampleout)
      for (j in 1:9)  climdata[[j]]<-c(climdata[[j]], climone[[j]])
    }
  }
  tc<-climdata[[1]]
  tmenc<-as.POSIXlt(time(tc),tz="UTC")
  if (length(tmenc) > length(tme)) {
    s<-which(tmenc >= tme[1] & tmenc <= tme[length(tme)])
    for (i in 1:9) {
      ro<-climdata[[i]]
      ro<-ro[[s]]
      climdata[[i]]<-ro
    }
  }
  if (out == "point") {
    if (class(r) != "logical") {
      ll<-.latlongfromrast(r)
      lat <- ll$lat
      long <- ll$long
    }
    xy<-data.frame(x=long,y=lat)
    m<-matrix(0,ncol=10,nrow=length(tme))
    for (i in 1:9) {
      v<-as.numeric(extract(climdata[[i]],xy))
      m[,i+1]<-v[-1]
    }
    climdata<-data.frame(m)
    names(climdata)<-c("obs_time","temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
    climdata$obs_time<-tme
  } else {
    for (i in 1:9)  climdata[[i]]<-wrap(climdata[[i]])
  }
  return(climdata)
}
