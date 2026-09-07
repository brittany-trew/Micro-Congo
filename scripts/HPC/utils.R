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
  ground <- filter_poi(las, Classification == 2)
  if(dem.algorithm == "tin") dtm <- rasterize_terrain(las, res = res, algorithm = tin())
  if(dem.algorithm == "knnidw") dtm <- rasterize_terrain(las, res = res, algorithm = knnidw())
  if(dem.algorithm == "kriging") dtm <- rasterize_terrain(las, res = res, algorithm = kriging())
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


#' @title Adjust LiDAR-derived PAI using MODIS-based Seasonal Scaling
#' @description
#' Adjusts a multi-layer Plant Area Index (PAI) raster derived from LiDAR using seasonal coefficients and intra-annual variation modeled from MODIS Leaf Area Index (LAI) data. The function rescales LiDAR-derived PAI to match MODIS maxima and applies harmonic regression-derived monthly variation for seasonal correction.
#' @param MODISpath Character. File path to a directory containing MODIS LAI raster layers (e.g., monthly composites).
#' @param pai_r SpatRaster. Multi-layer raster stack of LiDAR-derived PAI (e.g., one layer per canopy height band).
#' @details
#' The function reads all MODIS LAI rasters, scales them, models harmonic coefficients (`a0`, `a1`, `b1`), projects them to the resolution of `pai_r`, rescales PAI values using maximum values, and applies monthly adjustment. One raster per height band per month is written to disk.
#' @return None. Side effect: writes `.tif` files to disk in the `pai.path` directory.
#' @import terra
#' @noRd
.MODISAdjust <- function(MODISpath, pai_r, pai.path,dzd){
  files <- list.files(paste0(MODISpath), full.names = T)
  r.list <- lapply(files, rast) 
  r.stk <- rast(r.list)
  r.stk <- scale_modis_lai(r.stk) # Scale raw MODIS leaf area index data. 
  r.array <- as.array(r.stk)
  nyr <- length(r.list) # number of years
  nmons <- nlyr(r.stk) # number of months included
  
  #' Convert to MODIS LAI to PAI.
  pai_m <- model_PAI(r.array, nyr, nmons)
  pai_modis <- rast(pai_m$pai_MODIS, crs = terra::crs(r.stk), ext = ext(r.stk))
  
  #' **Derive intra-annual variation**
  coef.r <- rast(pai_m$coef_array, ext = ext(r.stk), crs = terra::crs(r.stk))
  names(coef.r) <- c("a0", "a1", "b1")
  
  coef.rf <- project(coef.r, terra::crs(pai_r))
  coef.rf <- terra::resample(coef.rf, pai_r, method = "bilinear")
  a0_fine <- coef.rf$a0
  a1_fine <- coef.rf$a1
  b1_fine <- coef.rf$b1
  omega <- 2 * pi / 12
  
  fine_predicted_stack <- rast()
  for (m in 1:12) {
    month_layer <- a0_fine + a1_fine * sin(omega * m) + b1_fine * cos(omega * m)
    fine_predicted_stack <- c(fine_predicted_stack, month_layer)
  }
  names(fine_predicted_stack) <- paste0("Month_", 1:12)
  
  #' **Additional range adjustment**
  pai_modis <- resample(pai_modis, pai_r)
  
  # Percentile for scaling
  modis_p95 <- app(pai_modis, fun = function(x) quantile(x, 0.99, na.rm=TRUE))
  pai_p95   <- app(pai_r,     fun = function(x) quantile(x, 0.99, na.rm=TRUE))
  scale_factor <- modis_p95 / pai_p95
  scale_factor[is.infinite(scale_factor)] <- NA
  scale_factor[scale_factor < 0] <- NA
  scale_factor <- clamp(scale_factor, 0.3, 3)   # guard against large multipliers
  
  pai_r_scaled <- pai_r * scale_factor
  
  #' **Model monthly PAI in LiDAR**
  pai_list <- list()
  for(i in 1:nlyr(pai_r_scaled)){
    pai_temp <- pai_r_scaled[[i]]
    scaled_stack <- scale_pai(pai_temp, coef.rf, input_month, predicted = fine_predicted_stack)
    pai_list[[i]] <- scaled_stack
    names(pai_list)[i] <- names(pai_temp)
  }
  
  for(i in 1:length(pai_list)){
    height <- sprintf("%02dm", i * dzd)
    writeRaster(pai_list[[i]], paste0(pai.path,"pai_",height,".tif"), overwrite = T)
  }
}


#' @title Scale MODIS Leaf Area Index (LAI) Data
#' @description
#' Scales raw MODIS LAI data from integer format to real values. This function also removes likely fill values above 100, which typically indicate invalid or missing data.
#' @param x A `SpatRaster` or numeric raster object representing MODIS LAI values in unscaled integer format (e.g., 0–100 or greater).
#' @return A scaled raster object (same class as input) with LAI values rescaled by a factor of 0.1 and invalid values set to `NA`.
#' @import terra
#' @noRd
scale_modis_lai <- function(x) {
  x[x > 100] <- NA  # remove likely fill values
  x * 0.1
}

#' @title Fit Sinusoidal Model to Time Series
#' @description Fits a harmonic regression model of the form y = a0 + a1*sin(ωt) + b1*cos(ωt) for a single pixel time series, given a precomputed design matrix cross-product.
#' @param y Numeric vector. Time series values (e.g., monthly LAI values) for a single pixel.
#' @param XtX Precomputed matrix cross-product (e.g., solve(t(X) %*% X) %*% t(X)) for sinusoidal regression.
#' @return Numeric vector of length 3 containing coefficients a0, a1, b1. Returns NA if all values in y are NA.
#' @noRd
fit_sinusoid <- function(y, XtX) {
  if (all(is.na(y))) return(rep(NA, 3))
  coef <- XtX %*% y  # [3 x 1]
  return(as.vector(coef))   # a0, a1, b1
}


#' @title Predict Seasonal Curve from Harmonic Coefficients
#' @description Uses harmonic regression coefficients (a0, a1, b1) to predict values across 12 months using a sinusoidal model.
#' @param a0 Numeric. Intercept term from harmonic regression.
#' @param a1 Numeric. Sine coefficient.
#' @param b1 Numeric. Cosine coefficient.
#' @return Numeric vector of length 12 representing predicted values for each month (January to December).
#' @details Requires that `omega <- 2 * pi / 12` and `month <- 1:12` are defined in the calling environment.
#' @noRd
predict_curve <- function(a0, a1, b1) {
  a0 + a1 * sin(omega * month) + b1 * cos(omega * month)
}

#' @title Model Intra-Annual Variation in Plant Area Index (PAI)
#' @description Estimates pixel-wise harmonic regression coefficients describing monthly PAI variation, using a 4D MODIS LAI array. Adds a 20% offset to account for PAI and fits a sinusoidal model to the mean annual cycle.
#' @param r.array 3D array (rows × cols × months × years) of scaled MODIS LAI values.
#' @param nyr Integer. Number of years in the time series.
#' @param nmons Integer. Total number of layers (e.g., months × years). Must equal `nyr * 12`.
#' @return A list with two elements: \code{coef_array}, a 3-layer array of regression coefficients (a0, a1, b1); and \code{pai_MODIS}, a 3D array of mean monthly PAI values across years.
#' @noRd
model_PAI <- function(r.array, nyr, nmons){
  nrow <- dim(r.array)[1]
  ncol <- dim(r.array)[2]
  nmon <- 12
  
  if (nmons != nyr*nmon) stop("Missing Months.")
  stk.array <- array(r.array, dim = c(nrow,ncol,nmon,nyr))
  
  #' Convert to PAI.
  maximums <- apply(stk.array, c(1,2),FUN = max)
  # calculate 20% of maximum LAI to account for PAI.
  pai_multi <- maximums * 0.2 # 20% of max lai
  pai_4D <- array(pai_multi, dim = dim(stk.array))
  PAI_array <- stk.array + pai_4D  
  
  mean_mon <- apply(PAI_array, c(1,2,3), mean, na.rm = TRUE)
  month <- 1:nmon
  omega <- 2 * pi / nmon  # frequency
  # Matrix for sinusoidal regression:
  X <- cbind(1, # Intercept
             sin(omega * month), # Sine term
             cos(omega * month)) # Cosine term
  
  XtX <- solve(t(X) %*% X) %*% t(X)  # Pre-compute (X'X)^-1 X' for fast OLS
  
  mean_mon_flat <- matrix(mean_mon, nrow = nrow * ncol, ncol = nmon)
  
  # Fit model to each pixel
  coefs <- t(apply(mean_mon_flat, 1, fit_sinusoid, XtX))
  
  # Reshape back to spatial
  coef_array <- array(coefs, dim = c(nrow, ncol, 3))
  
  return(list(coef_array = coef_array, pai_MODIS = mean_mon))
}

#' @title Scale Predicted Monthly PAI Using Observed Values
#' @description Scales a monthly PAI stack using the ratio of observed to expected PAI for a given input month, based on harmonic regression coefficients.
#' @param pai_temp SpatRaster. A single-layer raster of observed PAI for a specific height band and month.
#' @param coef.rf SpatRaster with three layers named "a0", "a1", and "b1", representing harmonic regression coefficients.
#' @param input_month Integer (1–12). Month corresponding to \code{pai_temp}.
#' @return SpatRaster. A 12-layer stack of monthly PAI values scaled by the observed-to-expected ratio.
#' @details Assumes that a global object \code{fine_predicted_stack} (12-month PAI predictions) is already defined in the environment.
#' @noRd
scale_pai <- function(pai_temp, coef.rf, input_month, predicted){
  a0_fine <- coef.rf$a0
  a1_fine <- coef.rf$a1
  b1_fine <- coef.rf$b1
  omega <- 2 * pi / 12
  
  expected <- a0_fine + a1_fine * sin(omega * input_month) + b1_fine * cos(omega * input_month)
  scaling_factor <- pai_temp / expected
  scaled_stack <- predicted * scaling_factor
  names(scaled_stack) <- month.name
  
  return(scaled_stack)
}


#' @title Identify Inflection Points in LAD Profiles
#' @description
#' Fits a GAM across interpolated LAD profiles from multiple LiDAR-derived vertical structures. Identifies heights where the LAD sharply declines, potentially indicating canopy base or ecotone transitions.
#' @param profile.list A named list of data frames, each with columns \code{z} (height) and \code{lad} (leaf area density). Null entries are filtered out.
#' @param plotInflec Logical. Whether to plot the smoothed LAD profile and detected inflection points. Default is \code{TRUE}.
#' @return A data frame with columns \code{z} (height) and \code{lad} (LAD value at inflection), representing significant LAD drop-off zones.
#' @importFrom mgcv gam
#' @importFrom tidyr pivot_longer
#' @importFrom stats approx
#' @noRd
.findInflec <- function(profile.list, plotInflec = TRUE) {
  stopifnot(requireNamespace("zoo", quietly = TRUE))
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  
  profile.list <- Filter(Negate(is.null), profile.list)
  
  if (length(profile.list) == 0) stop("No valid profiles in profile.list.")
  
  # Reference profile to align height bins
  lengths <- sapply(profile.list, function(df) length(df$z))
  z_ref <- profile.list[[which.max(lengths)]]$z
  
  # Interpolate LAD profiles to reference heights
  lad_matrix <- sapply(profile.list, function(df) {
    approx(x = df$z, y = df$lad, xout = z_ref, rule = 1, ties = mean)$y
  })
  
  lad_df <- data.frame(z = z_ref, lad_matrix)
  lad_df <- lad_df[complete.cases(lad_df), ]
  lad_long <- tidyr::pivot_longer(lad_df, cols = -z, names_to = "tile", values_to = "lad")
  
  # Fit GAM model
  gam_fit <- mgcv::gam(lad ~ s(z, k = 20), data = lad_long)
  
  # Predict over smooth height grid
  z_pred <- seq(min(lad_long$z), max(lad_long$z), by = 0.05)
  pred_df <- data.frame(z = z_pred)
  pred_df$lad <- as.numeric(predict(gam_fit, newdata = pred_df))
  
  # First derivative
  dz <- diff(z_pred)[1]
  pred_df$lad_d1 <- c(NA, diff(pred_df$lad) / dz)
  
  # Find local minima in dLAD/dz
  d1 <- pred_df$lad_d1
  drop_idx <- which(diff(sign(diff(d1))) == 2) + 1
  drop_idx_strong <- drop_idx[d1[drop_idx] < quantile(d1, 0.2, na.rm = TRUE)]
  
  drop_df <- data.frame(z = pred_df$z[drop_idx_strong],
                        lad = pred_df$lad[drop_idx_strong])
  
  if (plotInflec) {
    plot(pred_df$z, pred_df$lad, type = "l", col = "grey", lwd = 2,
         main = "LAD Profile with Inflection Zones",
         xlab = "Height (m)", ylab = "LAD")
    points(drop_df$z, drop_df$lad, col = "darkgreen", pch = 16)
  }
  
  return(drop_df)
}

#' Create SpatRaster object using a template
#' @import terra
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  terra::crs(r)<-terra::crs(tem)
  r
}

#' Used by albedo calculations - calculates wavelength specific spectral density
.Planck <- function(wavelength, temperature =  5504.85) {
  d <- wavelength * 1e-09
  h <- 6.6256e-34
  cc <- 299792458
  tt <- temperature + 273.15
  k <- 1.38054e-23
  b <- (2 * pi * h * cc^2) / (d^5 * (exp((h * cc) / (k * d * tt)) - 1))
  b
}

# fills NA values in a SpatRaster
.fillna<-function(ri,msk,zerotoNA=TRUE) {
  if (zerotoNA) ri[ri==0]<-NA
  mx<-min(dim(ri)[1:2])
  if (mx > 2) {
    r2<-.rastna(ri,2,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  if (mx > 4) {
    r2<-.rastna(ri,4,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  if (mx > 8) {
    r2<-.rastna(ri,8,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  if (mx > 16) {
    r2<-.rastna(ri,16,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  mm<-mean(as.vector(ri),na.rm=TRUE)
  ri[is.na(ri)]<-mm
  ri<-mask(ri,msk)
  return(ri)
}

#' Sorts out NAs in a raster dataset due to coarsening
.rastna<-function(r,af,msk) {
  r<-mask(r,msk)
  rc<-resample(aggregate(r,af,"mean",na.rm=TRUE),r)
  m<-.is(r)
  mc<-.is(rc)
  s<-which(is.na(m))
  m[s]<-mc[s]
  ro<-.rast(m,r)
  return(ro)
}

.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  if (class(r)[1] != "matrix") {
    if (dim(r)[3] > 1) {
      y<-as.array(r)
    } else y<-as.matrix(r,wide=TRUE)
  } else y<-r
  y
}

.extract_clima <- function(nc, r, resampleout = FALSE) {
  # get time
  nc_file <- ncdf4::nc_open(nc)
  tme <- as.POSIXlt(ncdf4::ncvar_get(nc_file, "valid_time"),
                    origin="1970-01-01 00:00", tz = "UTC")
  ncdf4::nc_close(nc_file)
  # extract variables
  varn <- c("t2m", "d2m", "sp", "u10" , "v10",  "tp", "avg_sdlwrf", "fdir", "ssrd", "lsm")
  rlst <- list()
  for (i in 1:9) rlst[[i]]<-rast(nc, subds=varn[i])
  rlst[[10]]<-rast(nc, subds=varn[10])[[1]]
  lsm <-  rlst[[10]]
  tc <- rlst[[1]]-273.15
  # coastal correction
  if (any(terra::values(lsm) < 1)) {
    # Calculate daily average
    # Indices to associate each layer with its yday
    ind <- rep(1:(dim(tc)[3]/24), each = 24)
    # Average across days
    tmean <- terra::tapp(tc, ind, fun = mean, na.rm = T)
    # Repeat the stack 24 times to expand back out to original timeseries
    tmean <- rep(tmean, 24)
    # Sort according to names so that the stack is now in correct order: each
    # daily mean, repeated 24 times
    # your pasted command properly sorts the names, X1 to X365 (or X366)
    tmean <- tmean[[paste0("X", sort(rep(seq(1:(dim(tc)[3]/24)), 24)))]]
    m <- (1 - lsm) * 1.285 + 1
    tdif <- (tc - tmean) * m
    tc<- tmean + tdif
  }
  rlst[[1]]<-tc
  # crop variables
  e<-ext(r)
  rr<-rast(e)
  terra::crs(rr)<-terra::crs(r)
  rll<-project(rr,"EPSG:4326")
  ell<-ext(rll)
  for (i in 1:10) rlst[[i]]<-crop(rlst[[i]],ell,snap='out')
  if (resampleout) {
    for (i in 1:10) {
      if (terra::crs(r) != terra::crs(rll)) {
        rlst[[i]]<-project(rlst[[i]],r)
      } else {
        rlst[[i]]<-resample(rlst[[i]],r)
      }
      rlst[[i]]<-mask(rlst[[i]],r)
    }
  }
  rlsto<-list()
  rte<-rlst[[1]]
  rte<-rte[[1]]
  # convert variables
  rlsto[[1]]<-rlst[[1]]
  rh <- .rast((satvapCpp(.is(rlst[[2]])-273.15) /  satvapCpp(.is(rlst[[1]]))) * 100, rte)
  rh[rh > 100]<-100
  rlsto[[2]]<-rh
  rlsto[[3]] <- rlst[[3]]/1000
  rlsto[[4]] <- rlst[[9]]/3600
  dni <- rlst[[8]]/3600
  ll<-.latslonsfromr(rte)
  si <- .rast(solarindexarray(tme$year+1900, tme$mon+1, tme$mday, tme$hour, ll$lats, ll$lons), rte)
  rlsto[[5]]  <- rlsto[[4]]  - (si * dni)
  rlsto[[6]]  <- rlst[[7]]
  rlsto[[7]]<-sqrt(rlst[[4]]^2+rlst[[5]]^2)*0.7477849 # Wind speed (m/s)
  rlsto[[8]]<-(atan2(rlst[[4]],rlst[[5]])*180/pi+180)%%360
  rlsto[[9]]<- rlst[[6]] * 1000
  names(rlsto)<-c("temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
  for (i in 1:9)  time(rlsto[[i]])<-as.POSIXct(tme)
  return(rlsto)
}

.latslonsfromr <- function(r) {
  lats<-.latsfromr(r)
  lons<-.lonsfromr(r)
  xy<-data.frame(x=as.vector(lons),y=as.vector(lats))
  xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = terra::crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2],
                   long = sf::st_coordinates(ll)[,1])
  lons<-array(ll$long,dim=dim(lons))
  lats<-array(ll$lat,dim=dim(lats))
  return(list(lats=lats,lons=lons))
}

#' Latitudes from SpatRaster object
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' Longitudes from SpatRaster object
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}


.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  terra::crs(r)<-terra::crs(tem)
  r
}

#' sharpens a coarse resolution raster using a fine-resolution categorical raster
.sharpen <- function(coarse, fine, msk = TRUE) {
  if (msk) {
    mskc <- resample(fine, coarse)
    coarse <- mask(coarse, mskc)
  }
  # Produce table of categorical fine scale variable
  ufine <- unique(as.vector(fine))
  ufine <- ufine[is.na(ufine) == FALSE]
  # Calculate the mean value of coarse for each unique value of fine
  mcoarse <- 0
  for (i in 1:length(ufine)) {
    # Calculate fraction of each fine type in each coarse grid cell
    ftype <- fine
    ftype[ftype != ufine[i]] <-0
    ftype[ftype == ufine[i]] <-1
    cfrac <- resample(ftype, coarse)
    # multiply coarse by fraction and sum
    mucoarse <- cfrac * coarse
    fsum <- sum(as.vector(mucoarse), na.rm=TRUE)
    # Sum the fractions across cells
    csum <- sum(as.vector(cfrac), na.rm=TRUE)
    # Calculate the mean coarse for each fine type
    mcoarse[i] <- fsum / csum
  }
  mx <- (max(mcoarse) / mean(mcoarse)) * 3
  # Calculate a fine-scale dataset of expected values of coarse
  finem <- as.matrix(fine, wide = TRUE)
  evals <- finem * 0
  for (i in 1:length(ufine)) {
    s <- which(finem == ufine[i])
    evals[s] <- mcoarse[i]
  }
  evals <- .rast(evals, fine)
  # Calculate multiplier by dividing actual by expected (coarse scale)
  mu <- coarse / resample(evals, coarse)
  mskc <- mskc * 0 +1
  mskc[is.na(mskc)] <- 1
  mu <- .fillna(mu, mskc)
  mu <- resample(mu, fine)
  mu[mu > mx] <- mx
  sharpened <- mu * evals
  return(sharpened)
}



.modalReplace <- function(r) {
  vals <- values(r)
  valid_vals <- vals[!is.na(vals) & vals != 0]
  
  if (length(valid_vals) == 0) return(r)
  
  # Calculate mode manually
  tbl <- table(valid_vals)
  mval <- as.numeric(names(tbl)[which.max(tbl)])
  
  r[r == 0] <- mval
  return(r)
}


.clusterInfl <- function(cp){
  c.heights <- .findInflec(cp)$z
  c.max <- round(max(sapply(cp, function(df) max(df$z, na.rm = TRUE))),1)
  return(list(gheights = c.heights, gmax = c.max))
}

#' Campbell leaf angle distribution approximation:
.Kcanopy <- function(x) {
  if (!is.numeric(x)) stop("x must be numeric.")
  if (any(x < 0.1 | x > 10)) warning("x should be between 0.1 and 10 for approximation to hold.")
  K <- x / (x + 1.774 * (x + 1.182)^(-0.733))
  return(K)
}


# Function to fill a single layer using CHM
.fillNas <- function(layer, chm.low) {
  df <- as.data.frame(c(layer, chm.low), na.rm = TRUE)
  colnames(df) <- c("micro", "chm")
  df <- df[is.finite(df$micro) & is.finite(df$chm), ]
  
  # Guard against empty or sparse data
  if (nrow(df) < 10) return(layer)
  
  # Try-catch in case lm() still fails (e.g., due to zero variance)
  tryCatch({
    mod <- lm(micro ~ chm, data = df)
    predicted <- terra::predict(chm.low, mod)
    filled <- layer
    na.idx <- is.na(values(layer))
    values(filled)[na.idx] <- values(predicted)[na.idx]
    return(filled)
  }, error = function(e) {
    warning("Skipping layer due to lm() failure: ", conditionMessage(e))
    return(layer)
  })
}


.precipcorrect <- function(hist_obs, hist_mod, fut_mod, rangelims = 1.05) {
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
  #fut_mod2<-mask(fut_mod,hist_obs) # this line as causing issues
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
