#' @title Load and Install Required R Packages
#' @description
#' Loads a vector of R package names. If a package is not already installed, the function installs it from CRAN and then loads it. Useful for reproducible workflows or scripts with many dependencies.
#' @param packages Character vector. Names of the packages to load.
#' @return No return value. Side effect: packages are loaded into the session (and installed if missing).
#' @export
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
      message(paste("Installing missing package:", pkg))
      install.packages(pkg, quiet = TRUE)
    }
    suppressPackageStartupMessages(
      suppressMessages(
        suppressWarnings(
          library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        )
      )
    )
  }
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


#' Process raw PAD output:
processPAD <- function(head.path, dzd, chm){
  all.pad <- list.files(paste0(head.path,"pad"), pattern = ".tif",full.names = TRUE)
  pad.list <- lapply(all.pad, rast)
  pad.r <- do.call(mosaic, pad.list)
  chm <- resample(chm, pad.r)
  
  # ---- heights ----
  n  <- nlyr(pad.r)
  z0 <- (0:(n - 1)) * dzd # height at the *bottom* of each layer (m), starting at 0
  keep <- which(z0 >= dzd) # keep layers starting at dzd and above
  pad.aboveG <- pad.r[[keep]] # Exclude the ground layer as above
  z0G <- z0[keep] # the corresponding height vector for pad.aboveG (same length as nlyr(pad.aboveG))
  
  chm_max <- as.numeric(global(chm, "max", na.rm = TRUE)[1, 1])  # maximum canopy height in the site (m)
  max_start <- floor(chm_max / dzd) * dzd # snap that max height down to the dzd grid (e.g., 23.87 -> 23.75)
  start_heights <- z0G[z0G <= max_start]  # allowable start heights (h0) up to the tallest canopy, on the dzd grid
  pad.aboveG <- pad.aboveG[[seq_along(start_heights)]]  # truncate PAD layers so they match the number of start heights (keeps indexing consistent)
  # Replace NA in PAD with 0 wherever CHM has data
  pad.aboveG[[1]] <- ifel(
    is.na(pad.aboveG[[1]]) & !is.na(chm),
    0,
    pad.aboveG[[1]]
  )
  return(list(pad.aboveG, start_heights, z0G))
}


# Extract height (metres) from parent folder name
# e.g. ".../dailyTemps/5.0m/DailyMax_2009.tif" -> 5
get_height <- function(f) {
  as.numeric(gsub("m$", "", basename(dirname(f))))
}

# Load a single height-band raster, resample to CHM grid, then mask pixels
# where canopy height is below this height band (i.e. above canopy)
load_masked <- function(f, chm) {
  hh <- get_height(f)
  r  <- resample(rast(f), chm)
  ifel(chm <= hh, NA, r)
}

collapse_vertical <- function(files, chm) {
  r.lst  <- lapply(files, load_masked, chm = chm)
  r.stk  <- rast(r.lst)
  r.stk[r.stk < 0] <- NA
  
  n_days <- nlyr(r.lst[[1]])
  index  <- rep(1:n_days, times = length(r.lst))
  colmax <- tapp(r.stk, index, fun = max, na.rm = TRUE)
  colmin <- tapp(r.stk, index, fun = min, na.rm = TRUE)
  colsd <- tapp(r.stk, index, fun = sd, na.rm = TRUE)
  return(list("colmax"=colmax, 
              "colmin" = colmin, 
              "colsd" = colsd))
  gc()
}
