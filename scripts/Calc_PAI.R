#' Manuscript: Vertical microclimate diversity shapes thermal refugia in Congolese rainforests.
#' Date: August 2026
#' Author: Brittany Trew (*brittany.trew@fas.harvard.edu*)
#' Description: XXXX
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R"))

#' Canopy Height Model (template)
chm_mosaic <- rast(paste0(head.path,"chm.tif"))

#' Estimate PAI (cumulative LAD) from Leaf Area Density (m2m3)
all.pad <- list.files(paste0(head.path,"pad"), pattern = ".tif", full.names = T)
pad.r <- lapply(all.pad, rast)

#' Tiles need to be padded with NAs at higher levels because tiles have different heights.
#' Maximum number of vertical layers across tiles
max_nlyr <- max(sapply(pad.r, terra::nlyr))

pad.list <- lapply(pad.r, function(r) {
  n_missing <- max_nlyr - terra::nlyr(r)
  if (n_missing > 0) {
    # Empty layer with exactly the same spatial geometry
    empty <- r[[1]]
    empty[] <- NA
    # Add missing layers to TOP of the vertical stack
    for (j in seq_len(n_missing)) {
      r <- c(r, empty)
    }
  }
  r
})

#' Mosaic the tiles together.
pad.r <- do.call(mosaic,pad.list)
pad.r <- terra::resample(pad.r, chm_mosaic)

start_heights <- seq(
  from = 0,
  by = dzd,
  length.out = nlyr(pad.r)
)

names(pad.r) <- sprintf(
  "PAD_%.1f_%.1fm",
  start_heights,
  start_heights + dzd
)


#' Distinguish no-data pixels from above canopy pixels
valid <- sum(!is.na(pad.r)) > 0
#' Then the above canopy pixels contribute 0
pad0 <- terra::ifel(
  is.na(pad.r),
  0,
  pad.r
)

#' Cumulative PAD (PAI)
ii <- nlyr(pad0):1
PAI.stk <- cumsum(
  pad0[[ii]]
)[[ii]] * dzd

#' Remove no data pixels
PAI.stk <- terra::mask(
  PAI.stk,
  valid,
  maskvalues = 0
)

names(PAI.stk) <- sprintf(
  "PAI_%.1fm_to_canopy",
  start_heights
)

pai.path <- paste0(head.path,"pai/")
dir.create(pai.path, showWarnings = F)
writeRaster(PAI.stk, 
            paste0(pai.path,"pai_noSeasonalAdjust.tif"),
            overwrite = T)


#' Perform a seasonal adjustment using coarse MODIS data.
MODISpath <- paste0(in.data,"MODIS_LAI/",site.location,"/",sample.name,"/")

pai_list <- .MODISAdjust(MODISpath, PAI.stk, dzd)

.MODISAdjust <- function(MODISpath, PAI, dzd){
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
  
  #' **Model monthly PAI in LiDAR**
  pai_list <- list()
  for(i in 1:nlyr(pai_r)){
    pai_temp <- pai_r[[i]]
    scaled_stack <- scale_pai(pai_temp, coef.rf, input_month, fine_predicted_stack)
    pai_list[[i]] <- scaled_stack
    names(pai_list)[i] <- names(pai_temp)
  }
  
  return(pai_list)
}









for(i in 1:length(pai_list)){
  hm <- names(pai_list)[[i]]
  print(hm)
  r <- pai_list[[i]]
  plot(r[[1]])
  writeRaster(pai_list[[i]], paste0(pai.path,hm,".tif"), overwrite = T)
}




