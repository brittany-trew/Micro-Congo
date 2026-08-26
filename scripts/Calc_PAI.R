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
plot(PAI.stk[[1]])

pai.path <- paste0(head.path,"pai/")
dir.create(pai.path, showWarnings = F)
writeRaster(PAI.stk, 
            paste0(pai.path,"pai_noSeasonalAdjust.tif"),
            overwrite = T)


#' Perform a seasonal adjustment using coarse MODIS data.
MODISpath <- paste0(mclidar.in,"MODIS_LAI/",site.location,"/",sample.name,"/")

pai_list <- .MODISAdjust(MODISpath, PAI.stk, input_month)

for(i in 1:length(pai_list)){
  hm <- names(pai_list)[i]
  print(hm)
  
  r <- pai_list[[i]]
  # Replace zeros with NA
  r[r == 0] <- NA
  
  # Check whether the entire raster stack is NA
  n_valid <- terra::global(
    !is.na(r),
    "sum",
    na.rm = TRUE
  )[, 1]
  
  # Skip if there are no valid values
  if (sum(n_valid, na.rm = TRUE) == 0) {
    message("Skipping ", hm, " — all values are NA")
    next
  }
  
  plot(r[[1]])
  
  terra::writeRaster(
    r,
    file.path(pai.path, paste0(hm, ".tif")),
    overwrite = TRUE
  )
}




