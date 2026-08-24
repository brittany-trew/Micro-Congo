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
})

#' Check all tils have the same number of layers.
sapply(pad.list, terra::nlyr)

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


#' Perform a seasonal adjustment.
MODISpath <- paste0(in.data,"MODIS_LAI/",site.location,"/",sample.name,"/")

pai_list <- .MODISAdjust(MODISpath, PAI.stk, dzd)

for(i in 1:length(pai_list)){
  hm <- names(pai_list)[[i]]
  print(hm)
  r <- pai_list[[i]]
  plot(r[[1]])
  writeRaster(pai_list[[i]], paste0(pai.path,hm,".tif"), overwrite = T)
}




