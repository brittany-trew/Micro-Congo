#' Manuscript: Vertical microclimate diversity shapes thermal refugia in Congolese rainforests.
#' Date: August 2026
#' Author: Brittany Trew (*brittany.trew@fas.harvard.edu*)
#' Description: Process land cover and ground metrics.
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

#' [Create fine-scale habitat layers]
file1 <- paste0(scripts.path, "land_cover/process_landcover_", sample.name, ".R")
file2 <- paste0(scripts.path, "land_cover/process_landcover_", site.location, ".R")

if (file.exists(file1)) {
  source(file1)
} else if (file.exists(file2)) {
  source(file2)
} else {
  warning("Neither landcover script exists.")
}

#' [Process soil-type.]
lc <- rast(paste0(head.path,"/landcover.tif"))
lc.wgs <- project(lc, "EPSG:4326")

if(!exists(paste0(head.path,"soiltype.tif"))){
  soil.temp <- paste0(head.path,"/soilgrid_temp/")
  dir.create(soil.temp, showWarnings = F)
  soildata <- soildata_download(r = lc.wgs, pathdir = soil.temp, deletefiles = TRUE)
  soildata <- project(soildata, terra::crs(lc), method = "near")
  soildd <- soildata_downscale(soildata, lc)
  soiltype <- soildata_gettype(soildd)
  plot(soiltype)
  writeRaster(soiltype, paste0(head.path,"soiltype.tif"), overwrite = T)
} else {soiltype <- rast(paste0(head.path,"soiltype.tif"))}


#' [Download Albedo and calculate Ground Reflectance.]
#' Uses MODIS albedo data
username = "bt302@exeter.ac.uk"
password = "Buddy250604"
r <- lc
pathout <- paste0(head.path,"/albedo/")
dir.create(pathout, showWarnings = F)

is_empty <- length(list.files(pathout, all.files = TRUE, no.. = TRUE)) == 0

if(is_empty){
  tme <- as.POSIXlt(0:(8760 - 1) * 3600, origin = "2023-01-01 00:00", tz = "UTC")
  credentials <- data.frame(
    username = username,
    password = password,
    stringsAsFactors = FALSE
  )
  albedo_download(r, tme, pathout, credentials)
}

albmodis <- albedo_process(r, pathout)
plot(albmodis)
writeRaster(albmodis, paste0(head.path,"/albedo_modis.tif"), overwrite = T)


albmodis <- rast(paste0(head.path,"/albedo_modis.tif"))
#' Needs sat imagery to downscale. 
sat.CIR <- list.files(paste0(head.path,"/sentinel-2/"), 
                      pattern = paste0("_CIR_.*masked.tif$"),full.names = T)
sat.RGB <- list.files(paste0(head.path,"/sentinel-2/"), 
                      pattern = paste0("_RGB_.*masked.tif$"),full.names = T)

CIR.r <- rast(sat.CIR)
RGB.r <- rast(sat.RGB)

albphoto <- albedo_fromaerial(RGB.r, CIR.r)
albadjusted <- albedo_adjust(albphoto, albmodis)
plot(albadjusted)

pai_ground <- rast(paste0(head.path,"pai/PAI_0.0m_to_canopy.tif"))
PAI_array <- as.array(pai_ground)
max_pai <- apply(PAI_array, c(1,2), FUN = max, na.rm = TRUE)
lai_addition <- max_pai * 0.2
lai_3D <- array(rep(lai_addition, times = 12), dim = dim(PAI_array))
LAI_array <- PAI_array - lai_3D
lai.r <- rast(LAI_array, crs = terra::crs(pai_ground), ext = ext(pai_ground))
lai.mean <- mean(lai.r)
plot(lai.mean[[1]])

albadjusted <- project(albadjusted, terra::crs(lai.mean))
albadjusted <- resample(albadjusted, lai.mean)
x <- lai.mean
values(x) <- 1
ref <- reflectance_calc(alb = albadjusted, lai = lai.mean, x)

dir.create(paste0(head.path,"/reflectance/"), showWarnings = F)
writeRaster(ref$gref, paste0(head.path,"/reflectance/groundR.tif"), overwrite = T)
writeRaster(ref$lref, paste0(head.path,"/reflectance/leafR.tif"), overwrite = T)

