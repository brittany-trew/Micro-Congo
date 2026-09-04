scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R"))

#' Tiling and climate download.
era5_temp <- paste0(out.data,site.location,"/","era5_temp/")
dir.create(era5_temp, showWarnings = F)

era5_out <- paste0(head.path,"/era5/") # processed ERA5 output folder
dir.create(era5_out, showWarnings = F)

r1 <- rast(paste0(head.path,"/landcover.tif"))
r <- terra::aggregate(r1,100, fun = "modal")

library(httr)
set_config(timeout(60000))
yr.seq <- seq(2024,2004,-1)
for(i in 1:length(yr.seq)){
  
  yr <- yr.seq[i]
  print(yr)
  
  era5_yr <- paste0(era5_temp,yr,"/")
  dir.create(era5_yr, showWarnings = F)
  
  start_time <- ymd_h(paste0(yr, "-01-01 00"))
  end_time   <- ymd_h(paste0(yr, "-12-31 23"))
  tme <- seq(from = start_time, to = end_time, by = "1 hour")
  
  
  file_prefix <- paste0("era5_",yr)
  req <- era5_download(r, tme, credentials = creds, file_prefix, pathout = era5_yr)
  req[[1]]$target <- sub("\\.zip$", ".nc", req[[1]]$target)


  #' Remove surplus folders and files (wish this junk didn't appear!).
  all.files <- list.files(era5_yr, full.names = TRUE)
  nc_files <- list.files(era5_yr, pattern = "\\.nc$", full.names = TRUE)
  to_remove <- setdiff(all.files, nc_files)
  is_dir <- file.info(to_remove)$isdir
  dirs_to_remove <- to_remove[is_dir]
  files_to_remove <- to_remove[!is_dir]
  if (length(files_to_remove) > 0) file.remove(files_to_remove)
  if (length(dirs_to_remove) > 0) unlink(dirs_to_remove, recursive = TRUE)

}

tile.path <- paste0(head.path,"tiles/")
dir.create(tile.path, showWarnings = F)

r <- crop(r, r1, snap="in")

expand_to_cover <- function(template, target){
  r <- template
  tx <- ext(target)
  
  # adjust x
  if (xmax(r) < xmax(tx)) xmax(r) <- xmax(r) + res(r)[1]
  if (xmin(r) > xmin(tx)) xmin(r) <- xmin(r) - res(r)[1]
  
  # adjust y
  if (ymin(r) > ymin(tx)) ymin(r) <- ymin(r) - res(r)[2]
  if (ymax(r) < ymax(tx)) ymax(r) <- ymax(r) + res(r)[2]
  
  extend(template, r)
}

r <- expand_to_cover(r, r1)
r[is.nan(r)] <- NA
r <- trim(r)
tf <- paste0(tile.path,"tile_1.tif")
if(!file.exists(tf)) makeTiles(x = r, 
                               y = 1, 
                               filename = paste0(tile.path,"tile_.tif"), 
                               buffer = 1,
                               na.rm = T)
all.tiles <- list.files(tile.path, full.names = T)

yr.seq <- seq(2024,2004,-1)
for(i in 1:length(yr.seq)){
  
  yr <- yr.seq[i]
  print(yr)
  
  era5_yr <- paste0(era5_temp,yr,"/")
  dir.create(era5_yr, showWarnings = F)
  
  era5_out_yr <- paste0(era5_out,yr,"/")
  dir.create(era5_out_yr, showWarnings = F)
  print(era5_out_yr)
  start_time <- ymd_h(paste0(yr, "-01-01 00"))
  end_time   <- ymd_h(paste0(yr, "-12-31 23"))
  tme <- seq(from = start_time, to = end_time, by = "1 hour")
  
  for(t in 1:length(all.tiles)){
    
    tile.file <- all.tiles[t]
    num <- as.numeric(gsub(".*tile_(\\d+)\\.tif$", "\\1", basename(tile.file)))
    num <- sprintf("%02d", num)
    print(paste0("Tile ", num, "."))
    
    r.tile <- rast(tile.file)
    #' Keep req as NULL but only have single year of files available in folder.
    era5climdata <- era5_process(req = NA, pathin = era5_yr, r.tile, tme, out = 'grid')
    
    saveRDS(era5climdata, paste0(era5_out_yr,"ERA5_tile_",num,".RDS"))
  }
}
