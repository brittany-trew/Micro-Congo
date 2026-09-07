scripts.path <- "scripts/"
source(paste0(scripts.path,"imbalanga/microParameters.R")) # loads worker functions

# -------------------------------
# Cluster batch command
# -------------------------------
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
d <- array_id + offset

# -------------------------------
# Functions
# -------------------------------

merge.tiles.day <- function(files, tile.path, chm){
  #' Available tile rasters
  tile.files <- list.files(
    tile.path,
    pattern = "tile_\\d+\\.tif$",
    full.names = TRUE
  )
  
  tile.ids <- as.integer(
    sub(".*tile_(\\d+)\\.tif$", "\\1", basename(tile.files))
  )
  
  m.list <- list()
  max.list <- list()
  min.list <- list()
  
  for(i in 1:length(files)){
    
    micro <- readRDS(files[[i]])
    
    #' Extract tile number from model output
    num <- as.integer(
      sub(".*tile_(\\d+)\\.RDS$", "\\1", basename(files[[i]]))
    )
    
    #' Match to corresponding tile raster
    tile.file <- tile.files[match(num, tile.ids)]
    
    tile <- rast(tile.file)
    template <- crop(chm, tile)
    
    mmean <- rast(micro$dailyMean, crs = crs(template), ext = ext(template))
    m.list[[i]] <- mmean
    
    mmax <- rast(micro$dailyMax, crs = crs(template), ext = ext(template))
    max.list[[i]] <- mmax
    
    mmn <- rast(micro$dailyMin, crs = crs(template), ext = ext(template))
    min.list[[i]] <- mmn
  }
  
  gc()
  m.list <- fix.gridcells(list.r = m.list, template = chm)
  valid.list <- Filter(Negate(is.null), m.list)
  mean.mos <- do.call(mosaic, c(valid.list, fun = mean))
  
  gc()
  max.list <- fix.gridcells(list.r = max.list, template = chm)
  valid.list <- Filter(Negate(is.null), max.list)
  max.mos <- do.call(mosaic, c(valid.list, fun = max))
  
  gc()
  min.list <- fix.gridcells(list.r = min.list, template = chm)
  valid.list <- Filter(Negate(is.null), min.list)
  min.mos <- do.call(mosaic, c(valid.list, fun = min))
  
  gc()
  
  return(list(
    mean = mean.mos,
    max = max.mos,
    min = min.mos
  ))
}


fix.gridcells <- function(list.r, template) {
  
  lapply(list.r, function(r) {
    
    #' Ensure CRS matches
    if(!same.crs(r, template)){
      r <- project(r, template)
    }
    
    #' Align grid
    r <- resample(r, template)
    
    return(r)
  })
}


# -------------------------------
# Set up
# -------------------------------

#' PAI layers are every 0.5 m
pai.layers <- list.files(
  paste0(head.path,"/pai"),
  pattern = "_to_canopy\\.tif$"
)

n.h <- length(pai.layers)

#' Maximum available height
max.height <- n.h * 0.5

#' Model heights: 0.2 m then every 1 m
model.heights <- c(0.2, seq(1, max.height, by = 1))
n_heights <- length(model.heights)

#' One Slurm task = one height
hm <- model.heights[d]
height <- sprintf("%02.1fm", hm)

print(paste0(
  "Running all years for height ", height,
  " (height ", d, " of ", n_heights, ")."
))


#' Paths
chm <- rast(paste0(head.path,"/chm.tif"))

in.path <- paste0(head.path,"micro_models/")
out.path <- paste0(head.path,"dailyTemps/")
tile.path <- paste0(head.path,"tiles/")

dir.create(out.path, showWarnings = FALSE)

#' All years
yr.seq <- seq(2004, 2024, 1)


# -------------------------------
# Process all years for this height
# -------------------------------

for(a in 1:length(yr.seq)){
  
  yr <- yr.seq[a]
  
  print(paste0(
    "Processing ", height,
    " for year ", yr
  ))
  
  model.path <- paste0(
    in.path,
    height,
    "/",
    yr,
    "/"
  )
  
  #' Find whatever daily tile files are available
  daily <- list.files(
    model.path,
    pattern = "microDaily_.*\\.RDS$",
    full.names = TRUE
  )
  
  #' Nothing exists for this year, so move on
  if(length(daily) == 0){
    print(paste0("No files found for ", height, " ", yr, ". Skipping."))
    next()
  }
  
  print(paste0(
    "Merging ", length(daily),
    " tiles for ", yr
  ))
  
  #' Merge available tiles
  daily.r <- merge.tiles.day(
    files = daily,
    tile.path = tile.path,
    chm = chm
  )

  #' Set temperatures to NA where model height exceeds the canopy
  below_height_mask <- chm < hm
  daily.r <- lapply(daily.r,
                    function(r) {
                      mask(r,
                           below_height_mask,
                           maskvalues = 1,
                           updatevalue = NA)})
  #' Output folder
  outpath <- paste0(out.path,height,"/")
  dir.create(outpath, showWarnings = FALSE)
  
  #' Save daily statistics
  writeRaster(
    daily.r$mean,
    paste0(outpath,"DailyMean_",yr,".tif"),
    overwrite = TRUE
  )
  
  writeRaster(
    daily.r$max,
    paste0(outpath,"DailyMax_",yr,".tif"),
    overwrite = TRUE
  )
  
  writeRaster(
    daily.r$min,
    paste0(outpath,"DailyMin_",yr,".tif"),
    overwrite = TRUE
  )
  
  rm(daily.r)
  gc()
}


print(paste0(
  "Finished all years for height ",
  height,
  "."
))