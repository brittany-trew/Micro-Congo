scripts.path <- "scripts/"
source(paste0(scripts.path,"imbalanga/parameters.R")) # loads worker functions

# -------------------------------
# Cluster batch command
# -------------------------------
# Slurm array ID
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
# Optional offset (default 0 if not set)
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
# Final job index
d <- array_id + offset

# Heights (0.2, then 1 m intervals
pai.layers <- list.files(
  paste0(head.path,"/pai"), 
  pattern = "to_canopy.tif")
n.h <- length(pai.layers)

# Maximum available PAI height
max_pai_height <- n.h * 0.5

# Model at 0.2 m, then every 1 m
model.heights <- c(
  0.2,
  seq(1, max_pai_height, by = 1)
)

n_heights <- length(model.heights)

# Paths
tile.path <- paste0(head.path, "tiles/")
era5_out  <- paste0(head.path, "/era5/")

# Ensure consistent sorting
all.tiles <- list.files(tile.path, pattern = "tile_\\d+\\.tif$", full.names = TRUE)
yr.seq <- seq(2024, 2004, -1)
#yr.seq <- c(2015,2007)

n_tiles  <- length(all.tiles)
n_years  <- length(yr.seq)

# -------------------------------
# Linear index: Tile × Year × Height
# -------------------------------
d0 <- d - 1
tile_index   <- (d0 %% n_tiles) + 1
year_index   <- ((d0 %/% n_tiles) %% n_years) + 1
height_index <- (d0 %/% (n_tiles * n_years)) + 1

tile   <- all.tiles[tile_index]
yr     <- yr.seq[year_index]
hm     <- model.heights[height_index]
height_label <- sprintf("%02.1fm", hm)
print(height_label)

# Extract tile number for loading ERA5 data
num <- sprintf("%02d", as.numeric(sub(".*tile_(\\d+)\\.tif$", "\\1", basename(tile))))
print(paste0("Running tile ", num,
             " (index ", tile_index, ") for year ", yr,
             " (index ", year_index, ") at height ", height_label,
             " (index ", height_index, ")."))

# -------------------------------
# Load input rasters
# -------------------------------
lc      <- rast(paste0(head.path, "/landcover.tif"))
dem     <- rast(paste0(head.path, "/dem.tif"))
chm     <- rast(paste0(head.path, "/chm.tif"))

# Height-dependent PAI
if(hm == 0.2){
  pai <- rast(paste0(head.path, "/pai/PAI_0.5m_to_canopy.tif"))
}else{
  pai     <- rast(paste0(head.path, "/pai/PAI_", height_label, "_to_canopy.tif"))
}
pai_0m  <- rast(paste0(head.path, "/pai/PAI_0.0m_to_canopy.tif"))

refl.g  <- rast(paste0(head.path, "/reflectance/groundR.tif"))
refl.c  <- rast(paste0(head.path, "/reflectance/leafR.tif"))
soiltype<- rast(paste0(head.path, "/soiltype.tif"))

# Geometry checks
if(compareGeom(pai,lc) == FALSE) stop("Geometry of PAI rasters do not match land use")
if(compareGeom(soiltype,lc)== FALSE) stop("Geometry of soil-type raster does not match land use")
if(compareGeom(refl.g,lc)== FALSE) stop("Geometry of ground reflectance raster does not match land use")
if(compareGeom(refl.c,lc)== FALSE) stop("Geometry of canopy reflectance raster does not match land use")
if(compareGeom(dem,lc)== FALSE) stop("Geometry of DEM raster does not match land use")
if(compareGeom(chm,lc)== FALSE) stop("Geometry of canopy height raster does not match land use")

all.params <- list(
  pai = pai,
  pai_0m = pai_0m,
  ground.r = refl.g,
  canopy.r = refl.c,
  soil.type = soiltype,
  canopy.height = chm,
  dem = dem,
  landuse = lc
)

soil.m <- c(soiltype = all.params$soil.type, groundr = all.params$ground.r)
soil.m$soiltype <- .modalReplace(soil.m$soiltype)
all.params$landuse <- .modalReplace(all.params$landuse)

# -------------------------------
# Crop to tile
# -------------------------------
tile.temp <- terra::rast(tile)
tile.params <- lapply(all.params, function(x) crop(x, tile.temp))
soil.n <- lapply(soil.m, function(x) crop(x, tile.temp))

# -------------------------------
# Prep veg parameters
# -------------------------------
veggie <- microclimf::vegpfromhab(
  habitats = tile.params$landuse, 
  hgts = tile.params$canopy.height,
  pai = .is(tile.params$pai_0m),
  clump0 = TRUE
)

# -------------------------------
# Time object
# -------------------------------
start_time <- lubridate::ymd_h(paste0(yr, "-01-01 00"))
end_time   <- lubridate::ymd_h(paste0(yr, "-12-31 23"))
tme <- seq(from = start_time, to = end_time, by = "1 hour")

# -------------------------------
# Load ERA5 reanalysis tile
# -------------------------------
clim.file <- paste0(era5_out, yr, "/", "ERA5_tile_", num, ".RDS")
print(clim.file)
climdata <- readRDS(clim.file)

#' Add floor to radiation values (can't be < 0)
swr <- unwrap(climdata$swdown)
swr[swr<0] <- 0
climdata$swdown <- wrap(swr)

dfr <- unwrap(climdata$difrad)
dfr[dfr<0] <- 0
climdata$difrad <- wrap(dfr)

# -------------------------------
# Run point and array models
# -------------------------------
mm <- microclimf:::runpointmodela(
  climarrayr = climdata, tme = tme,
  reqhgt = hm, dtm = tile.params$dem, vegp = veggie,
  soilc  = soil.n, zref = 2, windhgt = 2
)

print("point model finished.")

dtmc <- resample(tile.params$dem, unwrap(climdata$temp)[[1]])
micromodel <- microclimf:::runmicro(
  micropoint = mm, reqhgt = hm, vegp = veggie, soilc  = soil.n,
  dtm = tile.params$dem, dtmc = dtmc, pai_a = .is(tile.params$pai),
  altcorrect = 2
)
rm(mm)
print("array model finished.")

# -------------------------------
# Convert to Daily and Monthly Statistics
# -------------------------------
dailyTz   <- compute.daily(tme[-length(tme)], model_sub = micromodel$Tz)

# -------------------------------
# Save outputs
# -------------------------------
out.path <- paste0(head.path, "micro_models/", height_label, "/")
dir.create(out.path, showWarnings = FALSE)
out.path <- paste0(out.path, yr, "/")
dir.create(out.path, showWarnings = FALSE)

saveRDS(dailyTz,   paste0(out.path, "microDaily_tile_", num, ".RDS"))

print("model saved.")