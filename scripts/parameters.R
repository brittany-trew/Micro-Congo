#' Manuscript: Vertical microclimate diversity shapes thermal refugia in Congolese rainforests.
#' Date: August 2026
#' Author: Brittany Trew (*brittany.trew@fas.harvard.edu*)
#' Description: XXXX
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' [Load Helper Functions.]
scripts.path <- "scripts/"
#' Load helper functions:
source(paste0(scripts.path,"utils.R"))
source(paste0(scripts.path,"helpers_general.R"))
source(paste0(scripts.path,"helpers_MicroDrivers.R"))

#' [Required R packages]
packages <- c("dplyr", "tidyr", "readr", "microclimf",
              "myClim", "ggplot2", "sf", 'lidR',"RCSF", "RMCC", "Rcpp", 
              "rgl", "luna", "geodata","sf","xml2","jsonlite",
              "httr", "ecmwfr", "mcera5", "lubridate", "curl",
              "mgcv", "colorspace", "dendextend", "cluster", "factoextra","terra",
              "sfsmisc", "canopyLazR", "ggrepel", "broom", "purrr",
              "patchwork")
load_packages(packages)

#' [Set Pathways.]
#' (1) LiDAR Point Clouds.
lidar.data <- "/Volumes/HD01/LIDAR/"

#' (2) Output and input locations for the data driving the microclimate model
mclidar.path <- "../mclidar/"
mclidar.out <- paste0(mclidar.path,"output/")
mclidar.in <- paste0(mclidar.path,"data/")

#' (3) Location of ERA5 Reanalysis Data 
clim.path <- paste0(mclidar.in,"era5/")

#' (4) Inout and output locations for the manuscript analysis (post microclimate modelling)
in.data <- "data/"
out.data <- "output/"

#' (5) Location for temporary files (ideally an external hard drive)
terraOptions(
  tempdir = "/Volumes/MicroMaze2/tmp",
  memfrac = 0.7
)

volumes3 <- "/Volumes/MicroMaze2/"
volumes <- "/Volumes/MicroMaze/"

#' [Define Site Parameters]
#' (1) Spatial resolution; vertical interval; extinction coefficient
res <- 5; dzd <- 0.5

#' (2) Site information.
site.summary <- read.csv(paste0(in.data,"Site_Info/site-summary.csv"))
site.location <- unique(site.summary$SiteName)
print(paste0("Site Location: ", site.location))
all.samples <- unique(site.summary$Sample)
print(all.samples)
sample.name <- all.samples[1] #' *CHANGE AS NEEDED*
num <- which(site.summary$Sample == sample.name)

year <- 2021
num <- which(site.summary$Year == year)
site.sample <- site.summary[num,]
input_month <- site.sample$Month

#' (3) Projection system.
defined_proj <- paste0("EPSG:",site.sample$EPSG[[1]])
print(defined_proj)

#' (4) Site-specific Pathways.
head.path <- paste0(mclidar.out,site.location,"/",sample.name,"/",year,"/")
all.tiles <- list.files(paste0(lidar.data,site.location,"/",sample.name,"/",year), full.names = T)

#' LiDAR data
las.path <- paste0(lidar.data,site.location,"/",sample.name,"/",year)
las.files <- list.files(las.path, pattern = "\\.(las|laz)$", full.names = TRUE)

#' [Create output paths - site specific]
#' Output path for Plant Area Density.
pad.path <- paste0(head.path,"pad/")
dir.create(pad.path, showWarnings = F)


#' [Standardised Plotting Parameters]
colour.p <- c(
  "1" = "#ca763b",
  "2" = "#852e47",
  "3" = "#333333",
  "4" = "#a09c33"
)

heat_cols <- c(
  "#f4efe9",
  "#dfc9ae",
  "#ca763b",
  "#a95643",
  "#852e47",
  "#333333"
)

others <- c(
  "#671b27",
  "#852e47",
  "#ca763b"
)

option1 <- c(
  "1" = "#5a6b44",
  "2" = "#e6d7c3",
  "3" = "#0d384e",
  "4" = "#c85a3e"
)

option2 <- c(
  "1" = "#264a31",
  "2" = "#f27d16",
  "3" = "#333333",
  "4" = "#6b854a"
)

option3 <- c(
  "1" = "#c85a3e",
  "2" = "#eea83b",
  "3" = "#004d47",
  "4" = "#a09c33"
)