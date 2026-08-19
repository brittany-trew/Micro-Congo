#' This script sets up file locations and general parameters for 
#' microclimate analysis of LiDAR sites in Odzala NP, Congo.
#' Written by Brittany Trew in 2026 for the MicroMaze Project, funded by HSFP.
#' ----------------------------------------------------------------------------

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' [Inital Setup.]
scripts.path <- "scripts/"
source(paste0(scripts.path,"utils.R")) # loads worker functions

#' ~~ Required libraries.
packages <- c("dplyr", "tidyr", "readr", "microclimf",
              "myClim", "ggplot2", "sf", 'lidR',"RCSF", "RMCC", "Rcpp", 
              "rgl", "luna", "geodata","sf","xml2","jsonlite",
              "httr", "ecmwfr", "mcera5", "lubridate", "curl",
              "mgcv", "colorspace", "dendextend", "cluster", "factoextra","terra",
              "sfsmisc", "canopyLazR", "ggrepel", "broom", "purrr",
              "patchwork")
load_packages(packages)

#' ~~ CPP Code:
sourceCpp(paste0(scripts.path,"CppFunctions.cpp"))  # compiles and links required C++ code

#' ~~ Set Pathways.
#' External drive pathways
volumes <- "/Volumes/MicroMaze/"
volumes2 <- "/Volumes/HD01/"
volumes3 <- "/Volumes/MicroMaze2/"

lidar.data <- paste0(volumes2,"LIDAR/")

#' Pathways to mclidar git
mclidar.path <- "../mclidar/"
mclidar.out <- paste0(mclidar.path,"output/")
mclidar.in <- paste0(mclidar.path,"data/")

shp.path <- paste0(mclidar.in,"shapefiles/")
clim.path <- paste0(mclidar.in,"era5/")

#' Pathways to congo git
in.data <- "data/"
out.data <- "output/"

#' Define location of temporary files
terraOptions(
  tempdir = paste0(volumes3,"tmp"),
  memfrac = 0.7
)
#' ----------------------------------------------------------------------------

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' [Define Site and Model Parameters]
res <- 5; dzd <- 0.5

#' ~~ Site Parameters
site.summary <- read.csv(paste0(in.data,"Site_Info/site-summary.csv"))
site.location <- unique(site.summary$SiteName)
print(paste0("Site Location: ", site.location))

all.samples <- unique(site.summary$Sample)
print(all.samples)
year <- 2021
num <- which(site.summary$Year == year)

site.sample <- site.summary[num,]
input_month <- site.sample$Month
defined_proj <- paste0("EPSG:",site.sample$EPSG[[1]])
print(defined_proj)

head.path <- paste0(mclidar.out,site.location,"/",sample.name,"/",year,"/")
all.tiles <- list.files(paste0(lidar.data,site.location,"/",sample.name,"/",year), full.names = T)

#' Point cloud files:
las.path <- paste0(lidar.data,site.location,"/",sample.name,"/",year)
las.files <- list.files(las.path, pattern = "\\.(las|laz)$", full.names = TRUE)
#' ----------------------------------------------------------------------------

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' [Plotting parameters]
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