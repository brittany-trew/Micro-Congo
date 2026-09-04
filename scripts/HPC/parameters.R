#' PARAMETER SCRIPT FOR THE LIDAR to MICROCLIMATE WORKFLOW
#' Written by Brittany Trew in 2025 for the MicroMaze HFSP Project.
#####################################################################
#####################################################################

#' [Inital Setup.]
scripts.path <- "scripts/"
source(paste0(scripts.path,"utils.R")) # loads worker functions
source(paste0(scripts.path,"data_prep.R"))

#' ~~ Required libraries.
packages <- c("dplyr", "tidyr", "readr", "microclimf",
              "myClim", "ggplot2", 'lidR',"RCSF", "RMCC", "Rcpp", 
              "rgl", "luna", "geodata","sf","xml2","jsonlite",
              "httr", "ecmwfr", "mcera5", "lubridate", "curl",
              "mgcv","terra")
load_packages(packages)

sourceCpp(paste0(scripts.path,"CppFunctions.cpp"))  # compiles and links required C++ code

#' ~~ Set Pathways.
davies.cluster <- "/n/davies_lab/Users/btrew/"
setwd(davies.cluster)

in.data <- "data/"
out.data <- "output/"

####################################################################
####################################################################

#' [Define Site and Model Parameters]
#' ~~ Site Parameters
site.summary <- read.csv(paste0(in.data,"site-summary.csv"))
all.sites <- unique(site.summary$SiteName)
print(all.sites)
site.location <- all.sites[3] #' *CHANGE AS NEEDED*
print(paste0("Site Location: ", site.location))

num <- which(site.summary$SiteName == site.location)
site <- site.summary[num,]
all.samples <- unique(site$Sample)
print(all.samples)
sample.name <- all.samples[2] #' *CHANGE AS NEEDED*
num <- which(site.summary$Sample == sample.name)

site.sample <- site.summary[num,]
print(site.sample)
year <- 2021 #' *CHANGE AS NEEDED*
num <- which(site.sample$Year == year)

site.sample <- site.sample[num,]
input_month <- site.sample$Month
defined_proj <- paste0("EPSG:",site.sample$EPSG)
sample.name <- site.sample$Sample
print(paste0("Site Name: ",sample.name))

res <- 5; dzd <- 1

dir.create(paste0(out.data,site.location,"/"), showWarnings = F)
dir.create(paste0(out.data,site.location,"/",sample.name), showWarnings = F)
dir.create(paste0(out.data,site.location,"/",sample.name,"/",year), showWarnings = F)
head.path <- paste0(out.data,site.location,"/",sample.name,"/",year,"/")

#' Create or load the bounding box.
bbox.path <-paste0(in.data,"bounding_box/",site.location,"/",sample.name,"/")
bbox.file <<- paste0(bbox.path,sample.name,".shp")
bbox_vect <- vect(bbox.file)

#' Template raster.
bbox_raster <- rast(ext(bbox_vect), resolution = res, crs = terra::crs(bbox_vect))
r_bbox <- rasterize(bbox_vect, bbox_raster, field = 0)

