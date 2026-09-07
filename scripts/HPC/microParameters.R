#' Parameter script for microclimate modelling on a HPC.
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

#' ~~ Helper functions.
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

.modalReplace <- function(r) {
  vals <- values(r)
  valid_vals <- vals[!is.na(vals) & vals != 0]
  
  if (length(valid_vals) == 0) return(r)
  
  # Calculate mode manually
  tbl <- table(valid_vals)
  mval <- as.numeric(names(tbl)[which.max(tbl)])
  
  r[r == 0] <- mval
  return(r)
}

compute.daily <- function(tme, model_sub){
  
  udays <- unique(as.Date(tme))
  ndays <- length(udays)
  
  daily.min <- array(NA, dim = c(dim(model_sub)[c(1:2)], ndays))
  daily.max <- array(NA, dim = dim(daily.min))
  daily.mean <- array(NA, dim = dim(daily.min))
  
  # Loop over days and compute summaries
  for (i in seq_along(udays)) {
    day_indices <- which(as.Date(tme) == udays[i])
    daily_slice <- model_sub[,,day_indices]
    daily_slice[daily_slice==-Inf] <- 0
    
    daily.min[,,i] <- apply(daily_slice, c(1, 2), min, na.rm = TRUE)
    daily.max[,,i] <- apply(daily_slice, c(1, 2), max, na.rm = TRUE)
    daily.mean[,,i] <- apply(daily_slice, c(1, 2), mean, na.rm = TRUE)
  }
  daily.list <- list(
    dailyMin = daily.min, 
    dailyMax = daily.max, 
    dailyMean = daily.mean,
    day = udays)
  return(daily.list)
}

.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  if (class(r)[1] != "matrix") {
    if (dim(r)[3] > 1) {
      y<-as.array(r)
    } else y<-as.matrix(r,wide=TRUE)
  } else y<-r
  y
}

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

#' ~~ Required libraries.
packages <- c("dplyr", "tidyr", "readr", "microclimf",
              "myClim", "sf", "lubridate","terra")
load_packages(packages)

#' ~~ Set Pathways.
davies.cluster <- "/n/davies_lab/Users/btrew/"
setwd(davies.cluster)

in.data <- "data/"
out.data <- "output/"

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' [Define Site and Model Parameters]
res <- 5; dzd <- 1

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

dir.create(paste0(out.data,site.location,"/"), showWarnings = F)
dir.create(paste0(out.data,site.location,"/",sample.name), showWarnings = F)
dir.create(paste0(out.data,site.location,"/",sample.name,"/",year), showWarnings = F)
head.path <- paste0(out.data,site.location,"/",sample.name,"/",year,"/")
