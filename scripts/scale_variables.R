#' Scaling climate variables across sites.
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

#' Scaling variables across sites:
get_global_scale <- function(files) {
  
  total_sum <- 0
  total_sumsq <- 0
  total_n <- 0
  
  for (f in files) {
    r <- rast(f)
    r[!is.finite(r)] <- NA
    
    s  <- global(r, "sum", na.rm = TRUE)[, 1]
    ss <- global(r^2, "sum", na.rm = TRUE)[, 1]
    n  <- global(!is.na(r), "sum", na.rm = TRUE)[, 1]
    
    total_sum   <- total_sum + sum(s, na.rm = TRUE)
    total_sumsq <- total_sumsq + sum(ss, na.rm = TRUE)
    total_n     <- total_n + sum(n, na.rm = TRUE)
  }
  
  mu <- total_sum / total_n
  sd <- sqrt((total_sumsq - total_n * mu^2) / (total_n - 1))
  
  list(mean = mu, sd = sd)
}


inpath_all <- "/Volumes/MicroMaze2/"

#' Scale rasters (Needs to be done across sites)
files_tmin <- list.files(
  inpath_all,
  pattern = "yearlyMean_tmin_",
  recursive = TRUE,
  full.names = TRUE
)

files_tmean_sd <- list.files(
  inpath_all,
  pattern = "yearlySD_tmean_",
  recursive = TRUE,
  full.names = TRUE
)

files_q90_tmax <- list.files(
  inpath_all,
  pattern = "yearlyq90_tmax_",
  recursive = TRUE,
  full.names = TRUE
)

files_dr <- list.files(
  inpath_all,
  pattern = "DiurnalRange_",
  recursive = TRUE,
  full.names = TRUE
)

scale_tmin     <- get_global_scale(files_tmin)
scale_tmean_sd <- get_global_scale(files_tmean_sd)
scale_q90_tmax <- get_global_scale(files_q90_tmax)
scale_dr       <- get_global_scale(files_dr)

all.scales <- list("tmin" = scale_tmin,
                   "tmeanSD" = scale_tmean_sd,
                   "q90tmax" = scale_q90_tmax,
                   "meanDR" = scale_dr)
saveRDS(all.scales, paste0(out.data,"climate_distinctness/scaling.RDS"))
