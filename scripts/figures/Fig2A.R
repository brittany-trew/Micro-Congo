# =============================================================================
# Vertical Thermal Summaries
# Computes spatial maps of peak temperature and vertical thermal spread
# across the forest column for each year, then summarises across 2004-2024
# =============================================================================
terraOptions(
  tempdir = "/Volumes/MicroMaze2/tmp",
  memfrac = 0.7
)

# Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Vertical microclimate heterogeneity (SD of max temps)
inpath <- paste0(out.data, sample.name, "/verticalSummary/")
plot.out <- paste0(out.data,"plots/")

for(b in 1:length(all.samples)){
  sample.name <- all.samples[[b]]
  outpath <- paste0(out.data, sample.name, "/verticalSummary/")
  
  colsd_all <- rast(list.files(outpath, pattern = "^verticalSD_.*.tif$", full.names = TRUE))
  # MAP 1: Mean of vertical standard deviation in Daily Tmax
  mean_sd <- app(colsd_all, mean, na.rm = TRUE)

  # Remove values above the 95th percentile
  sd_q98 <- terra::global(
    mean_sd,
    fun = stats::quantile,
    probs = 0.95,
    na.rm = TRUE
  )[1, 1]
  mean_sd <- ifel(mean_sd > sd_q98, NA, mean_sd)
  
  plot(mean_sd)
  writeRaster(
    mean_sd,
    paste0(plot.out, "Fig2A_",sample.name,".tif"),
    overwrite = TRUE
  )
  message("  Fig 2A written for: ",sample.name,".")
  
}

for(b in 1:length(all.samples)){
  sample.name <- all.samples[[b]]
  outpath <- paste0(out.data, sample.name, "/verticalSummary/")
  
  colrange_all <- rast(list.files(outpath, pattern = "^verticalRange_.*.tif$", full.names = TRUE))
  # MAP 1: Mean of vertical standard deviation in Daily Tmax
  max_range <- app(colrange_all, max, na.rm = TRUE)
  # Remove values above the 98th percentile
  sd_q98 <- terra::global(
    max_range,
    fun = stats::quantile,
    probs = 0.95,
    na.rm = TRUE
  )[1, 1]
  max_range <- ifel(max_range > sd_q98, NA, max_range)
  
  plot(max_range)
  writeRaster(
    max_range,
    paste0(plot.out, "max_range_",sample.name,".tif"),
    overwrite = TRUE
  )
  message("  Max range written for: ",sample.name,".")
}

for (b in seq_along(all.samples)) {
  
  sample.name <- all.samples[[b]]
  
  outpath <- file.path(
    out.data,
    sample.name,
    "verticalSummary"
  )
  
  # Find and explicitly sort yearly files
  max_files <- sort(list.files(
    outpath,
    pattern = "^colmax_[0-9]{4}\\.tif$",
    full.names = TRUE
  ))
  
  range_files <- sort(list.files(
    outpath,
    pattern = "^verticalRange_[0-9]{4}\\.tif$",
    full.names = TRUE
  ))
  
  # Check that the same years are available for both variables
  max_years <- sub(
    "^colmax_([0-9]{4})\\.tif$",
    "\\1",
    basename(max_files)
  )
  
  range_years <- sub(
    "^verticalRange_([0-9]{4})\\.tif$",
    "\\1",
    basename(range_files)
  )
  
  stopifnot(
    length(max_files) > 0,
    length(range_files) > 0,
    identical(max_years, range_years)
  )
  
  # Load daily raster stacks
  colmax_all <- rast(max_files)
  colrange_all <- rast(range_files)
  
  # Confirm matching geometry and number/order of days
  stopifnot(nlyr(colmax_all) == nlyr(colrange_all))
  
  # Spatial mean column Tmax for each day
  daily_forest_tmax <- terra::global(
    colmax_all,
    fun = "mean",
    na.rm = TRUE
  )[, 1]
  
  # Overall hottest 5% of days
  hot_cutoff <- quantile(
    daily_forest_tmax,
    probs = 0.95,
    na.rm = TRUE
  )
  
  hot_day_indices <- which(
    daily_forest_tmax >= hot_cutoff
  )
  
  message(
    "  ", sample.name, ": ",
    length(hot_day_indices),
    " hot days selected from ",
    nlyr(colmax_all)
  )
  
  # Select matching vertical-range layers
  ranges_hot_days <- colrange_all[[hot_day_indices]]
  
  # Mean vertical range during hot days for each forest column
  mean_range_hot_days <- mean(
    ranges_hot_days,
    na.rm = TRUE
  )
  
  hot_range_q95 <- terra::global(
    mean_range_hot_days,
    fun = stats::quantile,
    probs = 0.95,
    na.rm = TRUE
  )[1, 1]
  
  mean_range_hot_days <- ifel(
    mean_range_hot_days > hot_range_q95,
    NA,
    mean_range_hot_days
  )
  names(mean_range_hot_days) <- "mean_range_hot_days"
  plot(mean_range_hot_days)
  
  writeRaster(
    mean_range_hot_days,
    file.path(
      plot.out,
      paste0(
        "mean_vertical_range_hottest5pct_",
        sample.name,
        ".tif"
      )
    ),
    overwrite = TRUE
  )
  
  message("  Hot-day vertical range written for: ", sample.name)
}
