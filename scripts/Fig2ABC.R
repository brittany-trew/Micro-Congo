# =============================================================================
# Vertical Thermal Summaries
# Computes spatial maps of peak temperature and vertical thermal spread
# across the forest column for each year, then summarises across 2004-2024
# =============================================================================
terraOptions(
  tempdir = "/Volumes/MicroMaze2/tmp",
  memfrac = 0.7
)

# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

plot.out <- paste0(out.data,"plots/")

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Vertical microclimate heterogeneity (SD of max temps)
sample.name <- all.samples[[4]]
head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")

inpath  <- paste0(volumes3, sample.name, "/")
outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

yr.seq        <- seq(2004, 2024, 1)

#' Loop for each year.
#' Loads daily max for each year. 
#' Calculates summaries of the vertical columns for each day. 
#' and saves as yearly files: eg., verticalSD_2018.tif.

for (i in 1:length(yr.seq)) {
  yr <- yr.seq[[i]]
  message("Processing ", yr, " ...")
  
  files_max <- list.files(
    inpath,
    pattern    = paste0("DailyMax_", yr, "\\.tif$"),
    recursive  = TRUE,
    full.names = TRUE
  )
  
  message("  Found ", length(files_max), " height bands: ",
          paste(sort(as.numeric(get_height(files_max))), collapse = ", "), "m")
  
  heights <- as.numeric(sub(".*/dailyTemps/([0-9.]+)m/.*", "\\1", files_max))
  files_max <- files_max[order(heights)]
  
  # Max, Min and SD of vertical columns.
  minmax <- collapse_vertical(files_max, chm)
  gc()
  # Daily vertical spread in peak temperatures
  r_range <- minmax$colmax - minmax$colmin

  writeRaster(minmax$colmax, paste0(outpath, "colmax_", yr, ".tif"), overwrite = TRUE)
  writeRaster(r_range,  paste0(outpath, "verticalRange_",  yr, ".tif"), overwrite = TRUE)
  writeRaster(minmax$colsd,  paste0(outpath, "verticalSD_",  yr, ".tif"), overwrite = TRUE)
  
  message("Output written for ", yr)
  rm(minmax, r_range)
  gc()
}

#' Fig. 2A: The distribution of vertical microclimate heterogeneity of each column 
#' as the standard deviation of daily maximum temperatures for each vertical column
yr.list <- list()
for (i in 1:length(yr.seq)) {
  
  yr <- yr.seq[[i]]
  print(yr)
  
  sample.list <- list()
  for(b in 1:length(all.samples)){
    sample.name <- all.samples[[b]]
    print(sample.name)
    outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
    
    if(yr == 2007){next()} # *FIX*
    colsd_all <- rast(list.files(outpath, pattern = paste0("^verticalSD_",yr,".tif$"), full.names = TRUE))
    
    colsd_mean <- mean(colsd_all, na.rm = TRUE)
    
    sample.list[[b]] <- data.frame(
      year = yr,
      sd_v = values(colsd_mean, mat = FALSE)
    ) %>%
      dplyr::filter(!is.na(sd_v))
    
  }
  
  yr.list[[i]] <- dplyr::bind_rows(sample.list)
}

sd_v <- dplyr::bind_rows(yr.list)

sd_98 <- quantile(sd_v$sd_v, 0.98, na.rm = TRUE)

bin_width <- 0.1
cell_area_ha <- (5 * 5) / 10000  # 0.0025 ha

sd_heat <- sd_v %>%
  dplyr::filter(sd_v <= sd_98) %>%
  dplyr::mutate(
    sd_bin = floor(sd_v / bin_width) * bin_width + bin_width / 2
  ) %>%
  dplyr::count(year, sd_bin) %>%
  dplyr::mutate(
    area_ha = n * cell_area_ha
  )

#' Fig 2A: Heatmap.
vh.plot <- ggplot(sd_heat, aes(
  x = year,
  y = sd_bin,
  fill = area_ha
)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = heat_cols,
    trans = "sqrt"
  ) +
  labs(
    x = "Year",
    y = "Vertical Heterogeneity (°C)",
    fill = "Area (ha)"
  ) +
  theme_classic()

ggsave(
  paste0(plot.out, "Fig2A.png"),
  plot = vh.plot,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 300
)


#' Derives final summary map: mean of each year.
message("Computing multi-year summaries...")

for(b in 1:length(all.samples)){
  sample.name <- all.samples[[b]]
  outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
  
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  chm  <- rast(paste0(head.path, "/chm.tif"))
  
  colsd_all <- rast(list.files(outpath, pattern = "^verticalSD_.*.tif$", full.names = TRUE))
  # MAP 1: Mean of vertical standard deviation in Daily Tmax
  mean_sd <- app(colsd_all, mean, na.rm = TRUE)
  mean_sd <- ifel(chm < 0.2, NA, mean_sd)
  
  # Remove values above the 98th percentile
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
    paste0(plot.out, "Fig2C_",sample.name,".tif"),
    overwrite = TRUE
  )
  message("  Fig 2C written for: ",sample.name,".")
  
}


temporal_all <- c(
  paste0(plot.out, "Fig2C_Imbalanga.tif"),
  paste0(plot.out, "Fig2C_MondoBai.tif"),
  paste0(plot.out, "Fig2C_MobaCluster.tif"),
  paste0(plot.out, "Fig2C_Lokoue.tif")
)

temporal_vals <- unlist(lapply(temporal_all, function(f) {
  values(rast(f), na.rm = TRUE)
}))

quantile(
  temporal_vals,
  probs = c(1/3, 2/3),
  na.rm = TRUE
)




#' Extras: 
colmax_all <- rast(list.files(outpath, pattern = "^colmax_.*.tif$", full.names = TRUE))
range_all  <- rast(list.files(outpath, pattern = "^verticalRange_.*.tif$",  full.names = TRUE))

# MAP 2: Absolute maximum Tmax in the column across all years
# "What is the worst-case peak heat exposure at this location?"
tmax_absolute <- app(colmax_all, max, na.rm = TRUE)
tmax_q90 <- app(colmax_all, fun = function(x) quantile(x, 0.90, na.rm = TRUE))
tmax_q90 <- ifel(chm < 0.2, NA, tmax_q90)
plot(tmax_q90)
writeRaster(
  tmax_q90,
  paste0(outpath, "Map2_Verticalq90DailyTmax_2004to2024.tif"),
  overwrite = TRUE
)
message("  Map 2 written: Q90 in Vertical Daily Tmax")

# MAP 3: Mean vertical daily range max temperatures across all years
mean_range <- app(range_all, mean, na.rm = TRUE)
mean_range <- ifel(chm < 0.2, NA, mean_range)
plot(mean_range)
writeRaster(
  mean_range,
  paste0(outpath, "Map3_MeanVerticalDailyRangeTmax_2004to2024.tif"),
  overwrite = TRUE
)
message("  Map 3 written: Mean Vertical Daily Range Tmax")

# MAP 5: Mean vertical spread on the hottest days (top 10% of column Tmax)
# "When heat load is greatest, how much escape does the column offer?"
# Threshold is computed per pixel across all days and years
hot_threshold        <- tmax_q90
hot_mask             <- colmax_all >= hot_threshold  # broadcasts across stack
range_hot            <- range_all
range_hot[!hot_mask] <- NA
range_on_hot_days    <- app(range_hot, mean, na.rm = TRUE)
writeRaster(
  range_on_hot_days,
  paste0(outpath, "Map5_RangeOnHottestDays_2004to2024.tif"),
  overwrite = TRUE
)
message("  Map 5 written: Vertical range on hottest days")
