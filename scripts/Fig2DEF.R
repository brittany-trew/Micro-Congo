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
#' Microclimate stability (SD of max temps through time)
for(b in 1:length(all.samples)){
  sample.name <- all.samples[[b]]
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  
  inpath  <- paste0(volumes3, sample.name, "/")
  outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
  
  chm           <- rast(paste0(head.path, "/chm.tif"))
  pai.layers    <- list.files(paste0(head.path, "/pai"))
  n.h           <- length(pai.layers)
  model.heights <- c(0.2,seq(1,(n.h-1)*1,by = 1))
  yr.seq        <- seq(2004, 2024, 1)
  
  # Use the first temperature raster as the working grid
  template <- rast(list.files(
    paste0(inpath, "dailyTemps/0.2m/"),
    pattern    = "DailyMax_.*\\.tif$",
    full.names = TRUE
  )[1])
  
  # Resample CHM once to the temperature grid
  chm_temp <- resample(chm, template)
  
  h.list <- vector("list", length(model.heights))
  
  for (i in 1:25) {
    
    hh <- model.heights[i]
    h  <- sprintf("%.1f", hh)
    
    files_max <- list.files(
      paste0(inpath, "dailyTemps/", h, "m/"),
      pattern    = "DailyMax_.*\\.tif$",
      full.names = TRUE
    )
    
    # No resampling, assuming temperature files share a grid
    r.stk <- rast(files_max)
    
    gridsd <- app(
      r.stk,
      function(x) sd(x[x >= 0], na.rm = TRUE)
    )
    
    gridsd <- ifel(chm_temp <= hh, NA, gridsd)
    
    h.list[[i]] <- gridsd
  }
  
  sd.stk <- rast(h.list)
  
  #' Save rasters.
  for(i in 1:nlyr(sd.stk)){
    hh <- model.heights[i]
    
    sd.yr <- sd.stk[[i]]
    writeRaster(sd.yr, paste0(outpath, "temporalSD_", hh, "m.tif"))
  }
}  


#' ----------------------------------------------------------------------------
#' Processing all site data for the heatmap
h.list <- list()
for (i in 1:length(model.heights)) {
  
  hh <- model.heights[i]
  print(hh)
  
  sample.list <- list()
  for(b in 1:length(all.samples)){
    sample.name <- all.samples[[b]]
    print(sample.name)
    outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")

    colsd_all <- rast(list.files(outpath, pattern = paste0("^temporalSD_",hh,"m.tif$"), full.names = TRUE))
    
    colsd_mean <- mean(colsd_all, na.rm = TRUE)
    
    sample.list[[b]] <- data.frame(
      height = hh,
      sd_v = values(colsd_mean, mat = FALSE)
    ) %>%
      dplyr::filter(!is.na(sd_v))
    
  }
  
  h.list[[i]] <- dplyr::bind_rows(sample.list)
}

sd_v <- dplyr::bind_rows(h.list)

sd_98 <- quantile(sd_v$sd_v, 0.98, na.rm = TRUE)

bin_width <- 0.1
cell_area_ha <- (5 * 5) / 10000  # 0.0025 ha

sd_heat <- sd_v %>%
  dplyr::filter(sd_v <= sd_98) %>%
  dplyr::mutate(
    sd_bin = floor(sd_v / bin_width) * bin_width + bin_width / 2
  ) %>%
  dplyr::count(height, sd_bin) %>%
  dplyr::mutate(
    area_ha = n * cell_area_ha
  )

#' Fig 2A: Heatmap Plot.
vh.plot <- ggplot(sd_heat, aes(
  x = sd_bin,
  y = height,
  fill = area_ha
)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = heat_cols,
    trans = "sqrt"
  ) +
  labs(
    x = "Temporal Stability (°C)",
    y = "Height",
    fill = "Area (ha)"
  ) +
  theme_classic()
vh.plot

ggsave(
  paste0(plot.out, "Fig2D.png"),
  plot = vh.plot,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 300
)

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

#' Derives final summary map: mean across heights.
message("Computing multi-year summaries...")

for(b in 1:length(all.samples)){
  sample.name <- all.samples[[b]]
  outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
  
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  template <- rast(paste0(plot.out, "Fig2C_",sample.name,".tif"))
  
  colsd_all <- rast(list.files(outpath, pattern = "^temporalSD_.*.tif$", full.names = TRUE))
  # MAP 1: Mean of vertical standard deviation in Daily Tmax
  mean_sd <- app(colsd_all, mean, na.rm = TRUE)
  mean_sd <- mask(mean_sd, template)

  # Remove values above the 98th percentile
  sd_q98 <- terra::global(
    mean_sd,
    fun = stats::quantile,
    probs = 0.999,
    na.rm = TRUE
  )[1, 1]
  mean_sd <- ifel(mean_sd > sd_q98, NA, mean_sd)
  
  plot(mean_sd)
  writeRaster(
    mean_sd,
    paste0(plot.out, "Fig2F_",sample.name,".tif"),
    overwrite = TRUE
  )
  message("  Fig 2C written for: ",sample.name,".")
  
}


temporal_all <- c(
  paste0(plot.out, "Fig2F_Imbalanga.tif"),
  paste0(plot.out, "Fig2F_MondoBai.tif"),
  paste0(plot.out, "Fig2F_MobaCluster.tif"),
  paste0(plot.out, "Fig2F_Lokoue.tif")
)

temporal_vals <- unlist(lapply(temporal_all, function(f) {
  values(rast(f), na.rm = TRUE)
}))

qs <- quantile(
  temporal_vals,
  probs = c(1/3, 2/3),
  na.rm = TRUE
)

temp_breaks <- c(-Inf, qs[[1]], qs[[2]], Inf)


for(b in seq_along(all.samples)) {
  
  sample.name <- all.samples[[b]]
  
  temp <- rast(paste0(plot.out, "Fig2F_", sample.name, ".tif"))

  # Reclassify to 1, 2, 3
  temp_class <- classify(
    temp,
    rcl = matrix(c(
      temp_breaks[[1]],      temp_breaks[[2]], 1,
      temp_breaks[[2]],  temp_breaks[[3]], 2,
      temp_breaks[[3]],  temp_breaks[[4]],      3
    ), ncol = 3, byrow = TRUE)
  )
  
  writeRaster(
    temp_class,
    paste0(plot.out, "TempClass_", sample.name, ".tif"),
    overwrite = TRUE
  )
  
}


