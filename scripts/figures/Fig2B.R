# =============================================================================
# Temporal Thermal Summaries
# Calculates mean temporal variation across forest heights and classifies
# pooled temporal heterogeneity into tertiles.
# =============================================================================
terraOptions(tempdir = "/Volumes/MicroMaze2/tmp",
             memfrac = 0.7)


# Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

plot.out <- paste0(out.data, "plots/")
dir.create(plot.out,
           showWarnings = FALSE)

#' ----------------------------------------------------------------------------
#' Temporal microclimate heterogeneity
#' Mean temporal SD across all modelled heights

message("Computing multi-year temporal summaries...")

for (b in seq_along(all.samples)) {
  
  sample.name <- all.samples[[b]]
  outpath <- paste0(out.data,
                    sample.name,
                    "/verticalSummary/")
  
  temporal_files <- sort(list.files(
    outpath,
    pattern = "^temporalSD_.*\\.tif$",
    full.names = TRUE
  ))
  
  # Skip sites that have not finished processing
  if (length(temporal_files) == 0) {
    message("  No temporal SD files found for: ",
            sample.name,
            ". Skipping.")
    next
  }
  
  message("  Processing ",
          sample.name,
          " using ",
          length(temporal_files),
          " height layers.")
  
  # Load temporal SD rasters from all heights
  colsd_all <- rast(temporal_files)
  
  # Calculate mean temporal SD across heights
  mean_sd <- app(colsd_all,
                 mean,
                 na.rm = TRUE)
  
  # Apply the existing forest mask
  template_file <- paste0(plot.out,
                          "Fig2A_",
                          sample.name,
                          ".tif")
  template <- rast(template_file)
  mean_sd <- mask(mean_sd,
                  template)
  names(mean_sd) <- "mean_temporal_sd"
  plot(mean_sd)
  
  writeRaster(mean_sd,
              paste0(plot.out,
                     "Fig2B_",
                     sample.name,
                     ".tif"),
              overwrite = TRUE)
  
  message("  Fig 2B written for: ",
          sample.name,
          ".")
  
  rm(colsd_all,
     mean_sd,
     template)
  gc()
}

