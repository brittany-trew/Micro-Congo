#' ----------------------------------------------------------------------------
#' Pool summary statistics across all sites.
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

plot.out <- paste0(out.data,"plots/")

table.out <- paste0(out.data,"tables/")

metric_files <- list(
  "Mean daily vertical SD" = setNames(
    file.path(plot.out,
              paste0("Fig2A_", all.samples, ".tif")),
    all.samples),
  "Maximum daily vertical range" = setNames(
    file.path(plot.out,
              paste0("max_range_", all.samples, ".tif")),
    all.samples),
  "Mean vertical range on hottest 5% of days" = setNames(
    file.path(plot.out,
              paste0("mean_vertical_range_hottest5pct_",
                     all.samples, ".tif")),
    all.samples))


# Extract and pool valid cell values
get_pooled_values <- function(files) {
  
  vals <- unlist(lapply(files, function(f) {
    terra::values(terra::rast(f),
                  mat = FALSE)
  }), use.names = FALSE)
  
  vals[is.finite(vals)]
}


# Calculate statistics for one metric
summarise_metric <- function(files, metric_name) {
  
  # Keep only completed sites
  available <- file.exists(files)
  
  if (any(!available)) {
    message(metric_name,
            " — missing sites: ",
            paste(names(files)[!available], collapse = ", "))
  }
  
  files <- files[available]
  
  # Skip metric if no sites have finished
  if (length(files) == 0) {
    message(metric_name,
            " — no completed sites found.")
    return(NULL)
  }
  
  x <- get_pooled_values(files)
  
  if (length(x) == 0) {
    message(metric_name,
            " — completed rasters contain no valid values.")
    return(NULL)
  }
  
  qs <- quantile(x,
                 probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
                 na.rm = TRUE,
                 names = FALSE)
  
  data.frame(metric = metric_name,
             sites_completed = length(files),
             columns_pooled = length(x),
             mean = mean(x),
             sd = sd(x),
             minimum = min(x),
             percentile_05 = qs[1],
             percentile_25 = qs[2],
             median = qs[3],
             percentile_75 = qs[4],
             percentile_95 = qs[5],
             maximum = max(x))
}


# Summarise all available metrics
pooled_results <- lapply(names(metric_files), function(metric_name) {
  summarise_metric(files = metric_files[[metric_name]],
                   metric_name = metric_name)
})

pooled_results <- Filter(Negate(is.null),
                         pooled_results)

if (length(pooled_results) == 0) {
  stop("No completed raster files were available.")
}

pooled_statistics <- do.call(rbind,
                             pooled_results)

rownames(pooled_statistics) <- NULL

# Round statistics but retain exact counts
stat_columns <- setdiff(names(pooled_statistics),
                        c("metric",
                          "sites_completed",
                          "columns_pooled"))

pooled_statistics[stat_columns] <- round(
  pooled_statistics[stat_columns],
  digits = 3)

print(pooled_statistics)

write.csv(pooled_statistics,
          file.path(table.out,
                    "pooled_vertical_thermal_statistics.csv"),
          row.names = FALSE)
