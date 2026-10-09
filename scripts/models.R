# PAI and thermal heterogeneity
# Four models, all fitted by lm with log-transformed responses.
# Uncertainty: spatial blocks resampled WITHIN sites, at 300 m and 500 m.
# All eligible observations are fitted. Only variogram diagnostics are sampled.
# No GLS, stratum-identity terms, or interaction models in this script.
#
# Model 1: whole-column PAI -> whole-profile vertical heterogeneity
# Model 2: whole-column PAI -> whole-profile temporal heterogeneity
# Model 3: whole-column PAI -> within-stratum vertical heterogeneity
# Model 4: matching stratum PAI -> within-stratum vertical heterogeneity
# Models 3 and 4 use exactly the same cell-stratum rows.
# PAI is an input to the temperature model: these are descriptive associations
# among model inputs and outputs, not independent empirical validation.

# 1. Setup ------------------------------------------------------------------
source("scripts/parameters.R")

block_sizes_m <- c(300, 500)
bootstrap_B <- 4999L
bootstrap_seed <- 123L
variogram_max_cells <- 2000L
# FALSE: omit stratum variograms only. Whole-profile variograms are still saved.
# Stratum variograms use separate site-stratum groups to avoid repeated locations.
run_stratum_variograms <- TRUE
plot_dir <- file.path(out.data, "plots", "pai_lm_final")
table_dir <- file.path(out.data, "tables", "pai_lm_final")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Extract paired raster values -------------------------------------------
# Preserves the supplied processPAD(pad.r, chm_mosaic) call.
# PAD layers must represent the dzd intervals given in their names.
# A single total-PAI raster is unchanged by mean(); a stack is averaged across
# its layers, preserving the existing workflow. It is not averaged over cells.
# Heterogeneity masks determine which structural-class cells are retained.

site.df <- function(mclidar.out, all.samples, inpath){
  site.list <- list()
  
  for(i in seq_along(all.samples)){
    sample.name <- all.samples[[i]]
    head.path <- paste0(mclidar.out,"OdzalaCongo/",sample.name,"/2021/")
    
    pai <- terra::rast(paste0(head.path,"pai/PAI_0.5m_to_canopy.tif"))
    
    mean.pai <- mean(pai)
    mean.pai[mean.pai==0]<- NA
    
    if(!file.exists(paste0(inpath,sample.name,".tif"))){
      next()
    }
    r <- terra::rast(paste0(inpath,sample.name,".tif"))
    terra::compareGeom(mean.pai, r, stopOnError = TRUE)
    paired <- c(mean.pai, r)
    names(paired) <- c("pai", "var")
    
    dat <- as.data.frame(paired, xy = TRUE, cells = TRUE, na.rm = FALSE) %>%
      dplyr::mutate(
        sample = as.character(sample.name),
        x_km = (x - terra::xmin(mean.pai)) / 1000,
        y_km = (y - terra::ymin(mean.pai)) / 1000
      )
    
    site.list[[sample.name]] <- dat |>
      dplyr::select(sample, cell, x_km, y_km, pai, var) |>
      dplyr::filter(is.finite(pai), is.finite(var), var > 0)
    
  }
  df <- dplyr::bind_rows(site.list)
  return(df)
}

site.df.strata <- function(mclidar.out, all.samples, inpath,
                           response = c("vertical", "temporal")) {
  response <- match.arg(response)
  raster_pattern <- switch(
    response,
    vertical = "^VerticalSDinDailyTmax_2004to2024_.*\\.tif$",
    temporal = "^TemporalSDin[[:space:]]*DailyTmax_2004to2024_.*\\.tif$"
  )
  # The temporal pattern accepts a space between 'in' and 'DailyTmax'.
  # Both responses must use the same class and height suffix, e.g. _f1_0to5.tif.
  
  site.list <- list()
  
  for (i in seq_along(all.samples)) {
    
    sample.name <- all.samples[[i]]
    head.path <- paste0(mclidar.out, "OdzalaCongo/", sample.name, "/2021/")
    
    files <- list.files(
      file.path(inpath, sample.name, "rasters"),
      pattern = raster_pattern,
      full.names = TRUE,
      recursive = TRUE
    )
    
    if (length(files) == 0L) {
      message("No ", response, " stratum rasters found for ", sample.name)
      next
    }
    
    pai_total <- mean(
      terra::rast(paste0(head.path, "pai/PAI_0.5m_to_canopy.tif"))
    )
    pai_total[pai_total == 0] <- NA
    
    chm_mosaic <- terra::rast(paste0(head.path, "chm.tif"))
    
    all.pad <- list.files(
      file.path(head.path, "pad"),
      pattern = "\\.tif$",
      full.names = TRUE)
    
    pad.r <- lapply(all.pad, terra::rast)
    pad0 <- processPAD(pad.r, chm_mosaic)
    
    # Height boundaries from PAD layer names
    lower <- as.numeric(
      sub("^PAD_([0-9.]+)_.*$", "\\1", names(pad0)))
    upper <- as.numeric(
      sub("^PAD_[0-9.]+_([0-9.]+)m$", "\\1", names(pad0)))
    
    if (any(!is.finite(lower)) || any(!is.finite(upper))) {
      stop("Could not parse PAD layer heights for ", sample.name)
    }
    strata.list <- list()
    
    for (j in seq_along(files)) {
      
      filename <- basename(files[[j]])
      heterogeneity <- terra::rast(files[[j]])
      if (terra::nlyr(heterogeneity) != 1L) stop("Expected one response layer: ", filename)
      
      # Structural class from filename, not the loop index
      structural.class <- as.integer(
        sub(".*_f([0-9]+)_.*$", "\\1", filename)
      )
      
      heights <- as.numeric(
        strsplit(
          sub(".*_([0-9.]+to[0-9.]+)\\.tif$", "\\1", filename),
          "to"
        )[[1]])
      
      if (length(heights) != 2L || any(!is.finite(heights)) ||
          !is.finite(structural.class) || heights[2] <= heights[1]) {
        stop("Could not parse structural class / height interval: ", filename)
      }
      # Select PAD layers entirely within this stratum
      zslice <- which(lower >= heights[1] & upper <= heights[2])
      if (length(zslice) == 0L) stop("No PAD layers selected: ", filename)
      if (any(abs((upper[zslice] - lower[zslice]) - dzd) > 1e-8) ||
          abs(sum(upper[zslice] - lower[zslice]) - diff(heights)) > 1e-8) {
        stop("Selected PAD layers do not cover the interval with dzd bins: ", filename)
      }
      pad_subset <- pad0[[zslice]]
      
      # Total PAI within the stratum
      PAI <- sum(pad_subset) * dzd
      PAI[PAI == 0] <- NA
      
      # Check that cell locations match
      terra::compareGeom(PAI, heterogeneity, stopOnError = TRUE)
      terra::compareGeom(PAI, pai_total, stopOnError = TRUE)
      
      paired <- c(PAI, heterogeneity, pai_total)
      names(paired) <- c("pai", "var", "pai_total")
      
      dat <- as.data.frame(paired, xy = TRUE, cells = TRUE) %>%
        dplyr::filter(is.finite(pai), is.finite(pai_total), is.finite(var), var > 0) %>%
        dplyr::mutate(
          sample = as.character(sample.name),
          structural_class = structural.class,
          zone_start = heights[1],
          zone_end = heights[2],
          zone_id = paste0(
            "f", structural.class, "_", zone_start, "to", zone_end
          ),
          x_km = (x - terra::xmin(paired)) / 1000,
          y_km = (y - terra::ymin(paired)) / 1000
        )
      
      strata.list[[j]] <- dat
    }
    
    site.list[[sample.name]] <- dplyr::bind_rows(strata.list)
  }
  
  return(dplyr::bind_rows(site.list))
}

# 3. Prepare observations ---------------------------------------------------
prepare_data <- function(dat, predictors, strata = FALSE) {
  required <- c("sample", "cell", "x_km", "y_km", "var", predictors)
  if (strata) required <- c(required, "zone_id")
  if (!all(required %in% names(dat))) stop("Missing extracted data columns")
  keep <- !is.na(dat$sample) & is.finite(dat$var) & dat$var > 0 &
    is.finite(dat$x_km) & is.finite(dat$y_km)
  for (p in predictors) keep <- keep & is.finite(dat[[p]])
  if (strata) keep <- keep & !is.na(dat$zone_id)
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) < 3L) stop("Not enough eligible observations")
  dat$sample <- factor(dat$sample)
  keys <- c("sample", "cell", if (strata) "zone_id")
  if (anyDuplicated(dat[keys])) stop("Duplicate extracted observations: check input files")
  dat <- dplyr::arrange(dat, sample, cell)
  dat
}

# 4. Bootstrap a simple PAI slope ------------------------------------------
# This helper supports ONLY log(var) ~ one predictor.
# Block summaries yield the same slope as refitting every resampled row.
# Multiple strata at a horizontal cell stay together in the same block.
bootstrap_slope <- function(dat, predictor, block_size_m, B, seed) {
  if (block_size_m <= 0 || B < 2L) stop("Invalid bootstrap settings")
  d <- dat
  d$.x <- d[[predictor]] - mean(d[[predictor]])
  d$.y <- log(d$var)
  d$block_id <- paste(d$sample,
                      floor(d$x_km * 1000 / block_size_m),
                      floor(d$y_km * 1000 / block_size_m), sep = "_")
  stats <- d %>%
    dplyr::group_by(sample, block_id) %>%
    dplyr::summarise(n = dplyr::n(), sx = sum(.x), sy = sum(.y),
                     sxx = sum(.x^2), sxy = sum(.x * .y), .groups = "drop") %>%
    dplyr::arrange(sample, block_id)
  by_site <- split(seq_len(nrow(stats)), stats$sample, drop = TRUE)
  if (any(lengths(by_site) < 2L)) stop("Need at least two occupied blocks per site")
  values <- as.matrix(stats[c("n", "sx", "sy", "sxx", "sxy")])
  calculate <- function(totals) {
    denom <- totals["sxx"] - totals["sx"]^2 / totals["n"]
    if (!is.finite(denom) || denom <= .Machine$double.eps * max(1, totals["sxx"])) {
      return(NA_real_)
    }
    unname((totals["sxy"] - totals["sx"] * totals["sy"] / totals["n"]) / denom)
  }
  set.seed(seed)
  draws <- replicate(B, {
    selected <- unlist(lapply(by_site, function(ids) {
      ids[sample.int(length(ids), length(ids), replace = TRUE)]
    }), use.names = FALSE)
    calculate(colSums(values[selected, , drop = FALSE]))
  })
  if (any(!is.finite(draws))) {
    stop("An unestimable slope occurred in a bootstrap resample; inspect predictor coverage")
  }
  counts <- d %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_observations = dplyr::n(),
      n_horizontal_cells = dplyr::n_distinct(cell),
      n_blocks = dplyr::n_distinct(block_id), .groups = "drop")
  counts$block_size_m <- block_size_m
  list(draws = draws, CI = as.numeric(quantile(draws, c(0.025, 0.975))),
       estimate = calculate(colSums(values)), block_counts = counts)
}

# 5. Diagnostics -----------------------------------------------------------
save_diagnostics <- function(model, id) {
  residual <- rstandard(model)
  png(file.path(plot_dir, paste0(id, "_QQ.png")),
      width = 6, height = 6, units = "in", res = 300)
  tryCatch({
    qqnorm(residual, pch = 16, cex = 0.3, main = id)
    qqline(residual, col = "#852e47", lwd = 2)
  }, finally = dev.off())
  png(file.path(plot_dir, paste0(id, "_residuals.png")),
      width = 6, height = 6, units = "in", res = 300)
  tryCatch({
    plot(fitted(model), residual, pch = 16, cex = 0.3,
         col = adjustcolor("black", alpha.f = 0.1),
         xlab = "Fitted log(thermal heterogeneity)",
         ylab = "Standardised residual", main = id)
    abline(h = 0, lty = 2)
  }, finally = dev.off())
}

save_variograms <- function(model, dat, id, strata = FALSE) {
  stopifnot(nobs(model) == nrow(dat))
  d <- dat
  d$.residual <- as.numeric(residuals(model))
  grouping <- if (strata) interaction(d$sample, d$zone_id, drop = TRUE) else d$sample
  parts <- split(d, grouping, drop = TRUE)
  set.seed(bootstrap_seed)
  for (nm in names(parts)) {
    v <- parts[[nm]]
    if (anyDuplicated(v[c("x_km", "y_km")])) {
      stop("Repeated locations in variogram group: ", nm)
    }
    v <- v[sample.int(nrow(v), min(variogram_max_cells, nrow(v))), , drop = FALSE]
    if (nrow(v) < 3L) next
    vg <- gstat::variogram(.residual ~ 1, locations = ~ x_km + y_km,
                           data = v, cutoff = 0.5, width = 0.01)
    filename <- paste0(id, "_variogram_", gsub("[^A-Za-z0-9_.-]", "_", nm))
    write.csv(as.data.frame(vg), file.path(table_dir, paste0(filename, ".csv")), row.names = FALSE)
    if (nrow(vg) == 0L) next
    png(file.path(plot_dir, paste0(filename, ".png")),
        width = 6, height = 6, units = "in", res = 300)
    tryCatch({
      plot(vg$dist, vg$gamma, type = "b", pch = 1,
           xlab = "Distance between cells (km)",
           ylab = "Residual semivariance", main = nm)
    }, finally = dev.off())
  }
}

# 6. Fit, summarise, and save one model -------------------------------------
analyse_model <- function(dat, predictor, id, response, predictor_label,
                          strata = FALSE) {
  # Dataset must already be prepared. na.fail prevents silently changing rows.
  formula <- stats::reformulate(predictor, response = "log(var)")
  model <- stats::lm(formula, data = dat, na.action = stats::na.fail)
  stopifnot(nobs(model) == nrow(dat))
  slope <- unname(coef(model)[predictor])
  if (!is.finite(slope)) stop("Full-data PAI slope cannot be estimated")
  summary_model <- summary(model)
  save_diagnostics(model, id)
  if (!strata || run_stratum_variograms) save_variograms(model, dat, id, strata)
  
  bootstraps <- lapply(block_sizes_m, function(size) {
    bootstrap_slope(dat, predictor, size, bootstrap_B, bootstrap_seed)
  })
  for (boot in bootstraps) {
    if (!isTRUE(all.equal(slope, boot$estimate, tolerance = 1e-8))) {
      stop("Bootstrap and fitted model disagree about the predictor / slope")
    }
  }
  n_cells <- nrow(dplyr::distinct(dat, sample, cell))
  results <- dplyr::bind_rows(lapply(seq_along(block_sizes_m), function(k) {
    ci <- bootstraps[[k]]$CI
    data.frame(
      model = id, response = response, predictor = predictor_label,
      n_observations = nobs(model), n_horizontal_cells = n_cells,
      n_sites = dplyr::n_distinct(dat$sample),
      R2_log = summary_model$r.squared,
      adjusted_R2_log = summary_model$adj.r.squared,
      slope_log = slope, lower_95_log = ci[1], upper_95_log = ci[2],
      percent_change = 100 * expm1(slope),
      lower_95_percent = 100 * expm1(ci[1]), upper_95_percent = 100 * expm1(ci[2]),
      block_size_m = block_sizes_m[k], bootstrap_replicates = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      CI_method = "95% percentile spatial-block bootstrap within sites"
    )
  }))
  counts <- dplyr::bind_rows(lapply(bootstraps, function(x) x$block_counts))
  counts$model <- id
  write.csv(results, file.path(table_dir, paste0(id, "_results.csv")), row.names = FALSE)
  write.csv(counts, file.path(table_dir, paste0(id, "_blocks.csv")), row.names = FALSE)
  print(results[c("model", "block_size_m", "percent_change",
                  "lower_95_percent", "upper_95_percent", "R2_log")])
  list(model = model, results = results, block_counts = counts,
       bootstrap_draws = lapply(bootstraps, function(x) x$draws))
}

# 7. Model 1: whole-column PAI -> whole-profile VERTICAL heterogeneity --------
dat_profile_vertical <- prepare_data(
  site.df(mclidar.out, all.samples, paste0(out.data, "/plots/Fig2A_")), "pai")
results_profile_vertical <- analyse_model(
  dat_profile_vertical, "pai", "model1_profile_vertical",
  "Whole-profile vertical thermal heterogeneity", "Whole-column PAI")
m_profile_vertical <- results_profile_vertical$model

# 8. Model 2: whole-column PAI -> whole-profile TEMPORAL heterogeneity --------
# Fig2B is the temporal heterogeneity raster prefix from the existing workflow.
dat_profile_temporal <- prepare_data(
  site.df(mclidar.out, all.samples, paste0(out.data, "/plots/Fig2B_")), "pai")
results_profile_temporal <- analyse_model(
  dat_profile_temporal, "pai", "model2_profile_temporal",
  "Whole-profile temporal thermal heterogeneity", "Whole-column PAI")
m_profile_temporal <- results_profile_temporal$model

# 9. Extract stratum data ONCE ---------------------------------------------
dat_strata <- prepare_data(
  site.df.strata(mclidar.out, all.samples, inpath = out.data, response = "vertical"),
  predictors = c("pai_total", "pai"), strata = TRUE)

# Repeated horizontal cells are intentional: a cell can contribute several
# strata. Both models use all these rows; the bootstrap keeps them together.
saveRDS(dat_strata, file.path(table_dir, "paired_stratum_data_vertical.rds"))

# 10. Model 3: whole-column PAI -> within-stratum VERTICAL heterogeneity ------
results_total <- analyse_model(
  dat_strata, "pai_total", "model3_strata_whole_column_PAI",
  "Within-stratum vertical thermal heterogeneity", "Whole-column PAI",
  strata = TRUE)
m_total <- results_total$model

# 11. Model 4: matching stratum PAI -> within-stratum VERTICAL heterogeneity --
results_stratum <- analyse_model(
  dat_strata, "pai", "model4_strata_matching_PAI",
  "Within-stratum vertical thermal heterogeneity", "Matching stratum PAI",
  strata = TRUE)
m_stratum <- results_stratum$model

# 12. Extract stratum data ONCE ---------------------------------------------
dat_strata_temporal <- prepare_data(
  site.df.strata(mclidar.out, all.samples, inpath = out.data, response = "temporal"),
  predictors = c("pai_total", "pai"), strata = TRUE)
saveRDS(dat_strata_temporal, file.path(table_dir, "paired_temporal_stratum_data.rds"))

# 13. Model 5: whole-column PAI -> within-stratum TEMPORAL heterogeneity -------
results_total_temporal <- analyse_model(
  dat_strata_temporal, "pai_total", "model5_temporal_strata_whole_column_PAI",
  "Within-stratum temporal thermal heterogeneity", "Whole-column PAI",
  strata = TRUE)
m_total_temporal <- results_total_temporal$model

# 14. Model 6: matching stratum PAI -> within-stratum TEMPORAL heterogeneity ---
results_stratum_temporal <- analyse_model(
  dat_strata_temporal, "pai", "model6_temporal_strata_matching_PAI",
  "Within-stratum temporal thermal heterogeneity", "Matching stratum PAI",
  strata = TRUE)
m_stratum_temporal <- results_stratum_temporal$model


# 12. Final tables and saved fitted models ----------------------------------
all_results <- list(
  profile_vertical = results_profile_vertical,
  profile_temporal = results_profile_temporal,
  vertical_strata_whole_column_PAI = results_total,
  vertical_strata_matching_PAI = results_stratum,
  temporal_strata_whole_column_PAI = results_total_temporal,
  temporal_strata_matching_PAI = results_stratum_temporal)
all_model_results <- dplyr::bind_rows(lapply(all_results, function(x) x$results))
all_block_counts <- dplyr::bind_rows(lapply(all_results, function(x) x$block_counts))
write.csv(all_model_results, file.path(table_dir, "all_model_results.csv"), row.names = FALSE)
write.csv(all_block_counts, file.path(table_dir, "all_model_block_counts.csv"), row.names = FALSE)

# ONLY Models 3 and 4 have the same response AND identical observations.
# Compare their descriptive R2 values directly. Models 1 and 2 describe
# different responses and their R2 values are not a like-for-like comparison.
comparison_pair <- function(total_model, stratum_model, dat, response) {
  stopifnot(nobs(total_model) == nrow(dat), nobs(stratum_model) == nrow(dat))
  data.frame(
    model = c("Whole-column PAI", "Matching stratum PAI"),
    response = response,
    n_observations = c(nobs(total_model), nobs(stratum_model)),
    R2_log = c(summary(total_model)$r.squared, summary(stratum_model)$r.squared),
    adjusted_R2_log = c(summary(total_model)$adj.r.squared,
                        summary(stratum_model)$adj.r.squared))
}
comparison <- dplyr::bind_rows(
  comparison_pair(m_total, m_stratum, dat_strata,
                  "Within-stratum vertical thermal heterogeneity"),
  comparison_pair(m_total_temporal, m_stratum_temporal, dat_strata_temporal,
                  "Within-stratum temporal thermal heterogeneity"))
write.csv(comparison, file.path(table_dir, "stratum_PAI_comparison.csv"), row.names = FALSE)
print(comparison)
saveRDS(all_results, file.path(table_dir, "PAI_models_and_bootstrap_results.rds"))
