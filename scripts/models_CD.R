scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

# Climate distinctness and warming: run after creating `slopes`.
# Required columns: estimate (degrees C/year), dist, type, zone, site, cell.
# Optional coordinates: x_km and y_km, in site-local kilometres.
# Otherwise reconstruct coordinates using cell IDs from the site's PCA raster.
# Cell IDs MUST refer to the same PCA grid used during extraction.
# All models use the same eligible cell-stratum observations.
# Predictor: centred log10(distinctness). Response: untransformed warming slope.
# Uncertainty: 4,999 spatial block bootstrap replicates within sites,
# using 300 m and 500 m blocks; all strata at a cell stay together.
# No ordinary lm p-values/F-tests are used for spatial inference.

block_sizes_m <- c(300, 500)
bootstrap_B <- 4999L
bootstrap_seed <- 123L
variogram_max_cells <- 2000L
run_variograms <- TRUE
period_years <- 2024 - 2004  # 20 elapsed years, despite 21 annual observations


## Prepare one common dataset ----
required <- c("estimate", "dist", "type", "zone", "site", "cell")
missing_cols <- setdiff(required, names(slopes))
if (length(missing_cols)) {
  stop("Missing slopes columns: ", paste(missing_cols, collapse = ", "))
}

slopes <- readRDS( paste0(out.data,"climate_distinctness/slopes_data.RDS"))

dat <- as.data.frame(slopes)
keep <- is.finite(dat$estimate) & is.finite(dat$dist) & dat$dist > 0 &
  !is.na(dat$type) & !is.na(dat$zone) & !is.na(dat$site) &
  is.finite(dat$cell)
message("Excluded ", sum(!keep),
        " rows with missing/invalid inputs or non-positive distinctness.")
dat <- dat[keep, , drop = FALSE]
if (nrow(dat) < 3L) stop("Not enough eligible observations")

if (!all(c("x_km", "y_km") %in% names(dat))) {
  dat$x_km <- dat$y_km <- NA_real_
  for (s in unique(as.character(dat$site))) {
    ii <- which(as.character(dat$site) == s)
    raster_file <- file.path(out.data, "pca_clusters", paste0(s, "_pca_C4.tif"))
    if (!file.exists(raster_file)) {
      stop("Need x_km/y_km columns or this PCA raster: ", raster_file)
    }
    r <- terra::rast(raster_file)
    if (terra::is.lonlat(r)) stop("PCA raster must use projected metre coordinates")
    # The project's PCA rasters are UTM, with coordinate units in metres.
    cells <- dat$cell[ii]
    if (any(cells != floor(cells) | cells < 1 | cells > terra::ncell(r))) {
      stop("Invalid PCA cell IDs for ", s)
    }
    xy <- terra::xyFromCell(r, cells)
    dat$x_km[ii] <- (xy[, 1] - terra::xmin(r)) / 1000
    dat$y_km[ii] <- (xy[, 2] - terra::ymin(r)) / 1000
  }
}
if (any(!is.finite(dat$x_km) | !is.finite(dat$y_km))) {
  stop("Missing/non-finite spatial coordinates: fix before bootstrapping")
}

dat <- dat %>%
  mutate(
    site = droplevels(factor(site)),
    type = droplevels(factor(type)),
    zone = droplevels(factor(zone)),
    type_zone = interaction(type, zone, drop = TRUE),
    dist_log10 = log10(dist),
    dist_c = dist_log10 - mean(dist_log10)
  ) %>%
  arrange(site, cell, type_zone)

if (nlevels(dat$type) < 2L || nlevels(dat$zone) < 2L) {
  stop("The requested model set requires multiple classes and strata")
}
if (anyDuplicated(dat[c("site", "cell", "type_zone")])) {
  stop("Duplicate site-cell-stratum observations: check slopes construction")
}
locations <- distinct(dat, site, cell, x_km, y_km)
if (anyDuplicated(locations[c("site", "cell")])) {
  stop("A cell has inconsistent coordinates across strata")
}
if (anyDuplicated(locations[c("site", "x_km", "y_km")])) {
  stop("Different cell IDs share coordinates within a site")
}
class_variation <- dat %>%
  group_by(type) %>%
  summarise(n_distinct_x = n_distinct(dist_c), .groups = "drop")
if (any(class_variation$n_distinct_x < 2L)) {
  stop("Distinctness does not vary within at least one structural class")
}

write.csv(
  data.frame(n_input = nrow(slopes), n_used = nrow(dat),
             n_excluded = nrow(slopes) - nrow(dat),
             mean_log10_distinctness = mean(dat$dist_log10)),
  file.path(out.data, "tables/tmax_model_prep.csv"), row.names = FALSE
)

## Fit the original five models, using log10 distinctness ----
m1 <- lm(estimate ~ dist_c, data = dat, na.action = na.fail)
m2 <- lm(estimate ~ dist_c * type, data = dat, na.action = na.fail)
m3 <- lm(estimate ~ dist_c * type + zone, data = dat, na.action = na.fail)
m_zone <- lm(estimate ~ type_zone, data = dat, na.action = na.fail)
m_final <- lm(estimate ~ type_zone + dist_c:type,
              data = dat, na.action = na.fail)

models <- list(m1 = m1, m2 = m2, m3 = m3,
               m_zone = m_zone, m_final = m_final)
stopifnot(all(vapply(models, nobs, numeric(1)) == nrow(dat)))

model_descriptions <- c(
  "Log10 climate distinctness only",
  "Log10 climate distinctness x structural class",
  "Log10 climate distinctness x structural class + stratum interval",
  "Class-specific stratum identity only",
  "Class-specific stratum identity + log10 distinctness slopes by class"
)

model_comp <- bind_rows(lapply(seq_along(models), function(i) {
  m <- models[[i]]
  sm <- summary(m)
  data.frame(
    model = names(models)[i], description = model_descriptions[i],
    formula = paste(deparse(formula(m)), collapse = " "),
    n_observations = nobs(m), n_horizontal_cells = nrow(locations),
    n_sites = nlevels(dat$site), R2 = sm$r.squared,
    adjusted_R2 = sm$adj.r.squared, rss = deviance(m),
    df_resid = df.residual(m), AIC = AIC(m),
    aliased_coefficients = sum(is.na(coef(m)))
  )
}))
model_comp$delta_AIC <- model_comp$AIC - min(model_comp$AIC)
model_comp$compared_with <- c(NA_character_, "m1", "m2", NA_character_, "m_zone")
model_comp$delta_R2 <- c(
  NA_real_,
  model_comp$R2[2] - model_comp$R2[1],
  model_comp$R2[3] - model_comp$R2[2],
  NA_real_,
  model_comp$R2[5] - model_comp$R2[4]
)
write.csv(model_comp, file.path(out.data, "tables/model_comparison.csv"),
          row.names = FALSE)

# Same layout as the previous manuscript table. delta_R2 is the increase
# over m1 for m2, over m2 for m3, and over m_zone for m_final.
model_comp_table <- model_comp %>%
  select(model, description, adjusted_R2, delta_R2, rss, df_resid, AIC) %>%
  mutate(across(c(adjusted_R2, delta_R2), ~ round(.x, 3)),
         rss = round(rss, 2), AIC = round(AIC, 1))
write.csv(model_comp_table,
          file.path(out.data, "tables/CD_warming_model_comparison_table.csv"),
          row.names = FALSE)

# These are descriptive changes in explained variance, not spatial F-tests.
# Explicit pairs avoid treating m_zone as the next step after m3.
comparison_pairs <- data.frame(
  reduced = c("m1", "m2", "m_zone"),
  expanded = c("m2", "m3", "m_final")
)
comparison_pairs$delta_R2 <- vapply(seq_len(nrow(comparison_pairs)), function(i) {
  summary(models[[comparison_pairs$expanded[i]]])$r.squared -
    summary(models[[comparison_pairs$reduced[i]]])$r.squared
}, numeric(1))
write.csv(comparison_pairs, file.path(out.data, "/tables/model_comparison_pairs.csv"),
          row.names = FALSE)
print(model_comp)

## Residual diagnostics, consistent with the PAI analysis ----
save_diagnostics <- function(model, id) {
  residual <- rstandard(model)
  png(file.path(paste0(out.data,"plots/"), paste0(id, "_QQ.png")),
      width = 6, height = 6, units = "in", res = 300)
  tryCatch({
    qqnorm(residual, pch = 16, cex = 0.3, main = id)
    qqline(residual, col = "#852e47", lwd = 2)
  }, finally = dev.off())
  png(file.path(paste0(out.data,"plots/"), paste0(id, "_residuals.png")),
      width = 6, height = 6, units = "in", res = 300)
  tryCatch({
    plot(fitted(model), residual, pch = 16, cex = 0.3,
         col = adjustcolor("black", alpha.f = 0.1),
         xlab = "Fitted Tmax warming (degrees C/year)",
         ylab = "Standardised residual", main = id)
    abline(h = 0, lty = 2)
  }, finally = dev.off())
}

save_variograms <- function(model, id) {
  residual <- as.numeric(residuals(model))
  groups <- split(seq_len(nrow(dat)),
                  interaction(dat$site, dat$type_zone, drop = TRUE))
  set.seed(bootstrap_seed)
  for (nm in names(groups)) {
    ii <- groups[[nm]]
    ii <- ii[sample.int(length(ii), min(variogram_max_cells, length(ii)))]
    if (length(ii) < 3L) next
    v <- dat[ii, c("x_km", "y_km")]
    v$.residual <- residual[ii]
    vg <- gstat::variogram(.residual ~ 1, locations = ~ x_km + y_km,
                           data = v, cutoff = 0.5, width = 0.01)
    filename <- paste0(id, "_variogram_", gsub("[^A-Za-z0-9_.-]", "_", nm))
    write.csv(as.data.frame(vg), file.path(paste0(out.data,"tables/"), paste0(filename, ".csv")),
              row.names = FALSE)
    if (nrow(vg) == 0L) next
    png(file.path(paste0(out.data,"plots/"), paste0(filename, ".png")),
        width = 6, height = 6, units = "in", res = 300)
    tryCatch({
      plot(vg$dist, vg$gamma, type = "b", pch = 1,
           xlab = "Distance between cells (km)",
           ylab = "Residual semivariance", main = nm)
    }, finally = dev.off())
  }
}

for (id in names(models)) {
  message("Diagnostics: ", id)
  save_diagnostics(models[[id]], id)
  if (run_variograms) save_variograms(models[[id]], id)
}

## Small design matrices for exact block-summary regression ----
# Each row's design is A[group, ] + dist_c * B[group, ].
# Block/group sums reproduce OLS coefficients without refitting millions of rows.
# No observations are sampled out of the fitted models or the bootstrap.
group_data <- dat[match(levels(dat$type_zone), dat$type_zone),
                  c("type", "zone", "type_zone", "dist_c")]
group_data$dist_c <- 0

prepare_design <- function(model, id) {
  tt <- delete.response(terms(model))
  A <- model.matrix(tt, group_data, contrasts.arg = model$contrasts)
  shifted <- group_data
  shifted$dist_c <- 1
  B <- model.matrix(tt, shifted, contrasts.arg = model$contrasts) - A
  # m3 may have redundant class/interval intercept columns. Keep the same
  # independent columns chosen by lm; predictions and slopes are preserved.
  independent <- model$qr$pivot[seq_len(model$rank)]
  A <- A[, independent, drop = FALSE]
  B <- B[, independent, drop = FALSE]
  if (id == "m1") {
    L <- B[1L, , drop = FALSE]
    effect_names <- "Overall"
  } else {
    effect_names <- levels(dat$type)
    L <- B[match(effect_names, as.character(group_data$type)), , drop = FALSE]
  }
  expected <- as.numeric(L %*% coef(model)[independent])
  if (any(rowSums(abs(L)) == 0) || any(!is.finite(expected))) {
    stop("A class-specific distinctness slope cannot be estimated in ", id)
  }
  list(A = A, B = B, L = L, effect_names = effect_names, expected = expected)
}

effect_models <- c("m1", "m2", "m3", "m_final")
designs <- setNames(lapply(effect_models, function(id) {
  prepare_design(models[[id]], id)
}), effect_models)

fit_from_sums <- function(totals, design) {
  z <- matrix(totals, ncol = 5L, byrow = TRUE)
  n <- z[, 1]
  safe_n <- pmax(n, 1)
  mx <- z[, 2] / safe_n
  my <- z[, 3] / safe_n
  vx <- z[, 4] - z[, 2]^2 / safe_n
  cov_xy <- z[, 5] - z[, 2] * z[, 3] / safe_n
  tol <- 1e-10 * pmax(1, abs(z[, 4]))
  if (any(vx < -tol)) stop("Numerical failure in within-group variance")
  vx <- pmax(vx, 0)
  # Two compressed design rows per group reproduce X'X and X'y exactly:
  # a weighted group mean, and a within-group predictor deviation.
  X <- rbind(
    (design$A + design$B * mx) * sqrt(n),
    design$B * sqrt(vx)
  )
  within_y <- numeric(length(n))
  good <- vx > 0
  within_y[good] <- cov_xy[good] / sqrt(vx[good])
  y <- c(my * sqrt(n), within_y)
  fitted <- lm.fit(X, y)
  if (fitted$rank < ncol(X) || any(!is.finite(fitted$coefficients))) {
    return(rep(NA_real_, nrow(design$L)))
  }
  as.numeric(design$L %*% fitted$coefficients)
}

bootstrap_models <- function(block_size_m) {
  message("Spatial bootstrap: ", block_size_m, " m blocks")
  d <- dat %>%
    transmute(site, cell, type_zone, dist_c, estimate,
              block_x = floor(x_km * 1000 / block_size_m),
              block_y = floor(y_km * 1000 / block_size_m))
  blocks <- distinct(d, site, block_x, block_y) %>%
    arrange(site, block_x, block_y) %>%
    mutate(block_id = row_number())
  d <- left_join(d, blocks, by = c("site", "block_x", "block_y"))
  by_site <- split(blocks$block_id, blocks$site, drop = TRUE)
  if (any(lengths(by_site) < 2L)) {
    stop("Need at least two occupied ", block_size_m, " m blocks per site")
  }
  
  stats <- d %>%
    group_by(block_id, type_zone) %>%
    summarise(n = n(), sx = sum(dist_c), sy = sum(estimate),
              sxx = sum(dist_c^2), sxy = sum(dist_c * estimate),
              .groups = "drop")
  values <- matrix(0, nrow = nrow(blocks), ncol = 5L * nlevels(dat$type_zone))
  gi <- match(as.character(stats$type_zone), levels(dat$type_zone))
  stat_names <- c("n", "sx", "sy", "sxx", "sxy")
  for (k in seq_along(stat_names)) {
    values[cbind(stats$block_id, (gi - 1L) * 5L + k)] <- stats[[stat_names[k]]]
  }
  totals <- colSums(values)
  for (id in effect_models) {
    reproduced <- fit_from_sums(totals, designs[[id]])
    if (!isTRUE(all.equal(reproduced, designs[[id]]$expected, tolerance = 1e-7))) {
      stop("Block-summary coefficients disagree with lm for ", id)
    }
  }
  
  draws <- lapply(designs, function(design) {
    matrix(NA_real_, nrow = bootstrap_B, ncol = length(design$effect_names),
           dimnames = list(NULL, design$effect_names))
  })
  set.seed(bootstrap_seed)
  for (b in seq_len(bootstrap_B)) {
    selected <- unlist(lapply(by_site, function(ids) {
      ids[sample.int(length(ids), length(ids), replace = TRUE)]
    }), use.names = FALSE)
    totals <- colSums(values[selected, , drop = FALSE])
    for (id in effect_models) {
      draws[[id]][b, ] <- fit_from_sums(totals, designs[[id]])
    }
    if (b %% 500L == 0L) message("  Replicate ", b, "/", bootstrap_B)
  }
  # Do not silently discard replicates in which rare identities disappear.
  failures <- vapply(draws, function(x) sum(!apply(is.finite(x), 1, all)), integer(1))
  if (any(failures > 0L)) {
    saveRDS(list(draws = draws, failures = failures),
            file.path(table_dir, paste0("failed_bootstrap_", block_size_m, "m.rds")))
    stop("Unestimable bootstrap replicates: ",
         paste(names(failures)[failures > 0], failures[failures > 0], collapse = "; "),
         ". Inspect class/stratum coverage; replicates have been saved.")
  }
  
  results <- bind_rows(lapply(effect_models, function(id) {
    ci <- apply(draws[[id]], 2, quantile, probs = c(0.025, 0.975))
    sm <- summary(models[[id]])
    data.frame(
      model = id, structural_class = designs[[id]]$effect_names,
      formula = paste(deparse(formula(models[[id]])), collapse = " "),
      predictor = "Centred log10 climate distinctness",
      response = "Tmax warming rate (degrees C/year)",
      n_observations = nrow(dat), n_horizontal_cells = nrow(locations),
      n_sites = nlevels(dat$site), R2 = sm$r.squared,
      adjusted_R2 = sm$adj.r.squared,
      effect_C_per_year = designs[[id]]$expected,
      effect_C_per_decade = designs[[id]]$expected * 10,
      conf_low_decade = ci[1, ] * 10,
      conf_high_decade = ci[2, ] * 10,
      effect_basis = "Per tenfold increase in climate distinctness",
      block_size_m = block_size_m, bootstrap_replicates = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      CI_method = "95% percentile spatial-block bootstrap within sites"
    )
  }))
  counts <- d %>%
    group_by(site) %>%
    summarise(n_observations = n(), n_horizontal_cells = n_distinct(cell),
              n_blocks = n_distinct(block_id), .groups = "drop") %>%
    mutate(block_size_m = block_size_m)
  list(results = results, draws = draws, block_counts = counts)
}

## Bootstrap confidence intervals and save results ----
bootstraps <- list()
for (size in block_sizes_m) {
  bootstraps[[as.character(size)]] <- bootstrap_models(size)
  saveRDS(bootstraps[[as.character(size)]],
          file.path(paste0(out.data,"tables/"), paste0("bootstrap_", size, "m.rds")))
}
all_effects <- bind_rows(lapply(bootstraps, function(x) x$results))
block_counts <- bind_rows(lapply(bootstraps, function(x) x$block_counts))
write.csv(all_effects, file.path(paste0(out.data,"tables/"), "all_model_distinctness_effects.csv"),
          row.names = FALSE)
write.csv(block_counts, file.path(paste0(out.data,"tables/"), "bootstrap_block_counts.csv"),
          row.names = FALSE)

dist_effects <- all_effects %>% dplyr::filter(model == "m_final")
write.csv(dist_effects, file.path(paste0(out.data,"tables/"), "dist_effects.csv"), row.names = FALSE)

supplementary_table <- all_effects %>%
  select(model, formula, structural_class, n_observations, n_horizontal_cells,
         n_sites, R2, adjusted_R2, effect_C_per_decade,
         conf_low_decade, conf_high_decade, effect_basis,
         block_size_m, bootstrap_replicates, CI_method) %>%
  mutate(across(c(R2, adjusted_R2, effect_C_per_decade,
                  conf_low_decade, conf_high_decade), ~ round(.x, 3)))
write.csv(supplementary_table,
          file.path(paste0(out.data,"tables/"), "CD_warming_supplementary_table.csv"),
          row.names = FALSE)

## Effects across the observed log10 distinctness distribution ----
# These spans are fixed descriptive features of the observed data.
# Their own sampling uncertainty is not bootstrapped.
iqr_by_type <- dat %>%
  group_by(type) %>%
  summarise(
    log10_dist_iqr = IQR(dist_log10),
    log10_dist_range = diff(range(dist_log10)),
    .groups = "drop"
  ) %>%
  mutate(structural_class = as.character(type)) %>%
  select(-type)

dist_effects_scaled <- dist_effects %>%
  left_join(iqr_by_type, by = "structural_class") %>%
  mutate(
    effect_IQR_decade = effect_C_per_decade * log10_dist_iqr,
    lower_IQR_decade = conf_low_decade * log10_dist_iqr,
    upper_IQR_decade = conf_high_decade * log10_dist_iqr,
    effect_IQR_20yrs = effect_IQR_decade * period_years / 10,
    lower_IQR_20yrs = lower_IQR_decade * period_years / 10,
    upper_IQR_20yrs = upper_IQR_decade * period_years / 10,
    effect_range_decade = effect_C_per_decade * log10_dist_range,
    lower_range_decade = conf_low_decade * log10_dist_range,
    upper_range_decade = conf_high_decade * log10_dist_range,
    effect_range_20yrs = effect_range_decade * period_years / 10,
    lower_range_20yrs = lower_range_decade * period_years / 10,
    upper_range_20yrs = upper_range_decade * period_years / 10
  )
write.csv(dist_effects_scaled, file.path(paste0(out.data,"tables/"), "dist_effects_scaled.csv"),
          row.names = FALSE)

## Plot final class-specific effects for both block sizes ----
class_colours <- c("1" = "#ca763b", "2" = "#852e47",
                   "3" = "#333333", "4" = "#a09c33")
dp <- ggplot(
  dist_effects,
  aes(y = factor(structural_class, levels = levels(dat$type)),
      x = effect_C_per_decade, xmin = conf_low_decade,
      xmax = conf_high_decade, colour = structural_class)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_pointrange(orientation = "y", linewidth = 0.5) +
  scale_colour_manual(values = class_colours) +
  facet_wrap(~ block_size_m, labeller = label_both) +
  labs(
    x = paste0("Change in Tmax warming rate (°C decade⁻¹)\n",
               "per tenfold increase in climate distinctness"),
    y = "Structural class", colour = "Structural class"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

ggsave(file.path(paste0(out.data,"plots/"), "CDWarming_effectSize.png"),
       plot = dp, width = 10, height = 6, units = "in", dpi = 300)

saveRDS(
  list(models = models, model_comparison = model_comp,
       effects = all_effects, block_counts = block_counts,
       scaled_effects = dist_effects_scaled,
       centring_mean_log10_distinctness = mean(dat$dist_log10),
       predictor = "log10(distinctness)", response_units = "degrees C/year",
       block_sizes_m = block_sizes_m, bootstrap_B = bootstrap_B,
       bootstrap_seed = bootstrap_seed),
  file.path(paste0(out.data,"tables/"), "CD_warming_models_and_results.rds")
)
print(supplementary_table)
print(dist_effects_scaled)
