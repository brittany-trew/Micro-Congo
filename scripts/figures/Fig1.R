#' Point cloud cross-section shaded by Tmax
#' ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

library(plotly)
library(signal)

las.path <- paste0(lidar.data,site.location,"/",sample.name,"/",year)
las.options <- list.files(
  paste0(las.path),
  full.names = TRUE
)

tmaxpath <- paste0(out.data, "MeanDailyTemps/")

#' ---------------------------------------------------------------------
# Read and normalise LAS
#' ---------------------------------------------------------------------
las2d <- readLAS(las.options[[1]])

dtm <- rasterize_terrain(
  las2d,
  res = 5,
  algorithm = tin()
)
las_norm <- normalize_height(las2d, dtm)

#' Extract the normalized point attributes as a data.table
las_dt <- as.data.table(las_norm@data)
#' Name & store normalized height above ground
las_dt[, z_h := Z]
#' Remove points below the estimated ground surface
las_dt <- las_dt[z_h >= 0]
#' Convert the projected X and Y coordinates to local coordinates
las_dt[, `:=`(
  norm_X = X - min(X, na.rm = TRUE),
  norm_Y = Y - min(Y, na.rm = TRUE)
)]

#' ---------------------------------------------------------------------
# Build Tmax raster stack
#' ---------------------------------------------------------------------
tmax_obj <- build_tmax_stack(tmaxpath)
r <- tmax_obj$r
heights <- tmax_obj$heights
height_check <- tmax_obj$height_check


# ============================================================
# Colour limits
# ============================================================

# Option 2: estimate limits from raster sample instead
r_sample <- terra::spatSample(
  r,
  size = 200000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)
global_lims <- quantile(unlist(r_sample), c(0.05, 0.95), na.rm = TRUE)


# ============================================================
# Make several cross-sections
# ============================================================

slice_ys <- as.numeric(quantile(
  las_dt$Y,
  probs = c(0.1, 0.3, 0.5, 0.7, 0.9),
  na.rm = TRUE
))

plots <- lapply(seq_along(slice_ys), function(i) {
  
  make_cross_section(
    las_dt = las_dt,
    r = r,
    heights = heights,
    slice_y = slice_ys[i],
    slice_width = 25,
    n_points = 2000000,
    lims = global_lims,
    label = paste0("section ", i),
    point_size = 0.01,
    point_alpha = 0.45
  )
  
})

names(plots) <- paste0("section_", seq_along(plots))


# View plots
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]


# ============================================================
# Save several cross-sections
# ============================================================

for (i in seq_along(plots)) {
  
  ggsave(
    filename = paste0(out.data, "plots/crossSection_tmaxV2_section_", i, ".png"),
    plot = plots[[i]],
    width = 8,
    height = 2,
    units = "in",
    dpi = 300
  )
  
}
