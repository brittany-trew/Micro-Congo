#' Point cloud cross-section shaded by Tmax
#' ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

las.options <- list.files(
  paste0(las.path),
  full.names = TRUE
)

tmaxpath <- paste0(out.data, sample.name,"/MeanDailyTemps/")

#' ---------------------------------------------------------------------
#' Read and normalise LAS
mem.maxVSize(vsize = 180000)
rm(las_dt)
las2d <- readLAS(las.options[[4]])
dtm <- rasterize_terrain(
  las2d,
  res = 5,
  algorithm = tin()
)
las_norm <- normalize_height(las2d, dtm)
rm(las2d)
rm(dtm)
gc()

#' Extract the normalized point attributes as a data.table
las_dt <- as.data.table(las_norm@data)
rm(las_norm)
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
#' Build Tmax raster stack
tmax_obj <- build_tmax_stack(tmaxpath)
r <- tmax_obj$r
heights <- tmax_obj$heights
height_check <- tmax_obj$height_check

# ============================================================
# Make several cross-sections
# ============================================================
slice_ys <- as.numeric(quantile(
  las_dt$Y,
  probs = c(0.1,0.2,0.5),
  na.rm = TRUE
))

p <- make_cross_section(
  las_dt = las_dt,
  r = r,
  heights = heights,
  slice_y = slice_ys[1],
  slice_width = 15,
  n_points = 10000000,
  point_size = 0.0001,
  point_alpha = 0.7,
  lims = c(28,36)
)

plt <- p + coord_fixed(ratio = 1)
plt

ggsave(filename = paste0(
  out.data, "plots/crossSection_tmax_section_", sample.name, ".png"),
  plot = plt,
  width = 8,
  height = 3,
  units = "in",
  dpi = 400
  )
