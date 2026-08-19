# ============================================================
# Point cloud cross-section shaded by Tmax
# ============================================================

library(lidR)
library(terra)
library(ggplot2)
library(data.table)
library(scales)

# Optional, if you still need them elsewhere
library(rgl)
library(plotly)
library(signal)
library(tidyr)
library(dplyr)


# ============================================================
# Colour palette
# ============================================================

cols_tmax <- c(
  "#192835",  # dark navy
  "#68AE9A",  # muted teal
  "#faa825",  # sand/peach
  "#d9416b",  # rose pink
  "#8A211B"   # dark red
)



# ============================================================
# Helper function: nearest height layer
# ============================================================

nearest_height_layer <- function(z, heights) {
  
  idx <- findInterval(z, heights)
  
  idx[idx == 0] <- 1
  idx[idx >= length(heights)] <- length(heights) - 1
  
  left  <- heights[idx]
  right <- heights[idx + 1]
  
  out <- ifelse(abs(z - left) <= abs(z - right), idx, idx + 1)
  
  out[z <= min(heights)] <- 1
  out[z >= max(heights)] <- length(heights)
  
  return(out)
}


# ============================================================
# Helper function: build Tmax stack by height
# ============================================================

build_tmax_stack <- function(tmaxpath) {
  
  all.tmax <- list.files(
    tmaxpath,
    pattern = "yearlyMean_tmax.*\\.tif$",
    full.names = TRUE
  )
  
  if (length(all.tmax) == 0) {
    stop("No Tmax files found in: ", tmaxpath)
  }
  
  # Extract numeric height from filename
  height_from_file <- basename(all.tmax)
  height_from_file <- tools::file_path_sans_ext(height_from_file)
  height_from_file <- sub("^.*tmax_", "", height_from_file)
  height_from_file <- as.numeric(height_from_file)
  
  if (anyNA(height_from_file)) {
    stop("Some heights could not be extracted from Tmax filenames.")
  }
  
  # Sort by real numeric height, not alphabetical filename order
  ord <- order(height_from_file)
  
  all.tmax <- all.tmax[ord]
  heights <- height_from_file[ord]
  
  height_check <- data.frame(
    layer = seq_along(all.tmax),
    file = basename(all.tmax),
    height = heights
  )
  
  print(height_check)
  
  # Build raster stack
  r_list <- vector("list", length(all.tmax))
  
  for (i in seq_along(all.tmax)) {
    r1 <- rast(all.tmax[[i]])
    r_list[[i]] <- mean(r1, na.rm = TRUE)
  }
  
  r <- rast(r_list)
  names(r) <- paste0("h_", gsub("\\.", "p", heights))
  
  stopifnot(nlyr(r) == length(heights))
  stopifnot(all(diff(heights) > 0))
  
  return(list(
    r = r,
    heights = heights,
    height_check = height_check
  ))
}


# ============================================================
# Helper function: make cross-section plot
# ============================================================

make_cross_section <- function(las_dt,
                               r,
                               heights,
                               slice_y,
                               slice_width = 25,
                               n_points = 2000000,
                               lims = NULL,
                               label = NULL,
                               cols = cols_tmax,
                               point_size = 0.01,
                               point_alpha = 0.45) {
  
  # Extract cross-section
  cs <- las_dt[
    Y >= slice_y - slice_width &
      Y <= slice_y + slice_width
  ]
  
  if (nrow(cs) == 0) {
    stop("No LiDAR points found in this cross-section.")
  }
  
  # Sample points for plotting
  set.seed(1)
  cs_plot <- copy(cs[sample(.N, min(.N, n_points))])
  
  # Match point height to nearest modelled Tmax height
  cs_plot[, z_match := pmin(pmax(z_h, min(heights)), max(heights))]
  cs_plot[, layer_id := nearest_height_layer(z_match, heights)]
  
  # Extract all Tmax layers at point X/Y locations
  temp_vals <- terra::extract(
    r,
    as.matrix(cs_plot[, .(X, Y)])
  )
  
  # Remove ID column if terra adds one
  if ("ID" %in% names(temp_vals)) {
    temp_vals <- temp_vals[, names(temp_vals) != "ID"]
  }
  
  temp_mat <- as.matrix(temp_vals)
  
  # For each point, select the Tmax layer matching its height
  cs_plot[, tmax := temp_mat[cbind(seq_len(.N), layer_id)]]
  cs_plot <- cs_plot[!is.na(tmax)]
  
  # Define colour limits if not supplied
  if (is.null(lims)) {
    lims <- quantile(cs_plot$tmax, c(0.05, 0.95), na.rm = TRUE)
  }
  
  if (is.null(label)) {
    label <- paste0("Y = ", round(slice_y, 1), " m")
  }
  
  # Plot
  p <- ggplot(cs_plot, aes(x = norm_X, y = z_h, colour = tmax)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_colour_gradientn(
      name = "Tmax",
      colours = cols,
      limits = lims,
      oob = scales::squish
    ) +
    labs(
      x = "Distance across tile (m)",
      y = "Height above ground (m)"
    ) +
    theme_classic()
  
  return(p)
}


# ============================================================
# User paths and site
# ============================================================

las.path <- "/Volumes/HD01/LiDAR/OdzalaCongo/"

site <- all.samples[[4]]

las.options <- list.files(
  paste0(las.path, site, "/2021/"),
  full.names = TRUE
)

tmaxpath <- paste0(volumes3, site, "/monthlyTemps/")


# ============================================================
# Read and normalise LAS
# ============================================================

las2d <- readLAS(las.options[[1]])

dtm <- rasterize_terrain(
  las2d,
  res = 5,
  algorithm = tin()
)

las_norm <- normalize_height(las2d, dtm)

las_dt <- as.data.table(las_norm@data)

las_dt[, z_h := Z]
las_dt <- las_dt[z_h >= 0]

las_dt[, `:=`(
  norm_X = X - min(X, na.rm = TRUE),
  norm_Y = Y - min(Y, na.rm = TRUE)
)]


# ============================================================
# Build Tmax raster stack
# ============================================================

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
