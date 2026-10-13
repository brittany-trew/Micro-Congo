#' Find nearest height layer
nearest_height_layer <- function(z, heights) {
  # Boundaries halfway between successive model heights
  midpoints <- (head(heights, -1) + tail(heights, -1)) / 2
  findInterval(z, midpoints) + 1L
}


#' Build Tmax stack by height
build_tmax_stack <- function(tmaxpath) {
  
  all.tmax <- list.files(
    tmaxpath,
    full.names = TRUE
  )
  
  if (length(all.tmax) == 0) {
    stop("No Tmax files found in: ", tmaxpath)
  }
  
  # Extract numeric height from each filename
  height_from_file <- basename(all.tmax)
  height_from_file <- tools::file_path_sans_ext(height_from_file)
  height_from_file <- sub("^.*MeanDailyMax_", "", height_from_file)
  height_from_file <- as.numeric(sub("m$", "", height_from_file))
  
  if (anyNA(height_from_file)) {
    stop("Some heights could not be extracted from Tmax filenames.")
  }
  
  # Apply the SAME ordering to both files and heights
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


#' Make cross-section plot
make_cross_section <- function(las_dt,
                               r,
                               heights,
                               slice_y,
                               slice_width = 5,
                               n_points = 100000,
                               lims = NULL,
                               label = NULL,
                               cols = cols_tmax,
                               point_size = 0.3,
                               point_alpha = 0.6) {
  
  # Find points within the strip (total thickness = 2 * slice_width)
  idx <- which(
    las_dt$Y >= slice_y - slice_width &
      las_dt$Y <= slice_y + slice_width &
      las_dt$norm_X >= 0 &
      las_dt$norm_X <= 250
  )
  
  if (length(idx) == 0L) {
    stop("No LiDAR points found in this cross-section.")
  }
  
  # Sample before copying, keeping only the necessary columns
  set.seed(1)
  if (length(idx) > n_points) {
    idx <- idx[sample.int(length(idx), n_points)]
  }
  
  cs_plot <- data.table::copy(
    las_dt[idx, .(X, Y, norm_X, z_h)]
  )
  
  # Match each point to its nearest modelled height
  cs_plot[, layer_id := nearest_height_layer(z_h, heights)]
  cs_plot[, tmax := NA_real_]
  
  # Extract only the matching temperature layer for each point
  for (j in unique(cs_plot$layer_id)) {
    rows <- which(cs_plot$layer_id == j)
    
    vals <- terra::extract(
      r[[j]],
      as.matrix(cs_plot[rows, .(X, Y)])
    )
    
    data.table::set(
      cs_plot,
      i = rows,
      j = "tmax",
      value = vals[[1]]
    )
  }
  
  cs_plot <- cs_plot[is.finite(tmax)]
  
  if (nrow(cs_plot) == 0L) {
    stop("No valid temperatures found in this cross-section.")
  }
  
  if (is.null(lims)) {
    lims <- quantile(cs_plot$tmax, c(0.1, 0.9), na.rm = TRUE)
  }
  
  ggplot(cs_plot, aes(x = norm_X, y = z_h, colour = tmax)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_colour_gradientn(
      name = "Tmax",
      colours = cols,
      limits = lims,
      oob = scales::squish
    ) +
    labs(
      x = "",
      y = ""
    ) +
    theme_classic()
}
