#' Find nearest height layer
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


#' Build Tmax stack by height
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


#' Make cross-section plot
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