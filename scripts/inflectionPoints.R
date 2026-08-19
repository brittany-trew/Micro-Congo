#' Inflection points.
#' CHANGE TO BE DONE SITE WIDE
in.path <- paste0(out.data,"pca_clusters/")

head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
pca <- rast(paste0(in.path, sample.name, "_pca_Clust5.tif"))
chm <- rast(paste0(head.path, "chm.tif"))

#' Process raw PAD output:
processPAD <- function(head.path, dzd, chm){
  all.pad <- list.files(paste0(head.path,"pad"), pattern = ".tif",full.names = TRUE)
  pad.list <- lapply(all.pad, rast)
  pad.r <- do.call(mosaic, pad.list)
  chm <- resample(chm, pad.r)
  
  # ---- heights ----
  n  <- nlyr(pad.r)
  z0 <- (0:(n - 1)) * dzd # height at the *bottom* of each layer (m), starting at 0
  keep <- which(z0 >= dzd) # keep layers starting at dzd and above
  pad.aboveG <- pad.r[[keep]] # Exclude the ground layer as above
  z0G <- z0[keep] # the corresponding height vector for pad.aboveG (same length as nlyr(pad.aboveG))
  
  chm_max <- as.numeric(global(chm, "max", na.rm = TRUE)[1, 1])  # maximum canopy height in the site (m)
  max_start <- floor(chm_max / dzd) * dzd # snap that max height down to the dzd grid (e.g., 23.87 -> 23.75)
  start_heights <- z0G[z0G <= max_start]  # allowable start heights (h0) up to the tallest canopy, on the dzd grid
  pad.aboveG <- pad.aboveG[[seq_along(start_heights)]]  # truncate PAD layers so they match the number of start heights (keeps indexing consistent)
  # Replace NA in PAD with 0 wherever CHM has data
  pad.aboveG[[1]] <- ifel(
    is.na(pad.aboveG[[1]]) & !is.na(chm),
    0,
    pad.aboveG[[1]]
  )
  return(list(pad.aboveG, start_heights, z0G))
}

pad.processed <- processPAD(head.path, dzd, chm)
pad.aboveG <- pad.processed[[1]]
start_heights <- pad.processed[[2]]

pca <- resample(pca,chm)
pad.aboveG <- resample(pad.aboveG, chm)

# ---- split categories (note: maskvalues masks those values OUT) ----
# If your PCA raster has values 1 and 2 and you want to KEEP one class:
pad_cat1 <- mask(pad.aboveG, pca, maskvalues = c(2, 3, 4, 5))
pad_cat2 <- mask(pad.aboveG, pca, maskvalues = c(1, 3, 4, 5))
pad_cat3 <- mask(pad.aboveG, pca, maskvalues = c(1, 2, 4, 5))
pad_cat4 <- mask(pad.aboveG, pca, maskvalues = c(1, 2, 3, 5))
pad_cat5 <- mask(pad.aboveG, pca, maskvalues = c(1, 2, 3, 4))

# ---- build DF ----
df1 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat1)), function(i) {
  data.frame(Forest.Type = 1, height_m = start_heights[i],
             pad = values(pad_cat1[[i]], mat = FALSE))
  }))
df2 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat2)), function(i) {
  data.frame(Forest.Type = 2, height_m = start_heights[i],
             pad = values(pad_cat2[[i]], mat = FALSE))
  }))
df3 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat3)), function(i) {
  data.frame(Forest.Type = 3, height_m = start_heights[i],
             pad = values(pad_cat3[[i]], mat = FALSE))
  }))
df4 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat4)), function(i) {
  data.frame(Forest.Type = 4, height_m = start_heights[i],
             pad = values(pad_cat4[[i]], mat = FALSE))}))
df5 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat5)), function(i) {
  data.frame(Forest.Type = 5, height_m = start_heights[i],
             pad = values(pad_cat4[[i]], mat = FALSE))
  }))
df6 <- do.call(rbind, lapply(seq_len(nlyr(pad.aboveG)), function(i) {
  data.frame(Forest.Type = 6, height_m = start_heights[i],
             pad = values(pad.aboveG[[i]], mat = FALSE))
  }))

dff <- bind_rows(df1, df2, df3, df4, df5, df6)

dff <- filter(dff, height_m < 30)
dff <- na.omit(dff)

# If a site is too empty, return an informative blank plot rather than erroring
if (nrow(dff) == 0) {
  return(
    ggplot() +
      theme_void() +
      ggtitle(sample.name) +
      annotate("text", x = 0, y = 0, label = "No PAD data after filters")
  )
}

mean_profiles <- aggregate(pad ~ Forest.Type + height_m,
                           dplyr::filter(dff, Forest.Type %in% c(1, 2, 3, 4, 5)),
                           mean, na.rm = TRUE)

mean_all <- aggregate(pad ~ height_m,
                      dplyr::filter(dff, Forest.Type == 6),
                      mean, na.rm = TRUE)

outpath <- paste0(out.data,sample.name,"/veg_zoning/")
dir.create(outpath)
all.types <- unique(dff$Forest.Type)
for(a in 1:length(all.types)){
  
  ft <- all.types[a]
  prof_all <- dplyr::filter(dff, Forest.Type == ft)  # all pixels, all heights
  require(mgcv)
  gam_fit <- gam(pad ~ s(height_m, k = 8, fx = TRUE), data = prof_all)
  z_pred <- seq(min(prof_all$height_m), max(prof_all$height_m), by = 0.05)
  pred <- data.frame(height_m = z_pred)
  pred$lad <- as.numeric(predict(gam_fit, newdata = pred))
  pred$lad <- pmax(pred$lad, 0)
  
  plot(pred$lad, pred$height_m, type = "l", col = "darkgreen", lwd = 2,
       xlab = "PAD", ylab = "Height (m)")
  
  lad <- pred$lad
  d1_sign <- sign(diff(lad))
  d1_sign[d1_sign == 0] <- NA
  d1_sign <- zoo::na.locf(d1_sign, na.rm = FALSE)
  sign_change <- diff(d1_sign)
  peak_idx   <- which(sign_change == -2) + 1
  trough_idx <- which(sign_change ==  2) + 1
  
  # Handle boundary peak - if profile starts descending, first point is a peak
  boundary_peak_idx <- if (lad[1] > lad[2]) 1 else NULL
  peak_idx_all <- sort(unique(c(boundary_peak_idx, peak_idx)))
  pred$height_m[peak_idx_all]
  pred$height_m[trough_idx]
  
  # Interleave peaks and troughs in height order
  all_boundaries <- sort(c(pred$height_m[peak_idx_all], pred$height_m[trough_idx]))
  
  boundaries <- zoo::rollmean(all_boundaries, 2)  # midpoint between each adjacent pair
  zone_types <- ifelse(seq_along(all_boundaries) %% 2 == 1, "vegetation", "open_gap")
  
  zones <- data.frame(
    centre_height = round(all_boundaries),
    type          = zone_types,
    zone_start    = round(c(min(pred$height_m), boundaries)),
    zone_end      = round(c(boundaries, max(pred$height_m)))
  )
  zones
  
  saveRDS(zones,paste0(outpath,"VegZoning_FT_",ft,".RDS"))
  
}
print("done")

