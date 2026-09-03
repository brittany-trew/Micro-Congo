#' Vertical structural stratification from PAD profiles
in.path <- paste0(out.data,"pca_clusters/")
outpath <- paste0(out.data,"strata/")
dir.create(outpath, showWarnings = F)


#' Calculate Max Heights.
max_heights <- c()
for(i in 1:length(all.samples)){
  
  sample.name <- all.samples[i]
  head.path <- paste0(
    mclidar.out,
    site.location,
    "/",
    sample.name,
    "/",
    year,
    "/"
  )
  
  chm <- rast(paste0(head.path, "chm.tif"))
  
  max_heights[i] <- global(
    chm,
    "max",
    na.rm = TRUE
  )[1, 1]
}

study_max_height <- max(max_heights)
study_max_height


#' Extract PAD profiles.
df_list <- list()
for(i in 1:length(all.samples)){
  
  #' Set datastream.
  sample.name <- all.samples[i]
  head.path <- paste0(mclidar.out,site.location,"/",sample.name,"/",year,"/")
  chm <- rast(paste0(head.path,"chm.tif"))
  
  all.pad <- list.files(paste0(head.path,"pad"), pattern = ".tif", full.names = T)
  pad.r <- lapply(all.pad, rast)
  pad0 <- processPAD(pad.r, chm)
  
  start_heights <- seq(
    from = 0,
    by = dzd,
    length.out = nlyr(pad0)
  )
  
  pca <- rast(paste0(in.path, sample.name, "_pca_C4.tif"))
  
  #' Use maskvalues to mask values OUT.
  pad_cat1 <- mask(pad0, pca, maskvalues = c(2, 3, 4))
  pad_cat2 <- mask(pad0, pca, maskvalues = c(1, 3, 4))
  pad_cat3 <- mask(pad0, pca, maskvalues = c(1, 2, 4))
  pad_cat4 <- mask(pad0, pca, maskvalues = c(1, 2, 3))
  
  df1 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat1)), function(j) {
    data.frame(Forest.Type = 1, height_m = start_heights[j],
               pad = values(pad_cat1[[j]], mat = FALSE))
  }))
  
  df2 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat2)), function(j) {
    data.frame(Forest.Type = 2, height_m = start_heights[j],
               pad = values(pad_cat2[[j]], mat = FALSE))
  }))
  
  df3 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat3)), function(j) {
    data.frame(Forest.Type = 3, height_m = start_heights[j],
               pad = values(pad_cat3[[j]], mat = FALSE))
  }))
  
  df4 <- do.call(rbind, lapply(seq_len(nlyr(pad_cat4)), function(j) {
    data.frame(Forest.Type = 4, height_m = start_heights[j],
               pad = values(pad_cat4[[j]], mat = FALSE))
  }))
  
  df_all <- do.call(rbind, lapply(seq_len(nlyr(pad0)), function(j) {
    data.frame(Forest.Type = 0, height_m = start_heights[j],
               pad = values(pad0[[j]], mat = FALSE))
  }))
  
  dff <- bind_rows(df1, df2, df3, df4, df_all)
  
  dff <- dff %>%
    dplyr::filter(
      !is.na(pad),
      height_m <= study_max_height
    )
  
  df_list[[sample.name]] <- dff
}

df <- bind_rows(df_list)


#' Calculate mean vertical PAD profiles.
mean_profiles <- aggregate(
  pad ~ Forest.Type + height_m,
  dplyr::filter(df, Forest.Type %in% c(1, 2, 3, 4)),
  mean,
  na.rm = TRUE
)

write_csv(mean_profiles, paste0(outpath,"mean_profiles.csv"))


#' Overall forest profile.
mean_all <- aggregate(
  pad ~ height_m,
  dplyr::filter(df, Forest.Type == 0),
  mean,
  na.rm = TRUE
)

write_csv(mean_all, paste0(outpath,"mean_profiles_all.csv"))


#' Parameters for structural zoning.
prominence_threshold <- 0.2

#' The ground boundary is defined where the initial decline
#' falls below this proportion of its maximum gradient.
ground_gradient_fraction <- 0.1


#' Fit GAMs to the mean profiles.
all.types <- sort(unique(mean_profiles$Forest.Type))

for(a in 1:length(all.types)){
  
  ft <- all.types[a]
  
  prof_all <- dplyr::filter(
    mean_profiles,
    Forest.Type == ft
  )
  
  prof_all <- dplyr::arrange(
    prof_all,
    height_m
  )
  
  
  #' Square-root transform PAD so that the large ground-level
  #' values do not dominate the fitted vertical profile.
  prof_all$sqrt_pad <- sqrt(prof_all$pad)
  
  
  #' Fit penalised GAM.
  gam_fit <- mgcv::gam(
    sqrt_pad ~ s(height_m, k = 15),
    data = prof_all,
    method = "REML"
  )
  
  
  #' Evaluate the continuous fitted curve every 0.01 m.
  prediction_heights <- seq(
    min(prof_all$height_m),
    max(prof_all$height_m),
    by = 0.01
  )
  
  pred <- data.frame(height_m = prediction_heights)
  
  pred$sqrt_pad <- as.numeric(
    predict(gam_fit, newdata = pred)
  )
  
  
  #' Negative square-root PAD values are physically impossible,
  #' so constrain them to zero before back-transforming.
  pred$sqrt_pad_pos <- pmax(pred$sqrt_pad, 0)
  pred$pad <- pred$sqrt_pad_pos^2
  
  
  #' Diagnostic plot.
  plot(pred$pad, pred$height_m, type = "l", col = "darkgreen", lwd = 2,
       xlab = "PAD", ylab = "Height (m)",
       main = paste("Structural Class", ft))
  
  pad <- pred$pad
  
  
  #' ------------------------------------------------------------
  #' STEP 1: Identify local maxima and minima.
  #' ------------------------------------------------------------
  
  #' Calculate direction of change between consecutive points.
  d1_sign <- sign(diff(pad))
  
  #' Replace flat sections with NA.
  d1_sign[d1_sign == 0] <- NA
  
  #' Fill internal flat sections.
  d1_sign <- zoo::na.locf(
    d1_sign,
    na.rm = FALSE
  )
  
  #' Fill any remaining leading NAs.
  d1_sign <- zoo::na.locf(
    d1_sign,
    fromLast = TRUE,
    na.rm = FALSE
  )
  
  #' Calculate changes in direction.
  sign_change <- diff(d1_sign)
  
  #' Increasing to decreasing = local maximum.
  peak_idx <- which(sign_change == -2) + 1
  
  #' Decreasing to increasing = local minimum.
  trough_idx <- which(sign_change == 2) + 1
  
  #' If PAD immediately declines from ground level,
  #' treat the lower endpoint as a density maximum.
  boundary_peak_idx <- if(
    pad[1] > pad[2]
  ) 1 else NULL
  
  peak_idx <- sort(unique(c(
    boundary_peak_idx,
    peak_idx
  )))
  
  
  #' Combine maxima and minima.
  extrema <- data.frame(
    idx = c(peak_idx, trough_idx),
    height_m = c(
      pred$height_m[peak_idx],
      pred$height_m[trough_idx]
    ),
    type = c(
      rep("dense", length(peak_idx)),
      rep("open", length(trough_idx))
    )
  )
  
  #' Sort from ground upwards.
  extrema <- extrema[
    order(extrema$height_m),
  ]
  
  
  #' PAD associated with each fitted extremum.
  extrema$pad <- pred$pad[extrema$idx]
  
  
  #' Difference from preceding extremum.
  extrema$contrast_prev <- c(
    NA,
    abs(diff(extrema$pad))
  )
  
  
  #' ------------------------------------------------------------
  #' STEP 2: Calculate prominence of interior extrema.
  #' ------------------------------------------------------------
  
  #' Prominence is the smaller PAD contrast between an extremum
  #' and its two neighbouring extrema.
  extrema$prominence <- NA_real_
  
  if(nrow(extrema) > 2){
    
    for(i in 2:(nrow(extrema) - 1)){
      
      extrema$prominence[i] <- min(
        abs(extrema$pad[i] - extrema$pad[i - 1]),
        abs(extrema$pad[i] - extrema$pad[i + 1])
      )
    }
  }
  
  
  #' Maximum observed above-ground PAD provides the reference
  #' against which prominence is expressed.
  reference_pad <- max(
    prof_all$pad[
      prof_all$height_m >= dzd
    ],
    na.rm = TRUE
  )
  
  extrema$prominence_pct <-
    100 * extrema$prominence / reference_pad
  
  print(extrema)
  
  
  #' ------------------------------------------------------------
  #' STEP 3: Retain structurally prominent extrema.
  #' ------------------------------------------------------------
  
  extrema_keep <- extrema %>%
    dplyr::filter(
      !is.na(prominence_pct),
      prominence_pct >= prominence_threshold
    )
  
  
  #' Ground is retained as the lower structural anchor.
  lower_anchor <- data.frame(
    idx = 1,
    height_m = min(pred$height_m),
    type = "dense",
    pad = pred$pad[1]
  )
  
  
  strata_extrema <- dplyr::bind_rows(
    lower_anchor,
    extrema_keep %>%
      dplyr::select(idx, height_m, type, pad)
  ) %>%
    dplyr::arrange(height_m)
  
  
  #' If the final retained environment is dense, add an upper
  #' open anchor so the transition into sparse upper canopy
  #' can be located.
  upper_anchor_added <- FALSE
  
  if(strata_extrema$type[nrow(strata_extrema)] != "open"){
    
    upper_anchor <- data.frame(
      idx = nrow(pred),
      height_m = max(pred$height_m),
      type = "open",
      pad = pred$pad[nrow(pred)]
    )
    
    strata_extrema <- dplyr::bind_rows(
      strata_extrema,
      upper_anchor
    )
    
    upper_anchor_added <- TRUE
  }
  
  
  #' ------------------------------------------------------------
  #' STEP 4: Calculate rate of PAD change with height.
  #' ------------------------------------------------------------
  
  pred$dPAD <- NA_real_
  
  #' Centred finite-difference approximation of first derivative.
  pred$dPAD[2:(nrow(pred) - 1)] <-
    (
      pred$pad[3:nrow(pred)] -
        pred$pad[1:(nrow(pred) - 2)]
    ) /
    (
      pred$height_m[3:nrow(pred)] -
        pred$height_m[1:(nrow(pred) - 2)]
    )
  
  
  #' ------------------------------------------------------------
  #' STEP 5: Locate structural boundaries.
  #' ------------------------------------------------------------
  
  boundary_idx <- c()
  
  
  #' FIRST BOUNDARY:
  #' The ground peak is an endpoint, so using the point of
  #' maximum gradient would place the boundary almost at 0 m.
  #' Instead, find where the initial PAD decline has levelled off.
  
  first_upper_idx <- strata_extrema$idx[2]
  
  #' Search between the ground and first retained environment.
  ground_section <- 2:first_upper_idx
  
  ground_section <- ground_section[
    is.finite(pred$dPAD[ground_section])
  ]
  
  
  #' Identify maximum initial PAD gradient.
  max_ground_gradient_idx <- ground_section[
    which.max(abs(pred$dPAD[ground_section]))
  ]
  
  max_ground_gradient <- abs(
    pred$dPAD[max_ground_gradient_idx]
  )
  
  
  #' Define when the decline has substantially levelled off.
  ground_gradient_threshold <-
    ground_gradient_fraction * max_ground_gradient
  
  
  #' Search only after the point of maximum initial decline.
  ground_search <- max_ground_gradient_idx:first_upper_idx
  
  ground_search <- ground_search[
    is.finite(pred$dPAD[ground_search])
  ]
  
  
  #' Find points where the gradient has fallen below the threshold.
  ground_candidates <- ground_search[
    abs(pred$dPAD[ground_search]) <= ground_gradient_threshold
  ]
  
  
  #' Use the first point where the decline has levelled off.
  if(length(ground_candidates) > 0){
    
    ground_boundary_idx <- ground_candidates[1]
    
  } else {
    
    #' Fallback: use the weakest remaining gradient.
    ground_boundary_idx <- ground_search[
      which.min(abs(pred$dPAD[ground_search]))
    ]
  }
  
  
  boundary_idx <- c(
    boundary_idx,
    ground_boundary_idx
  )
  
  
  #' Remaining boundaries:
  #' Between subsequent retained environments, locate the
  #' point where PAD changes most rapidly.
  if(nrow(strata_extrema) > 2){
    
    for(i in 2:(nrow(strata_extrema) - 1)){
      
      #' Lower and upper retained environments.
      idx1 <- strata_extrema$idx[i]
      idx2 <- strata_extrema$idx[i + 1]
      
      #' Search between the two environments.
      section_idx <- idx1:idx2
      
      #' Remove positions where derivative is unavailable.
      section_idx <- section_idx[
        is.finite(pred$dPAD[section_idx])
      ]
      
      #' Find the point of maximum absolute PAD change.
      change_idx <- section_idx[
        which.max(abs(pred$dPAD[section_idx]))
      ]
      
      boundary_idx <- c(
        boundary_idx,
        change_idx
      )
    }
  }
  
  
  #' Convert boundary indices to heights.
  boundary_heights <- pred$height_m[boundary_idx]
  
  
  #' ------------------------------------------------------------
  #' STEP 6: Convert boundaries into vertical strata.
  #' ------------------------------------------------------------
  
  zones <- data.frame(
    type = strata_extrema$type,
    centre_height = strata_extrema$height_m,
    zone_start = c(
      min(pred$height_m),
      boundary_heights
    ),
    zone_end = c(
      boundary_heights,
      max(pred$height_m)
    )
  )
  
  
  #' Snap final boundaries to native 0.5-m PAD resolution.
  zones <- zones %>%
    dplyr::mutate(
      zone_start = round(zone_start / dzd) * dzd,
      zone_end = round(zone_end / dzd) * dzd
    )
  
  
  #' Artificial upper anchor is not a true extremum.
  if(upper_anchor_added){
    zones$centre_height[nrow(zones)] <- NA
  }
  
  
  cat("\nStructural Class:", ft, "\n")
  print(zones)
  
  
  saveRDS(
    zones,
    paste0(outpath,"StrataZones_ft_",ft,".RDS")
  )
}

print("done")