#' Warming of tmax - timeseries
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

yr.seq <- 2004:2024

# Read a height file and align its layers to all 21 years
read_tmax <- function(file, template, years) {
  
  r <- tryCatch(
    rast(file),
    error = function(e) {
      warning("Skipping unreadable file: ", file, call. = FALSE)
      NULL
    }
  )
  
  if (is.null(r)) return(NULL)
  
  # Read actual years from layer names, e.g. "2004" or "X2004"
  layer_years <- suppressWarnings(
    as.integer(sub("^X", "", names(r)))
  )
  
  valid <- !is.na(layer_years) & layer_years %in% years
  
  if (!any(valid)) {
    warning(
      "Skipping file: no recognisable year names in ", file,
      call. = FALSE
    )
    return(NULL)
  }
  
  if (anyDuplicated(layer_years[valid])) {
    warning(
      "Skipping file: duplicate year names in ", file,
      call. = FALSE
    )
    return(NULL)
  }
  
  # Do not silently combine rasters with different grids
  compareGeom(template, r, stopOnError = TRUE)
  
  r[!is.finite(r)] <- NA
  
  empty <- template * NA_real_
  aligned <- rast()
  
  for (yr in years) {
    idx <- which(layer_years == yr)
    
    if (length(idx) == 0L) {
      aligned <- c(aligned, empty)
    } else {
      aligned <- c(aligned, r[[idx]])
    }
  }
  
  names(aligned) <- as.character(years)
  aligned
}


site_list <- list()

outpath <- paste0(out.data, "tmax_warming/")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

strata.path <- paste0(out.data, "/strata/")
pca.path <- paste0(out.data, "pca_clusters/")

for (a in seq_along(all.samples)) {
  
  sample.name <- all.samples[[a]]
  inpath <- paste0(out.data, sample.name, "/yearlyTemps")
  
  pca.file <- paste0(pca.path, sample.name, "_pca_C4.tif")
  
  if (!file.exists(pca.file)) {
    warning("Skipping site: missing PCA raster for ", sample.name,
            call. = FALSE)
    next
  }
  
  pca <- rast(pca.file)
  
  forest.types <- sort(unique(as.vector(values(pca))))
  forest.types <- forest.types[is.finite(forest.types)]
  
  all.tmax <- list.files(
    inpath,
    pattern = "^yearlyMean_tmax_[0-9.]+\\.tif$",
    full.names = TRUE
  )
  
  if (length(all.tmax) == 0L) {
    warning("Skipping site: no Tmax files for ", sample.name,
            call. = FALSE)
    next
  }
  
  heights_found <- as.numeric(
    sub(".*yearlyMean_tmax_([0-9.]+)\\.tif$", "\\1", all.tmax)
  )
  
  ordered.h <- order(heights_found)
  all.tmax <- all.tmax[ordered.h]
  heights <- heights_found[ordered.h]
  
  df_new <- list()
  
  for (f in seq_along(forest.types)) {
    
    ft <- forest.types[f]
    
    zones.file <- paste0(
      strata.path, "StrataZones_ft_", ft, ".RDS"
    )
    
    if (!file.exists(zones.file)) {
      warning("Skipping structural class: missing ", zones.file,
              call. = FALSE)
      next
    }
    
    zones <- readRDS(zones.file)
    
    for (z in seq_len(nrow(zones))) {
      
      z1 <- zones[z, ]
      
      files_max <- all.tmax[
        heights >= z1$zone_start &
          heights < z1$zone_end
      ]
      
      if (length(files_max) == 0L) next
      
      # Align each available height to the same 21 years
      height_stacks <- lapply(
        files_max,
        read_tmax,
        template = pca,
        years = yr.seq
      )
      
      height_stacks <- Filter(Negate(is.null), height_stacks)
      
      if (length(height_stacks) == 0L) next
      
      tmax.r <- do.call(c, height_stacks)
      
      # Mask to structural class
      tmax.r[pca != ft] <- NA
      
      # Safe because each height now has exactly 21 aligned layers
      year_index <- rep(
        yr.seq,
        times = length(height_stacks)
      )
      
      zone_mean <- tapp(
        tmax.r,
        index = year_index,
        fun = mean,
        na.rm = TRUE
      )
      
      zone_mean[!is.finite(zone_mean)] <- NA
      names(zone_mean) <- as.character(yr.seq)
      
      vals <- as.data.frame(values(zone_mean))
      names(vals) <- as.character(yr.seq)
      
      n.years <- length(yr.seq)
      
      # Load distinctness
      cd.file <- paste0(
        out.data, "climate_distinctness/", sample.name,
        "/ft_", ft, "_Mean_allyears_",
        z1$zone_start, "_", z1$zone_end, ".tif"
      )
      
      cd_vec <- rep(NA_real_, nrow(vals))
      
      if (file.exists(cd.file)) {
        
        cd <- tryCatch(
          rast(cd.file),
          error = function(e) NULL
        )
        
        if (!is.null(cd)) {
          compareGeom(pca, cd, stopOnError = TRUE)
          
          # Retains your original use of the first distinctness layer
          cd_vec <- values(mean(cd, na.rm = TRUE), mat = FALSE)
          cd_vec[!is.finite(cd_vec)] <- NA
        } else {
          warning("Unreadable distinctness file: ", cd.file,
                  call. = FALSE)
        }
        
      } else {
        warning("Missing distinctness file: ", cd.file,
                call. = FALSE)
      }
      
      df <- data.frame(
        cell = rep(seq_len(nrow(vals)), times = n.years),
        year = rep(yr.seq, each = nrow(vals)),
        mean = as.vector(as.matrix(vals)),
        zone = paste0(z1$zone_start, "-", z1$zone_end, "m"),
        type = ft,
        dist = rep(cd_vec, times = n.years)
      )
      
      df <- df[is.finite(df$mean), , drop = FALSE]
      
      df_new[[length(df_new) + 1L]] <- df
    }
  }
  
  df_grid <- dplyr::bind_rows(df_new)
  
  if (nrow(df_grid) == 0L) next
  
  df_grid$site <- sample.name
  site_list[[sample.name]] <- df_grid
}


df_allsites <- bind_rows(site_list) %>%
  mutate(
    # Cell numbers restart within each site, so create a globally unique ID
    cell_id = interaction(site, cell, drop = TRUE)
  )

# Calculate one warming slope per cell × forest type × zone
# All sites now belong to the same pooled analysis
slopes <- df_allsites %>%
  dplyr::filter(
    is.finite(mean),
    is.finite(dist),
    !is.na(year),
    !is.na(site),
    !is.na(cell),
    !is.na(cell_id),
    !is.na(zone),
    !is.na(type)
  ) %>%
  dplyr::group_by(site, type, zone, cell_id, cell) %>%
  dplyr::filter(dplyr::n_distinct(year) > 1) %>%
  dplyr::summarise(
    estimate = coef(lm(mean ~ year))[["year"]],
    dist = mean(dist, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(slopes, paste0(out.data,"climate_distinctness/slopes_data.RDS"))

# One pooled point for each forest-type × zone combination
slopes_sum <- slopes %>%
  dplyr::group_by(type, zone) %>%
  dplyr::summarise(
    slope_med = median(estimate, na.rm = TRUE),
    slope_lo  = quantile(estimate, 0.25, na.rm = TRUE),
    slope_hi  = quantile(estimate, 0.75, na.rm = TRUE),
    dist_med  = median(dist, na.rm = TRUE),
    n_cells   = n(),
    n_sites   = n_distinct(site),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    zone_start = as.numeric(sub("-.*", "", zone)),
    zone_end   = as.numeric(sub(".*-", "", sub("m", "", zone))),
    med_height = (zone_start + zone_end) / 2,
    slope_med_decade = slope_med * 10,
    slope_lo_decade  = slope_lo * 10,
    slope_hi_decade  = slope_hi * 10,
    type = factor(type)
  )

label_df <- slopes_sum %>%
  arrange(desc(dist_med)) %>%
  slice_head(n = 12)

p <- ggplot(
  slopes_sum,
  aes(
    x = dist_med,
    y = slope_med_decade,
    colour = factor(type)
  )
) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    colour = "grey30",
    fill = "grey80",
    linewidth = 0.4,
    alpha = 0.6,
    se = TRUE
  ) +
  scale_x_log10() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey70"
  ) +
  geom_point(
    aes(size = med_height),
    alpha = 0.8
  ) +
  geom_point(
    colour = "black",
    size = 0.3,
    show.legend = FALSE
  ) +
  scale_size_continuous(
    name = "Stratum midpoint (m)",
    range = c(5, 1)
  ) +
  scale_colour_manual(
    values = c(
      "1" = "#ca763b",
      "2" = "#852e47",
      "3" = "#333333",
      "4" = "#a09c33"
    )
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = gsub("-", " - ", zone)),
    colour = "grey30",
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.4,        # space around labels
    point.padding = 0.4,      # space around points
    force = 3,               # stronger repulsion
  ) +
  labs(
    x = "Vertical climate distinctness (log scale)",
    y = expression(
      paste("Warming trend (", degree, "C decade"^-1, ")")
    ),
    colour = "Structural class") +
  theme_bw()

p

ggsave(
  paste0(out.data, "figures/Fig4.png"),
  plot = p,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)
