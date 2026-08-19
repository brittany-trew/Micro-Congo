#' Vertical Decoupling. 
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

df.list <- list()
for(i in 4:length(all.samples)){
  sample.name <- all.samples[[i]]
  print(sample.name)
  
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  inpath  <- paste0(volumes3, sample.name, "/")
  outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
  dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
  
  chm           <- rast(paste0(head.path, "/chm.tif"))
  pai.layers    <- list.files(paste0(head.path, "/pai"))
  n.h           <- length(pai.layers)
  model.heights <- seq(0, (n.h - 1) * dzd, by = dzd)
  yr.seq        <- seq(2004, 2024, 1)
  
  height_dirs <- list.dirs(paste0(inpath, "/dailyTemps/"), recursive = FALSE, full.names = TRUE)
  heights     <- as.numeric(gsub("m$", "", basename(height_dirs)))
  height_dirs <- height_dirs[order(heights)]
  heights     <- sort(heights)
  
  df_years <- list()
  #' FOR EACH YEAR.
  for(b in 5:length(yr.seq)){
    print(yr.seq[[b]])
    yr=yr.seq[[b]]
    dts <- seq(as.Date(paste0(yr,"-01-01")), as.Date(paste0(yr,"-12-31")), by = "day")
    idx <- format(dts, "%Y-%m")
    
    #' ALL HEIGHTS.
    #' Monthly Max Temps
    f <- list.files(
      height_dirs,
      pattern = paste0("DailyMax_", yr, "\\.tif$"),
      full.names = TRUE
    )
    h <- as.numeric(sub("m$", "", basename(dirname(f))))
    ord <- order(h)
    f <- f[ord]
    h <- h[ord]
    
    #' FOR EACH HEIGHT...
    #' Derive Mean Monthly Maxima
    r.list <- list()
    for(a in 1:length(f)){
      
      r <- load_masked(f[[a]], chm)
      monthly_max <- tapp(r, index = idx, fun = mean, na.rm = TRUE)
      names(monthly_max) <- idx[!duplicated(idx)]
      r.list[[a]] <- monthly_max
    }
    
    month.list <- lapply(1:12, function(i) {
      s <- rast(lapply(r.list, function(x) x[[i]]))
      names(s) <- paste0("h_", h, "m")
      s
    })
    names(month.list) <- names(r.list[[1]])
    
    df_mons <- list()
    #' Now calc refuge for each month (each item in list)
    for(c in 1:length(month.list)){
      max.stk <- mask(month.list[[c]], chm > 0.2, maskvalues = FALSE) # Mask to forested pixels
      max.stk[max.stk > 55] <- NA
      
      mn <- names(month.list)[c]
      # Extract values
      vals_mat <- values(max.stk)

      # Remove all-NA pixels first
      row_min <- apply(vals_mat, 1, min, na.rm = TRUE)
      valid_pixels <- is.finite(row_min)
      vals_mat <- vals_mat[valid_pixels, , drop = FALSE]
      
      # Global min from valid data
      global_min <- min(vals_mat, na.rm = TRUE)
      # Site minimum per height band
      height_min <- apply(vals_mat, 2, min, na.rm = TRUE)
     
      df_refuge <- do.call(rbind, lapply(seq_along(h), function(j) {
        current_temp   <- vals_mat[, j]
        horiz_refuge   <- current_temp - height_min[j]
        full_3d_refuge <- current_temp - global_min
        
        data.frame(
          date           = mn,
          height         = heights[j],
          full_3d_refuge = full_3d_refuge,
          horiz_refuge   = horiz_refuge
        )
      }))
      
      df_refuge <- df_refuge |>
        mutate(
          relative_gain = ifelse(
            horiz_refuge == 0, NA,
            100 * (full_3d_refuge - horiz_refuge) / horiz_refuge
          ),
          absolute_gain = full_3d_refuge - horiz_refuge
        )
      
      df_refuge <- na.omit(df_refuge)
      df_mons[[mn]] <- df_refuge
    } # END MONTH LOOP
    df_years[[as.character(yr)]] <- bind_rows(df_mons)    
  
    }# END YEAR LOOP

  df.list[[sample.name]] <- bind_rows(df_years, .id = "year")
  
} #END SITE LOOP

df_plot <- bind_rows(df.list, .id = "site") %>%
  distinct(site, year, date, height, absolute_gain, .keep_all = TRUE)

sum_df <- df_plot %>%
  group_by(site, height) %>%
  summarise(
    med = median(absolute_gain, na.rm = TRUE),
    q25 = quantile(absolute_gain, 0.25, na.rm = TRUE),
    q75 = quantile(absolute_gain, 0.75, na.rm = TRUE),
    q05 = quantile(absolute_gain, 0.05, na.rm = TRUE),
    q95 = quantile(absolute_gain, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

plot_list <- lapply(unique(sum_df$site), function(s) {
  
  ggplot(filter(sum_df, site == s), aes(height, med)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
    geom_smooth(linewidth = 0.8, colour = "#db4655", se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#237891") +
    labs(
      title = s,
      x = "Height above ground (m)",
      y = "Cooling gain (°C)"
    ) +
    theme_bw()
})

require(patchwork)
p <- wrap_plots(plot_list, ncol = 2)

ggsave(
  paste0(out.data, "plots/buffer_all_sites.png"),
  plot = p,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)


sum_df <- df.list$Imbalanga %>%
  group_by(height) %>%
  summarise(
    med = median(absolute_gain, na.rm = TRUE),
    q25 = quantile(absolute_gain, 0.25, na.rm = TRUE),
    q75 = quantile(absolute_gain, 0.75, na.rm = TRUE),
    q05 = quantile(absolute_gain, 0.05, na.rm = TRUE),
    q95 = quantile(absolute_gain, 0.95, na.rm = TRUE),
    .groups = "drop"
  )
p <-ggplot(sum_df, aes(height, med)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
  geom_smooth(linewidth = 0.8, colour = "#db4655", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#237891") +
  labs(
    x = "Height above ground (m)",
    y = "Cooling gain (°C)"
  ) +
  theme_bw()
ggsave(
  paste0(out.data, "plots/buffer_Imbalanga.png"),
  plot = p,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)
