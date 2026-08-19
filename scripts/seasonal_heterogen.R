#' Plot the two SD values, looking at variation in seasons, not over the entire twenty years. 

# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

terraOptions(
  tempdir = "/Volumes/MicroMaze2/tmp",
  memfrac = 0.7
)

yr.seq        <- seq(2004, 2024, 1)

sample_cells <- function(vals, prop = 0.10) {
  ok_cells <- which(rowSums(is.finite(vals)) > 0)
  sample(ok_cells, size = floor(length(ok_cells) * prop))
}

SD_vertical <- function(files, chm) {
  r.lst  <- lapply(files, load_masked, chm = chm)
  r.stk <- do.call(c, r.lst)
  r.stk[r.stk < 0] <- NA
  
  n_days <- nlyr(r.lst[[1]])
  index  <- rep(1:n_days, times = length(r.lst))
  colsd <- tapp(r.stk, index, fun = sd, na.rm = TRUE)
  gc()
  return(colsd)
}

final.list <- list()
for(b in 1:length(all.samples)){
  
  #' Site set up.
  sample.name <- all.samples[[b]]
  print(sample.name)
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  inpath  <- paste0(volumes3, sample.name, "/")
  outpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
  dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
  
  chm       <- rast(paste0(head.path, "/chm.tif"))
  
  pai.layers    <- list.files(paste0(head.path, "/pai"))
  n.h           <- length(pai.layers)
  model.heights <- c(0.2, seq(dzd, (n.h - 1) * dzd, by = dzd))
  
  height_dirs <- list.dirs(
    paste0(inpath, "/dailyTemps/"),
    recursive = FALSE,
    full.names = TRUE
  )
  
  # Monthly SD of Vertical Column
  df_years <- list()
  for(i in 1:length(yr.seq)){
    
    yr <- yr.seq[i]
    print(yr)
    
    #' (1) VERTICAL.
    r <- rast(paste0(outpath, "verticalSD_",  yr, ".tif"))

    #' Convert Days to Months
    dts <- seq(
      as.Date(paste0(yr, "-01-01")),
      as.Date(paste0(yr, "-12-31")),
      by = "day"
    )
    idx    <- format(dts, "%Y-%m")
    months <- unique(idx)
    
    monthly_vsd <- tapp(
      r,
      index = idx,
      fun = mean,
      na.rm = TRUE
    )
    
    names(monthly_vsd) <- months
    
    #' (1) TIME.
    #' Daily rasters:
    f_max <- list.files(
      height_dirs,
      pattern = paste0("DailyMax_", yr, "\\.tif$"),
      full.names = TRUE
    )
    #' Ordered by height:
    h_max <- as.numeric(sub("m$", "", basename(dirname(f_max))))
    f_max <- f_max[order(h_max)]

    # Time SD: for each height, monthly SD through daily Tmax
    time_sd_list <- list()
    for(a in seq_along(f_max)) {
      r <- load_masked(f_max[[a]], chm)
      r[r < 0] <- NA
      # SD through time within each month at this height
      time_sd_list[[length(time_sd_list) + 1]] <- tapp(
        r,
        index = idx,
        fun = sd,
        na.rm = TRUE
      )
    }
    
    # stack all height-specific monthly time SD rasters
    time_sd_stk <- do.call(c, time_sd_list)
    # Mean across heights for each month
    time_sd_monthly <- list()
    for(j in seq_along(months)) {
      month_layers <- seq(
        from = j,
        to = nlyr(time_sd_stk),
        by = length(months)
      )
      time_sd_monthly[[j]] <- app(time_sd_stk[[month_layers]], mean, na.rm = TRUE)
    }
    
    time_sd <- do.call(c, time_sd_monthly)
    names(time_sd) <- months
    gc()
    
    #' Build data frame.
    vals_sd   <- values(monthly_vsd, na.rm = T)
    vals_time <- values(time_sd, na.rm = T)
    
    keep_cells <- sample_cells(vals_sd, prop = 0.25)
    
    vals_vertical <- vals_sd[keep_cells, , drop = FALSE]
    vals_time     <- vals_time[keep_cells, , drop = FALSE]
    
    df_long <- data.frame(
      cell       = rep(keep_cells, times = ncol(vals_vertical)),
      VerticalSD = as.vector(vals_vertical),
      TimeSD     = as.vector(vals_time),
      month      = rep(months, each = length(keep_cells)),
      year       = yr,
      site = rep(sample.name, times = ncol(vals_vertical))
    )

    df_years[[as.character(yr)]] <- df_long
    rm(monthly_vsd, time_sd, time_sd_stk, time_sd_list, vals_sd, vals_time)
    gc()
  }
  
  df <- do.call(rbind, df_years)
  
  df <- df %>%
    mutate(
      month_num = as.integer(substr(month, 6, 7)),
      season = seasons[as.character(month_num)]
    )
  
  final.list[[sample.name]] <- df
  
  rm(df, df_years)
  gc()
}

df_all <- dplyr::bind_rows(final.list)

season_order <- c("Dec/Jan/Feb", "Mar/Apr/May", "Jun/Jul/Aug", "Sep/Oct/Nov")

require(patchwork)
plotlist <- lapply(season_order, function(s) {
  
  season.df <- df_all |> 
    filter(season == s)
  
  ggplot(season.df) +
    geom_point(
      aes(x = TimeSD, y = VerticalSD, colour = site),
      alpha = 0.1,
      size = 0.2
    ) +
    scale_colour_manual(values = c("#769c43", "#5bbda7", "#db4655", "#db9346")) +
    theme_bw() +
    labs(
      x = "Temporal SD - °C",
      y = "Vertical SD - °C",
      title = s,
      colour = "Site"
    )
})

final_grid <- wrap_plots(plotlist, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(final_grid)
print("done.")

ggsave(
  paste0(out.data, "plots/ThermalSpaceSeasonal.png"),
  plot = final_grid,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)


p1 <- ggplot(df_all,
             aes(x = season, y = VerticalSD, fill = season)) +
  geom_violin(trim = TRUE, alpha = 0.8, colour = NA) +
  geom_boxplot(
    aes(colour = season),
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    colour = "grey40",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = c("#769c43", "#5bbda7", "#db4655", "#db9346")) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Mean vertical SD of daily Tmax - °C"
  ) +
  theme(legend.position = "none")


p2 <- ggplot(df_all,
             aes(x = season, y = TimeSD, fill = season)) +
  geom_violin(trim = TRUE, alpha = 0.8, colour = NA) +
  geom_boxplot(
    aes(colour = season),
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    colour = "grey40",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = c("#769c43", "#5bbda7", "#db4655", "#db9346")) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Temporal SD of vertical mean of daily Tmax - °C"
  ) +
  theme(legend.position = "none")

final.p <- p1 / p2
final.p

ggsave(
  paste0(out.data, "plots/HeteroSeasonalBoxplot.png"),
  plot = final.p,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)


df_all %>%
  group_by(season) %>%
  summarise(
    median_vertical = median(VerticalSD, na.rm = TRUE),
    median_time = median(TimeSD, na.rm = TRUE)
  )
