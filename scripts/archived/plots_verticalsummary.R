#' Kernel density plots of the monthly tmin/tmax environment
#' By Season: Dec-Feb, Mar-May, Jun-Aug, Sep-Nov.
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

#' Helper function:
#' Helper: sample cells early, before long dataframe is created
sample_cells <- function(vals, prop = 0.10) {
  ok_cells <- which(rowSums(is.finite(vals)) > 0)
  sample(ok_cells, size = floor(length(ok_cells) * prop))
}

# ── Setup ---------------------------------------------------------------------
dzd <- 1
set.seed(123)
final.list <- list()

for(b in 2:length(all.samples)){
  
  sample.name <- all.samples[[b]]
  print(sample.name)
  
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  chm       <- rast(paste0(head.path, "/chm.tif"))
  
  pai.layers    <- list.files(paste0(head.path, "/pai"))
  n.h           <- length(pai.layers)
  model.heights <- c(0.2, seq(dzd, (n.h - 1) * dzd, by = dzd))
  
  height_dirs <- list.dirs(
    paste0(inpath, "/dailyTemps/"),
    recursive = FALSE,
    full.names = TRUE
  )
  
  df_years <- list()
  
  for(i in 1:length(yr.seq)){
    
    yr <- yr.seq[i]
    print(yr)
    if(yr==2007){next()}
    dts <- seq(
      as.Date(paste0(yr, "-01-01")),
      as.Date(paste0(yr, "-12-31")),
      by = "day"
    )
    
    idx    <- format(dts, "%Y-%m")
    months <- idx[!duplicated(idx)]
    
    # ── Monthly max -----------------------------------------------------------
    f_max <- list.files(
      height_dirs,
      pattern = paste0("DailyMax_", yr, "\\.tif$"),
      full.names = TRUE
    )
    
    h_max <- as.numeric(sub("m$", "", basename(dirname(f_max))))
    f_max <- f_max[order(h_max)]
    
    r.list <- list()
    
    for(a in seq_along(f_max)){
      r <- rast(f_max[[a]])
      monthly_max <- tapp(r, index = idx, fun = mean, na.rm = TRUE)
      names(monthly_max) <- months
      r.list[[a]] <- monthly_max
    }
    
    max.stk <- do.call(c, r.list)
    rm(r.list, r, monthly_max)
    gc()
    
    # ── Monthly min -----------------------------------------------------------
    f_min <- list.files(
      height_dirs,
      pattern = paste0("DailyMin_", yr, "\\.tif$"),
      full.names = TRUE
    )
    
    h_min <- as.numeric(sub("m$", "", basename(dirname(f_min))))
    f_min <- f_min[order(h_min)]
    
    template <- rast(f_min[[1]])
    
    r.list <- list()
    
    for(a in seq_along(f_min)){
      r <- rast(f_min[[a]])
      r <- resample(r, template)
      
      monthly_min <- tapp(r, index = idx, fun = mean, na.rm = TRUE)
      names(monthly_min) <- months
      r.list[[a]] <- monthly_min
    }
    
    min.stk <- do.call(c, r.list)
    rm(r.list, r, monthly_min)
    gc()
    
    heights <- sort(as.numeric(sub("m$", "", basename(dirname(f_min)))))
    
    # ── Early sample ----------------------------------------------------------
    vals_max <- values(max.stk)
    vals_min <- values(min.stk)
    
    keep_cells <- sample_cells(vals_max, prop = 0.10)
    
    vals_max <- vals_max[keep_cells, , drop = FALSE]
    vals_min <- vals_min[keep_cells, , drop = FALSE]
    
    # ── Build long dataframe ------------------------------------------------------
    month_vec  <- rep(months, times = length(heights))
    height_vec <- rep(heights, each = length(months))
    
    df_long <- data.frame(
      cell   = rep(seq_len(nrow(vals_max)), times = ncol(vals_max)),
      tmax   = as.vector(vals_max),
      tmin   = as.vector(vals_min),
      month  = rep(month_vec, each = nrow(vals_max)),
      height = rep(height_vec, each = nrow(vals_max)),
      year   = yr
    )
    
    df_long <- df_long[
      complete.cases(df_long) &
        is.finite(df_long$tmax) &
        is.finite(df_long$tmin),
    ]
    
    # ── Sample early, balanced by height -----------------------------------------
    height_counts <- table(df_long$height)
    target_n <- floor(quantile(height_counts, 0.10))
    
    df_long <- df_long %>%
      group_by(height) %>%
      slice_sample(n = min(n(), target_n)) %>%
      ungroup()
    
    df_years[[as.character(yr)]] <- df_long
    
    rm(max.stk, min.stk, vals_max, vals_min, df_long)
    gc()
  }
  
  df <- do.call(rbind, df_years)
  
  df <- df %>%
    mutate(
      month_num = as.integer(substr(month, 6, 7)),
      season = seasons[as.character(month_num)]
    )
  
  all.3d <- list(
    SD = df |> filter(season == "Dec/Jan/Feb"),
    SR = df |> filter(season == "Mar/Apr/May"),
    LD = df |> filter(season == "Jun/Jul/Aug"),
    LR = df |> filter(season == "Sep/Oct/Nov")
  )
  
  final.list[[sample.name]] <- all.3d
  
  rm(df, df_years, all.3d)
  gc()
}

#' Combine all sites
df_all <- dplyr::bind_rows(final.list)
saveRDS(df_all, "allSites_seasons.RDS")

plotlist <- list()
for(a in 1:length(df_all)){
  
  season.df <- df_all[[a]]
  season <- unique(season.df$season)

  p <- ggplot() +
    geom_point(
      data = season.df,
      aes(tmin, tmax, colour = height),
      alpha = 0.1,
      size = 0.2
    ) +
    scale_colour_gradientn(
      colours = c("#769c43","#5bbda7","#237891","#864ead","#7c456c","#94375f","#db4655","#db9346")
    ) +
    theme_bw()+
    labs(
      x = "Monthly Minimum Temperature - °C",
      y = "Monthly Maximum Temperature - °C",
      title = season,
      colour = "Height"
    )
  
  plotlist[[a]] <- p
  }

require(patchwork)
final_grid <- wrap_plots(plotlist, ncol = 2, nrow = 2)+
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")  # change ncol to taste
print(final_grid)
print("done.")

ggsave(paste0(out.data,"plots/ThermalSpaceSeasonal.png"), plot = final_grid,
       width = 10, height = 10, units = "in", dpi = 300)

