#' Extract forest types by grid cells using PCAs
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

in.path <- paste0(out.data, "pca_clusters/")

#' Function to extract PAD data for a LiDAR site.
extract_pad_profiles <- function(sample.name,
                                 site.location,
                                 year,
                                 out.data.old,
                                 in.path,
                                 dzd) {
  
  head.path <- paste0(
    out.data.old, site.location, "/", sample.name, "/", year, "/"
  )
  
  pca <- rast(paste0(in.path, sample.name, "_pca_Clust4.tif"))
  chm <- rast(paste0(head.path, "chm.tif"))
  
  pad.processed <- processPAD(head.path, dzd, chm)
  pad.aboveG   <- pad.processed[[1]]
  start_heights <- pad.processed[[2]]
  
  pca <- resample(
    pca,
    pad.aboveG[[1]],
    method = "near"
  )
  
  pca_values <- values(pca, mat = FALSE)
  
  # PAD values assigned to PCA structural class
  type_data <- bind_rows(
    lapply(seq_len(nlyr(pad.aboveG)), function(i) {
      
      data.frame(
        Site = sample.name,
        Forest.Type = pca_values,
        height_m = start_heights[i],
        pad = values(pad.aboveG[[i]], mat = FALSE)
      )
    })
  ) %>% filter(
      Forest.Type %in% 1:4,
      height_m < 30,
      !is.na(pad)
    )
  
  # All PAD values, irrespective of forest type
  overall_data <- bind_rows(
    lapply(seq_len(nlyr(pad.aboveG)), 
           function(i) {
      data.frame(
        Site = sample.name,
        height_m = start_heights[i],
        pad = values(pad.aboveG[[i]], mat = FALSE))
    })) %>% 
    filter(
      height_m < 30,
      !is.na(pad))
  
  list(
    type_data = type_data,
    overall_data = overall_data
  )
}


#' Apply function to all sites.
site_profiles <- lapply(all.samples, function(s) {
  message("Extracting: ", s)
  extract_pad_profiles(
    sample.name = s,
    site.location = site.location,
    year = year,
    out.data.old = out.data.old,
    in.path = in.path,
    dzd = dzd
  )
})

type_data_all <- bind_rows(
  lapply(site_profiles, `[[`, "type_data")
)

overall_data_all <- bind_rows(
  lapply(site_profiles, `[[`, "overall_data")
)


#' Mean and summary for all the sites & classes:
mean_profiles <- type_data_all %>%
  dplyr::group_by(Forest.Type, height_m) %>%
  dplyr::summarise(
    pad = mean(pad, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(Forest.Type, height_m)

overall_summary <- overall_data_all %>%
  dplyr::group_by(height_m) %>%
  dplyr::summarise(
    mean_pad = mean(pad, na.rm = TRUE),
    lo = quantile(pad, 0.1, na.rm = TRUE),
    hi = quantile(pad, 0.9, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(height_m)


#" Figure 4A.
final_plot <- ggplot() +
  # 10-90th percentile ribbons for all sites.
  geom_ribbon(
    data = overall_summary,
    aes(
      y = height_m,
      xmin = lo,
      xmax = hi),
    fill = "grey80",
    alpha = 0.4) +
  # Overall mean profile for all sites
  geom_path(
    data = overall_summary,
    aes(
      x = mean_pad,
      y = height_m),
    linewidth = 0.5,
    colour = "grey40",
    alpha = 0.8,
    linetype = "dashed") +
  # Mean profile for each structural class
  geom_path(
    data = mean_profiles,
    aes(
      x = pad,
      y = height_m,
      colour = factor(Forest.Type),
      group = Forest.Type),
    linewidth = 0.7) +
  scale_colour_manual(
    values = c(
      "1" = "#ca763b",
      "2" = "#852e47",
      "3" = "#333333",
      "4" = "#a09c33"),
    labels = c(
      "1" = "1",
      "2" = "2",
      "3" = "3",
      "4" = "4")) +
  labs(
    x = "Plant Area Density",
    y = "Height (m)",
    colour = "Structural Class") +
  theme_classic() +
  theme(
    legend.position = "bottom")

#' Plot
final_plot

ggsave(
  paste0(out.data, "plots/Vertical_Profiles_All_Sites.png"),
  plot = final_plot,
  width = 4.5,
  height = 6,
  units = "in",
  dpi = 300)


#" Figure 4B (Schematic only).
filt.data <- filter(type_data_all,
                    Forest.Type == 1)
filt.profile <- filt.data %>%
  dplyr::group_by(height_m) %>%
  dplyr::summarise(
    pad = mean(pad, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(height_m)

#' Add inflection points.
outpath <- paste0(out.data,"MondoBai/veg_zoning/")
strata <- readRDS(paste0(outpath,"VegZoning_FT_1.RDS"))
boundary_points <- strata %>%
  summarise(
    height_m = list(unique(c(zone_start, zone_end)))) %>%
  tidyr::unnest(height_m) %>%
  arrange(height_m) %>%
  left_join(
    filt.profile %>%
      dplyr::select(height_m, pad),
    by = "height_m")

#' Plot mean line with inflection points and derived-strata.
filt_plot <- ggplot() +
  geom_path(
    data = filt.profile,
    aes(
      x = pad,
      y = height_m
    ),
    colour = "#333333",
    linewidth = 0.7
  ) +
  geom_point(
    data = boundary_points,
    aes(
      x = pad,
      y = height_m
    ),
    size = 2.5
  ) +
  geom_text(
    data = boundary_points,
    aes(
      x = pad,
      y = height_m,
      label = paste0(height_m, " m")
    ),
    nudge_x = 0.03,
    nudge_y = 0.1,
    hjust = 0,
    size = 3
  ) +
  labs(
    x = "Plant Area Density",
    y = "Height (m)"
  ) +
  theme_classic()
filt_plot
ggsave(
  paste0(out.data, "plots/Vertical_Profiles_class1example.png"),
  plot = filt_plot,
  width = 4.5,
  height = 6,
  units = "in",
  dpi = 300
)


count.data <- filter(filt.data,
                    height_m == 0.5)
nrow(count.data)
