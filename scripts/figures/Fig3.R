#' Extract forest types by grid cells using PCAs
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

in.path <- paste0(out.data, "pca_clusters/")

profiles <- read.csv(paste0(outpath,"mean_profiles.csv"))

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
  
  chm <- rast(
    paste0(head.path, "chm.tif")
  )
  
  max_heights[i] <- global(
    chm,
    "max",
    na.rm = TRUE
  )[1, 1]
}
study_max_height <- max(max_heights)
study_max_height

profiles <- filter(profiles, height_m <= study_max_height)

#' Mean and summary for all the sites & classes:
mean_profiles <- profiles %>%
  dplyr::group_by(Forest.Type, height_m) %>%
  dplyr::summarise(
    pad = mean(pad, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(Forest.Type, height_m)

overall_summary <- profiles %>%
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
    legend.position = "right")+
  ylim(0,40)

#' Plot
final_plot

ggsave(
  paste0(out.data, "plots/Vertical_Profiles_All_Sites.png"),
  plot = final_plot,
  width = 5,
  height = 4.5,
  units = "in",
  dpi = 300)


#" Figure 4B (Schematic only).
filt.data <- filter(profiles,
                    Forest.Type == 4)
filt.profile <- filt.data %>%
  dplyr::group_by(height_m) %>%
  dplyr::summarise(
    pad = mean(pad, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(height_m)

#' Add inflection points.
outpath <- paste0(out.data,"strata/")
strata <- readRDS(paste0(outpath,"StrataZones_ft_4.RDS"))
strata <- strata %>%
  dplyr::mutate(
    zone_start = ifelse(
      zone_start == 0,
      0,
      pmax(0.5, floor(zone_start / dzd + 0.5) * dzd)
    ),
    zone_end = pmax(
      0.5,
      floor(zone_end / dzd + 0.5) * dzd
    )
  )
boundary_points <- strata %>%
  summarise(
    height_m = list(unique(c(zone_start, zone_end)))) %>%
  tidyr::unnest(height_m) %>%
  arrange(height_m) %>%
  left_join(
    filt.profile %>%
      dplyr::select(height_m, pad),
    by = "height_m") %>%
  dplyr::slice(-c(1, n()))
  

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
  theme_classic()+
  ylim(0,40)
filt_plot
ggsave(
  paste0(out.data, "plots/Vertical_Profiles_class4example.png"),
  plot = filt_plot,
  width = 4.5,
  height = 6,
  units = "in",
  dpi = 300
)

