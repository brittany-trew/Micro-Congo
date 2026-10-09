#' Fig 2X. Boxplots of outputs.
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

plot.out <- paste0(out.data,"plots/")

list.sd <- list()
for(b in 1:length(all.samples)){
  sample.name <- all.samples[[b]]
  inpath <- paste0(out.data, sample.name, "/verticalSummary/")
  
  files <- list.files(
    inpath,
    pattern = "^verticalSD_.*\\.tif$",
    full.names = TRUE
  )
  
  if (length(files) == 0L) next
  
  colsd_all <- terra::rast(files)
  colsd_mean <- mean(colsd_all, na.rm = T)
  vals <- as.vector(terra::values(colsd_mean, mat = TRUE))
  vals <- vals[is.finite(vals)]
  
  df1 <- data.frame(
    val = vals,
    sample = sample.name,
    variable = "Daily Max (SD)"
  )
  
  files <- list.files(
    inpath,
    pattern = "^colmax95_.*\\.tif$",
    full.names = TRUE
  )
  
  colmax_all <- terra::rast(files)
  colmax_mean <- mean(colmax_all, na.rm = T)
  vals <- as.vector(terra::values(colmax_mean, mat = TRUE))
  vals <- vals[is.finite(vals)]
  
  df2 <- data.frame(
    val = vals,
    sample = sample.name,
    variable = "Daily Max (95th Percentile)"
  )
  
  hot_daysl <- terra::rast(paste0(plot.out,
    "mean_vertical_range_hottest5pct_",
    sample.name,
    ".tif"
  ))
  vals <- as.vector(terra::values(hot_daysl, mat = TRUE))
  vals <- vals[is.finite(vals)]
  
  df3 <- data.frame(
    val = vals,
    sample = sample.name,
    variable = "Range of Daily Max (Hottest Days)"
  )
  
  df <- rbind(df1,df2,df3)
  
  list.sd[[b]] <- df
}

# Assign the first three colours to your categories
plot.colours <- setNames(
  unname(colour.p[1:3]),
  levels(df$variable)
)

# Sample points separately within each category
set.seed(123)

point.data <- df %>%
  group_by(variable) %>%
  slice_sample(n = 1000) %>%
  ungroup()

ggplot(df, aes(x = variable, y = val)) +
  geom_boxplot(
    aes(fill = variable),
    width = 0.5,
    alpha = 0.8,
    colour = "#303030",
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_point(
    data = point.data,
    aes(colour = variable),
    position = position_jitter(width = 0.16, height = 0, seed = 123),
    size = 0.4,
    alpha = 0.3
  ) +
  scale_fill_manual(values = plot.colours) +
  scale_colour_manual(values = plot.colours) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = NULL, y = "Temperature (°C)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "#303030"),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.line = element_line(colour = "#505050", linewidth = 0.4),
    axis.ticks = element_line(colour = "#505050"),
    plot.margin = margin(15, 20, 15, 15)
  )

ggsave(
  filename = paste0(plot.out,"summary_boxplots.png"),
  width = 8,
  height = 4,
  units = "in",
  dpi = 300,
  bg = "white"
)
