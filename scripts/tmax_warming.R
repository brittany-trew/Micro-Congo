#' Warming of tmax - timeseries
# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

all.samples <- all.samples[-1]
dzd = 1
yr.seq        <- seq(2004, 2024, 1)

site_list <- list()
for(a in 1:length(all.samples)){
  sample.name <- all.samples[[a]]
  if(sample.name=="MobaCluster"){
    yr.seq <- yr.seq[-4]
    }else{yr.seq <- seq(2004, 2024, 1)}
  #' Site set up.
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  inpath  <- paste0(volumes3, sample.name, "/")
  outpath <- paste0(inpath, "monthlyTemps/")
  dir.create(outpath, showWarnings = F)
  
  in.path <- paste0(out.data,sample.name,"/veg_zoning/")
  
  chm           <- rast(paste0(head.path, "/chm.tif"))
  pca.path <- paste0(out.data,"pca_clusters/")
  pca <- rast(paste0(pca.path, sample.name, "_pca_Clust4.tif"))

  forest.types <- unique(values(pca, na.rm = T))
  all.tmax <- list.files(outpath, pattern = "yearlyMean_tmax.*.tif",
                         full.names = T)
  heights_found <- as.numeric(
    sub(".*yearlyMean_tmax_([0-9.]+)\\.tif$", "\\1", all.tmax)
  )
  
  #' --------------------------------------------------------------------------
  #' --------------------------------------------------------------------------
  #' Main loop: Forest type & Zone Combinations.
  df_new <- list()
  for(f in 1:length(forest.types)){
    #' For forest type f:
    ft <- forest.types[f]
    zones <- readRDS(paste0(in.path,"VegZoning_FT_",ft,".RDS"))
    for(z in 1:nrow(zones)){
      z1 <- zones[z,]
      files_max <- all.tmax[heights_found >= z1$zone_start & heights_found < z1$zone_end]
      h.names <- as.numeric(
        sub(".*yearlyMean_tmax_([0-9.]+)\\.tif$", "\\1", files_max)
      )
      if (length(files_max) == 0) {
        next}    
      
      tmax.r <- rast(files_max)
      tmax.r[is.infinite(tmax.r)] <- NA
      
      #' Mask to forest type.
      tmax.r[pca != ft] <- NA
      
      year_index <- rep(yr.seq, times = length(h.names))
      
      zone_mean <- tapp(
        tmax.r,
        index = year_index,
        fun = mean,
        na.rm = TRUE
      )
      
      names(zone_mean) <- yr.seq
      
      vals <- as.data.frame(values(zone_mean))
      names(vals) <- yr.seq
      
      n.years <- length(yr.seq)
      
      # Load in distinctness
      cd <- rast(paste0(out.data, "climate_distinctness/", sample.name, "/FT_",
                        ft, "_Mean_allyears_",
                        z1$zone_start, "_", z1$zone_end, ".tif"))
      
      cd_vec <- values(cd)[, 1]
      
      df <- data.frame(
        cell = rep(1:nrow(vals), times = n.years),
        year = rep(yr.seq, each = nrow(vals)),
        mean = as.vector(as.matrix(vals)),
        zone = paste0(z1$zone_start, "-", z1$zone_end, "m"),
        type = ft,
        dist = rep(cd_vec, times = n.years)
      )
      df <- df[is.finite(df$mean), ]
      
      df_new <- append(df_new, list(df))
    }
    
  }
  #' Main loop END
  #' --------------------------------------------------------------------------
  #' --------------------------------------------------------------------------

  #' Bind all type/zone combos:
  df_grid <- bind_rows(df_new)
  #' Add sample name label.
  df_grid$site <- sample.name
  
  site_list <- append(site_list, list(df_grid))
}

df_allsites <- bind_rows(site_list) %>%
  mutate(
    # Cell numbers restart within each site, so create a globally unique ID
    cell_id = interaction(site, cell, drop = TRUE)
  )

# Calculate one warming slope per cell × forest type × zone
# All sites now belong to the same pooled analysis
slopes <- df_allsites %>%
  filter(
    is.finite(mean),
    is.finite(dist),
    !is.na(year),
    !is.na(cell_id),
    !is.na(zone),
    !is.na(type)
  ) %>%
  group_by(type, zone, cell_id) %>%
  filter(n_distinct(year) > 1) %>%
  dplyr::summarise(
    estimate = coef(lm(mean ~ year))[["year"]],
    dist     = mean(dist, na.rm = TRUE),
    site     = first(site),  # retained for checking/model structure
    .groups  = "drop"
  )

# One pooled point for each forest-type × zone combination
slopes_sum <- slopes %>%
  group_by(type, zone) %>%
  dplyr::summarise(
    slope_med = median(estimate, na.rm = TRUE),
    slope_lo  = quantile(estimate, 0.25, na.rm = TRUE),
    slope_hi  = quantile(estimate, 0.75, na.rm = TRUE),
    dist_med  = median(dist, na.rm = TRUE),
    n_cells   = n(),
    n_sites   = n_distinct(site),
    .groups   = "drop"
  ) %>%
  mutate(
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
    x = slope_med_decade,
    y = dist_med,
    colour = type
  )
) +
  geom_errorbarh(
    aes(
      xmin = slope_lo_decade,
      xmax = slope_hi_decade
    ),
    alpha = 0.4,
    height = 0.08
  ) +
  geom_point(
    aes(size = med_height),
    alpha = 0.85
  ) +
  scale_size_continuous(
    name = "Zone midpoint (m)",
    range = c(5, 1)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey70"
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
    aes(label = zone),
    colour = "grey30",
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  labs(
    x = expression(
      paste("Warming trend (", degree, "C decade"^-1, ")")
    ),
    y = "Vertical climate distinctness",
    colour = "Structural class"
  ) +
  theme_bw()

p

ggsave(
  paste0(out.data, "plots/CDWarming_byZone.png"),
  plot = p,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)


#' Model.
slopes <- slopes %>%
  mutate(
    type = factor(type),
    zone = factor(zone),
    type_zone = interaction(type, zone, drop = TRUE),
    dist_c = dist - mean(dist, na.rm = TRUE)
  )

#' (1) are cells with higher climate distinctness warming faster or slower
m1 <- lm(estimate ~ dist_c, data = slopes)
summary(m1)

m2 <- lm(estimate ~ dist_c * type, data = slopes)
summary(m2)
anova(m1, m2)

m3 <- lm(estimate ~ dist_c * type + zone, data = slopes)
summary(m3)
anova(m2, m3)

slopes <- slopes %>%
  mutate(
    type_zone = interaction(type, zone, drop = TRUE)
  )
m_zone <- lm(
  estimate ~ type_zone,
  data = slopes
)
m_final <- lm(
  estimate ~ type_zone + dist_c:type,
  data = slopes
)
summary(m_final)
anova(m_zone, m_final)


coef_tab <- as.data.frame(coef(summary(m_final))) %>%
  tibble::rownames_to_column("term")

ci_tab <- as.data.frame(confint(m_final)) %>%
  tibble::rownames_to_column("term") %>%
  rename(
    conf_low = `2.5 %`,
    conf_high = `97.5 %`
  )

dist_effects <- coef_tab %>%
  left_join(ci_tab, by = "term") %>%
  filter(grepl("^dist_c:type", term)) %>%
  mutate(
    type = gsub("dist_c:type", "", term),
    effect_C_per_year = Estimate,
    effect_C_per_decade = Estimate * 10,
    conf_low_decade = conf_low * 10,
    conf_high_decade = conf_high * 10,
    p_value = `Pr(>|t|)`
  ) %>%
  select(
    type,
    effect_C_per_year,
    effect_C_per_decade,
    conf_low_decade,
    conf_high_decade,
    p_value
  )

dist_effects_final <- dist_effects %>%
  mutate(
    effect_C_per_year = round(effect_C_per_year, 4),
    effect_C_per_decade = round(effect_C_per_decade, 3),
    conf_low_decade = round(conf_low_decade, 3),
    conf_high_decade = round(conf_high_decade, 3),
    p_value = ifelse(p_value < 0.001, "<0.001", as.character(round(p_value, 3)))
  ) %>%
  rename(
    `Forest type` = type,
    `Effect, °C yr-1` = effect_C_per_year,
    `Effect, °C decade-1` = effect_C_per_decade,
    `Lower 95% CI, °C decade-1` = conf_low_decade,
    `Upper 95% CI, °C decade-1` = conf_high_decade,
    `P value` = p_value
  )

write_csv(dist_effects_final, paste0(out.data, "tables/dist_effects.csv"))

dp <- ggplot(
  dist_effects,
  aes(y = factor(type),
      x = effect_C_per_decade,
      xmin = conf_low_decade,
      xmax = conf_high_decade,
      colour = factor(type))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_pointrange(size = 0.8) +
  scale_colour_manual(
    values = c("#769c43", "#5bbda7", "#db4655", "#db9346"),
    labels = c("1", "2", "3", "4")
  ) +
  theme_bw() +
  labs(
    x = expression(paste("Effect of vertical climate distinctness on Tmax warming (", degree, "C decade"^-1, ")")),
    y = "Forest type",
    colour = "Forest type"
  )

ggsave(
  paste0(out.data, "plots/CDWarming_effectSize.png"),
  plot = dp,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)


models <- list(
  m1 = m1,
  m2 = m2,
  m3 = m3,
  m_zone = m_zone,
  m_final = m_final
)

model_comp <- tibble::tibble(
  model = names(models),
  description = c(
    "Climate distinctness only",
    "Climate distinctness × forest type",
    "Climate distinctness × forest type + vertical zone",
    "Forest type × vertical zone only",
    "Forest type × vertical zone + climate distinctness by forest type"
  ),
  r2 = sapply(models, function(x) summary(x)$r.squared),
  adj_r2 = sapply(models, function(x) summary(x)$adj.r.squared),
  rss = sapply(models, deviance),
  df_resid = sapply(models, df.residual),
  AIC = sapply(models, AIC)
) %>%
  mutate(
    delta_r2 = r2 - lag(r2)
  )

model_comp <- model_comp %>%
  mutate(
    r2 = round(r2, 3),
    adj_r2 = round(adj_r2, 3),
    delta_r2 = round(delta_r2, 3),
    rss = round(rss, 2),
    AIC = round(AIC, 1)
  )

write_csv(model_comp, paste0(out.data, "tables/model_comparison.csv"))



iqr_by_type <- slopes %>%
  group_by(type) %>%
  summarise(
    dist_iqr = IQR(dist, na.rm = TRUE),
    dist_range = diff(range(dist, na.rm = TRUE)),
    .groups = "drop"
  )

dist_effects_scaled <- dist_effects %>%
  left_join(iqr_by_type, by = "type") %>%
  mutate(
    effect_IQR_decade = effect_C_per_decade * dist_iqr,
    effect_IQR_20yrs  = effect_IQR_decade * 2,
    effect_range_decade = effect_C_per_decade * dist_range,
    effect_range_20yrs  = effect_range_decade * 2
  )

dist_effects_scaled
write_csv(dist_effects_scaled, paste0(out.data, "tables/dist_effects_scaled.csv"))

