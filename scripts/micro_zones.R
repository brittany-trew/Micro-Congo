#' Microclimate by Zones.
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R")) # loads worker functions

# ── Helper functions ----------------------------------------------------------
# Extract height (metres) from parent folder name
# e.g. ".../dailyTemps/5.0m/DailyMax_2009.tif" -> 5
get_height <- function(f) {
  as.numeric(gsub("m$", "", basename(dirname(f))))
}

# Load a single height-band raster, resample to CHM grid, then mask pixels
# where canopy height is below this height band (i.e. above canopy)
load_masked <- function(f, chm) {
  hh <- get_height(f)
  r  <- resample(rast(f), chm)
  ifel(chm < hh, NA, r)
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

#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Background Prep.
for(i in 1:length(all.samples)){
  
  sample.name <- all.samples[[i]]
  
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  chm <- rast(paste0(head.path, "chm.tif"))
  pca.path <- paste0(out.data,"pca_clusters/")
  pca <- rast(paste0(pca.path, sample.name, "_pca_Clust4.tif"))
  
  in.path <- paste0(out.data,sample.name,"/veg_zoning/")
  outpath <- paste0(out.data,sample.name,"/micro_zoned/")
  dir.create(outpath, showWarnings = FALSE)
  outpath <- paste0(outpath,"/all_years/")
  dir.create(outpath, showWarnings = FALSE)
  
  pai.layers    <- list.files(paste0(head.path, "/pai"))
  n.h           <- length(pai.layers)
  model.heights <- seq(0, (n.h - 1) * dzd, by = dzd)
  yr.seq        <- seq(2004, 2024, 1)
  
  inpath  <- paste0(volumes3, sample.name, "/")
  
  forest.types <- unique(values(pca, na.rm = T))
  
  for(f in 1:length(forest.types)){
    #' For forest type f:
    ft <- forest.types[f]
    zones <- readRDS(paste0(in.path,"VegZoning_FT_",ft,".RDS"))
    
    outpath <- paste0(out.data,sample.name,"/micro_zoned/")
    dir.create(outpath, , showWarnings = FALSE)
    outpath <- paste0(outpath,"/all_years/")
    dir.create(outpath, showWarnings = FALSE)
    outpath <- paste0(outpath,"/ForestType_",ft,"/")
    dir.create(outpath, showWarnings = FALSE)
    
    for(z in 1:nrow(zones)){
      
      z1 <- zones[z,]
      
      for (a in 1:length(yr.seq)) {
        
        yr <- yr.seq[[a]]
        if(sample.name=="MobaCluster"& yr == 2007){next()}
        message("Processing ", yr, " ...")
        if(!file.exists(paste0(
          outpath, "VerticalSDinDailyTmax_2004to2024_",z1$zone_start,"to",z1$zone_end,".tif"))
          ){
          
          files_max <- list.files(
            inpath,
            pattern    = paste0("DailyMax_", yr, "\\.tif$"),
            recursive  = TRUE,
            full.names = TRUE
          )
          
          heights_found <- as.numeric(get_height(files_max))
          files_max <- files_max[heights_found >= z1$zone_start & heights_found < z1$zone_end]
          
          # Hottest point in vertical column on each day
          # Coolest peak temperature in vertical column on each day
          # min of DailyMax across height bands = coolest a layer gets at its own peak
          sdVertical <- SD_vertical(files_max, chm)
          writeRaster(sdVertical, paste0(
            outpath, "sdVertical_", yr, "_",ft, "_",z1$zone_start,"to",z1$zone_end,".tif"),
            overwrite = TRUE)
          gc()
        }
        
      }
    }
  }
  
  #' Take the mean overall.
  for(f in 1:length(forest.types)){
    ft <- forest.types[f]
    zones <- readRDS(paste0(in.path,"VegZoning_FT_",ft,".RDS"))
    
    outpath <- paste0(out.data,sample.name,"/micro_zoned/")
    outpath <- paste0(outpath,"/all_years/")
    outpath <- paste0(outpath,"/ForestType_",ft,"/")
    
    for(z in 1:nrow(zones)){
      
      z1 <- zones[z,]
      
      if(!file.exists(paste0(outpath, "VerticalSDinDailyTmax_2004to2024_",z1$zone_start,"to",z1$zone_end,".tif"))
      ){
        sd_all  <- rast(list.files(outpath, pattern = paste0("^sdVertical.*",z1$zone_start,"to",z1$zone_end,".tif$"),  full.names = TRUE))
        mean_sd <- app(sd_all, mean, na.rm = TRUE)
        mean_sd <- ifel(chm < 0.2, NA, mean_sd)
        mean_sd[pca != ft] <- NA
        plot(mean_sd)
        writeRaster(
          mean_sd,
          paste0(outpath, "VerticalSDinDailyTmax_2004to2024_",z1$zone_start,"to",z1$zone_end,".tif"),
          overwrite = TRUE
        )
        file.remove(list.files(outpath, pattern = paste0("^sdVertical.*",z1$zone_start,"to",z1$zone_end,".tif$"),  full.names = TRUE))
      }
      
    }
  }
}


#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Collecting all Site Metrics:
df.list <- list()
for(i in 1:length(all.samples)){
  
  sample.name <- all.samples[[i]]
  
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  chm <- rast(paste0(head.path, "chm.tif"))
  pca.path <- paste0(out.data,"pca_clusters/")
  pca <- rast(paste0(pca.path, sample.name, "_pca_Clust4.tif"))
  
  in.path <- paste0(out.data,sample.name,"/veg_zoning/")
  outpath <- paste0(out.data,sample.name,"/micro_zoned/")
  dir.create(outpath, showWarnings = FALSE)
  outpath <- paste0(outpath,"/all_years/")
  dir.create(outpath, showWarnings = FALSE)
  
  pai.layers    <- list.files(paste0(head.path, "/pai"))
  n.h           <- length(pai.layers)
  model.heights <- seq(0, (n.h - 1) * dzd, by = dzd)
  yr.seq        <- seq(2004, 2024, 1)
  
  inpath  <- paste0(volumes3, sample.name, "/")
  
  forest.types <- unique(values(pca, na.rm = T))
  
  outpath <- paste0(out.data,sample.name,"/micro_zoned/")
  outpath <- paste0(outpath,"/all_years/")
  
  files <- list.files(outpath, 
                      pattern = paste0("VerticalSDinDailyTmax_2004to2024_"),
                      full.names = T,
                      recursive = T)
  
  require(stringr)
  df_zone <- purrr::map_dfr(files, function(f) {
    
    r <- rast(f)
    
    forest_type <- str_extract(f, "ForestType_\\d+")
    zone <- str_match(f, "VerticalSDinDailyTmax_2004to2024_(.*)\\.tif")[,2]
    
    vals <- values(r, mat = FALSE)
    
    data.frame(
      site = sample.name,
      forest_type = forest_type,
      zone = zone,
      cell = seq_along(vals),
      therm_het_zone = vals
    ) %>%
      filter(is.finite(therm_het_zone))
  })
  
  #' Generate PAI for slices.
  pad.processed <- processPAD(head.path, dzd, chm)
  pad.aboveG <- pad.processed[[1]]
  pad.aboveG <- resample(pad.aboveG,chm)
  start_heights <- pad.processed[[2]]
  z0G <- pad.processed[[3]]
  df_pai <- purrr::map_dfr(forest.types, function(ft) {
    
    zones <- readRDS(paste0(in.path, "VegZoning_FT_", ft, ".RDS"))
    
    purrr::map_dfr(seq_len(nrow(zones)), function(z) {
      
      z1 <- zones[z, ]
      
      zslice <- which(z0G >= z1$zone_start & z0G < z1$zone_end)
      
      PAI <- app(pad.aboveG[[zslice]], sum, na.rm = TRUE) * dzd
      
      # mask to this forest type
      PAI[pca != ft] <- NA
      
      vals <- values(PAI, mat = FALSE)
      
      data.frame(
        site        = sample.name,
        forest_type = paste0("ForestType_", ft),
        zone        = paste0(z1$zone_start, "to", z1$zone_end),
        cell        = seq_along(vals),
        PAI_zone    = vals,
        PAI_total = values(pad.aboveG[[1]], mat = FALSE)
      ) %>%
        filter(is.finite(PAI_zone))
    })
  })
  
  df_zone <- df_zone %>%
    left_join(
      df_pai,
      by = c("site", "forest_type", "zone", "cell")
    )
  df.list[[i]] <- df_zone
  
}

df_all <- bind_rows(df.list)

#' ---------------------------------------------------------------------------
#' ---------------------------------------------------------------------------
df_zone <- df_all %>%
  dplyr::filter(PAI_zone > 0) %>%   # remove empty vegetation
  dplyr::mutate(
    zone_id = paste(forest_type, zone, sep = "_")
  )

#' Model 1
#' Does total PAI predict thermal heterogeneity in any zone?
m_total <- lm(therm_het_zone ~ PAI_total, data = df_zone)
summary(m_total)$r.squared


#' Model 2
#' Does PAI in a zone predict thermal heterogeneity in that zone, ignoring which zone it is?
#' i.e. does vegetation amount matter overall.
m_totalZ <- lm(therm_het_zone ~ PAI_zone, data = df_zone)
summary(m_totalZ)$r.squared

#' Model 3
#' fits a different relationship for each zone_id
#' vegetation amount matters differently depending on where in the forest column that vegetation occurs.
m_zone <- lm(
  therm_het_zone ~ PAI_zone * zone_id,
  data = df_zone
)
summary(m_zone)$r.squared

AIC(m_total, m_totalZ, m_zone)

model_results <- data.frame(
  model = c(
    "Total PAI",
    "Zone PAI",
    "Zone-specific PAI"
  ),
  interpretation = c(
    "Whole-column vegetation amount",
    "Vegetation amount within matching vertical zone",
    "Zone-specific PAI–thermal heterogeneity relationship"
  ),
  R2 = c(
    summary(m_total)$r.squared,
    summary(m_totalZ)$r.squared,
    summary(m_zone)$r.squared
  ),
  AIC = AIC(m_total, m_totalZ, m_zone)$AIC
)

write_csv(model_results, paste0(out.data,"pai_models.csv"))



#' Plot...
df_plot <- df_zone %>%
  mutate(
    zmin = as.numeric(str_extract(zone, "^[0-9]+")),
    zmax = as.numeric(str_extract(zone, "(?<=to)[0-9]+")),
    zone_mid = (zmin + zmax) / 2
  ) %>%
  group_by(site, forest_type, zone_id, zone, zmin, zmax, zone_mid) %>%
  summarise(
    PAI_med   = median(PAI_zone, na.rm = TRUE),
    PAI_q25   = quantile(PAI_zone, 0.25, na.rm = TRUE),
    PAI_q75   = quantile(PAI_zone, 0.75, na.rm = TRUE),
    
    therm_med = median(therm_het_zone, na.rm = TRUE),
    therm_q25 = quantile(therm_het_zone, 0.25, na.rm = TRUE),
    therm_q75 = quantile(therm_het_zone, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(site, forest_type, zone_mid)

p <- ggplot(df_plot,
       aes(x = therm_med, y = zone_mid, colour = forest_type)) +
  geom_errorbar(
    aes(xmin = therm_q25, xmax = therm_q75),
    width = 0.4,
    alpha = 0.8
  ) +
  geom_point(aes(size = PAI_med), alpha = 0.8) +
  scale_colour_manual(
    values = c("#769c43", "#5bbda7", "#db4655", "#db9346"),
    labels = c("1", "2", "3", "4")
  )+
  scale_size_continuous(name = "Median PAI") +
  facet_wrap(~ site) +
  labs(
    y = "Midpoint height of vertical zone (m)",
    x = "Median thermal heterogeneity within vertical zone",
    colour = "Forest type"
  ) +
  theme_bw()
p

ggsave(
  paste0(out.data, "plots/ZoneHetero.png"),
  plot = p,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)


#' ---------------------------------------------------------------------------
#' ---------------------------------------------------------------------------
#' Supplementary plots.
#' ---------------------------------------------------------------------------
#' Model Diagnostics
par(mfrow = c(2,2))
plot(m_zone)
png("Figure_S4_modelD_Full.png", width = 2000, height = 1800, res = 300)

par(mfrow = c(2, 2))
plot(m_zone)
dev.off()


#' Model comparison
model_results <- bind_rows(
  glance(m_total)  %>% mutate(model = "Whole-column total PAI"),
  glance(m_totalZ) %>% mutate(model = "Zone PAI only"),
  glance(m_zone)   %>% mutate(model = "Zone-specific PAI")
) %>%
  select(model, r.squared, adj.r.squared, AIC, BIC, df.residual, p.value)
model_results


