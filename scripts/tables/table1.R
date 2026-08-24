#' Table 1: Summary of results for Forest Types.

# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

all.samples <- all.samples
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
  in.path <- paste0(out.data,sample.name,"/veg_zoning/")
  chm           <- rast(paste0(head.path, "/chm.tif"))
  pai <- mean(rast(paste0(head.path, "pai/pai_0.5m_to_canopy.tif")), na.rm = T)
  
  pad.processed <- processPAD(head.path, dzd, chm)
  pad.aboveG <- pad.processed[[1]]
  pad.aboveG <- resample(pad.aboveG, chm)
  pad_heights <- c(0.5,seq(1, (nlyr(pad.aboveG) - 1) * dzd, by = dzd))
  
  peak_pad_q95 <- app(pad.aboveG, fun = function(x) {
    if (all(is.na(x))) return(NA_real_)
    quantile(x, 0.95, na.rm = TRUE)
  })
  peak_height <- app(pad.aboveG, fun = function(x) {
    if (all(is.na(x))) return(NA_real_)
    pad_heights[which.max(x)]
  })
  
  sdPAD <- terra::app(
    pad.aboveG,
    fun = function(x) {
      if (sum(is.finite(x)) < 2) return(NA_real_)
      sd(x, na.rm = TRUE)
    }
  )  
  
  #' Forest Type Area & cells#' Foreschmt Type Area & cells
  pca.path <- paste0(out.data,"pca_clusters/")
  pca <- rast(paste0(pca.path, sample.name, "_pca_Clust4.tif"))
  
  cell_area_m2 <- prod(res(pca))
  forest_area <- terra::freq(pca) %>%
    as.data.frame() %>%
    rename(
      forest_type = value,
      n_cells = count
    ) %>%
    mutate(
      area_m2 = n_cells * cell_area_m2,
      area_ha = area_m2 / 10000,
      percent = 100 * area_m2 / sum(area_m2)
    ) %>%
    select(
      forest_type, n_cells, area_ha, percent
    )
  
  summary_df <- forest_area
  forest.types <- forest_area$forest_type
  
  zone_summary_list <- list()
  for (f in seq_along(forest.types)) {
    
    ft <- forest.types[f]
    zones <- readRDS(paste0(in.path, "VegZoning_FT_", ft, ".RDS"))
    
    #' CHM
    chm_ft <- terra::ifel(pca == ft, chm, NA)
    median_chm <- terra::global(
      chm_ft,
      FUN = "median",
      na.rm = TRUE
    )[1, 1]
    
    #' PAI
    pai_ft <- terra::ifel(pca == ft, pai, NA)
    median_pai <- terra::global(
      pai_ft,
      FUN = "median",
      na.rm = TRUE
    )[1, 1]
    
    #' PAD.
    peak_ft <- terra::ifel(pca == ft, peak_pad_q95, NA)
    median_padPeak <- terra::global(
      peak_ft,
      FUN = "median",
      na.rm = TRUE
    )[1, 1]
    
    heightPAD_ft <- terra::ifel(pca == ft, peak_height, NA)
    median_padPeakH <- terra::global(
      heightPAD_ft,
      FUN = "median",
      na.rm = TRUE
    )[1, 1]
    
    sdPAD_ft <- terra::ifel(pca == ft, sdPAD, NA)
    median_sdPAD <- terra::global(
      sdPAD_ft,
      FUN = "median",
      na.rm = TRUE
    )[1, 1]
    
    zone_stats <- zones %>%
      mutate(
        zone_mid = (zone_start + zone_end) / 2,
        zone_height = zone_end - zone_start
      ) %>%
      summarise(
        forest_type = ft,
        n_zones = n(),
        mean_zone_height = mean(zone_height, na.rm = TRUE),
        median_chm = median_chm,
        median_pai = median_pai,
        median_Max_PAD = median_padPeak,
        median_Max_PAD_Height = median_padPeakH,
        median_SD_PAD = median_sdPAD
      )
    zone_summary_list[[f]] <- zone_stats
  }
  
  zone_summary <- bind_rows(zone_summary_list)
  summary_df <- summary_df %>%
    left_join(zone_summary, by = "forest_type")
  summary_df

  summary_df <- summary_df %>%
    mutate(site = sample.name)
  
  site_list[[sample.name]] <- summary_df
  }  

site_summaries <- bind_rows(site_list)

summary_all_sites <- site_summaries %>%
  group_by(forest_type) %>%
  summarise(
    n_sites = n_distinct(site),
    
    # Area: summed across sites
    total_n_cells = sum(n_cells, na.rm = TRUE),
    total_area_ha = sum(area_ha, na.rm = TRUE),
    
    # Zones: site-level summaries
    mean_n_zones = mean(n_zones, na.rm = TRUE),
    sd_n_zones = sd(n_zones, na.rm = TRUE),
    
    mean_zone_height = mean(mean_zone_height, na.rm = TRUE),
    sd_zone_height = sd(mean_zone_height, na.rm = TRUE),
    
    # Structure: mean of site-level medians
    mean_median_chm = mean(median_chm, na.rm = TRUE),
    sd_median_chm = sd(median_chm, na.rm = TRUE),
    
    mean_median_pai = mean(median_pai, na.rm = TRUE),
    sd_median_pai = sd(median_pai, na.rm = TRUE),
    
    mean_median_Max_PAD = mean(median_Max_PAD, na.rm = TRUE),
    sd_median_Max_PAD = sd(median_Max_PAD, na.rm = TRUE),
    
    mean_median_Max_PAD_Height = mean(median_Max_PAD_Height, na.rm = TRUE),
    sd_median_Max_PAD_Height = sd(median_Max_PAD_Height, na.rm = TRUE),
    
    mean_median_SD_PAD = mean(median_SD_PAD, na.rm = TRUE),
    sd_median_SD_PAD = sd(median_SD_PAD, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    total_percent = 100 * total_area_ha / sum(total_area_ha, na.rm = TRUE)
  )

summary_all_sites
zone_height_weighted <- site_summaries %>%
  group_by(forest_type) %>%
  summarise(
    weighted_mean_zone_height =
      sum(mean_zone_height * n_zones, na.rm = TRUE) /
      sum(n_zones, na.rm = TRUE),
    total_zone_count = sum(n_zones, na.rm = TRUE),
    .groups = "drop"
  )
summary_all_sites <- summary_all_sites %>%
  left_join(zone_height_weighted, by = "forest_type")

summary_all_sites_clean <- summary_all_sites %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  )

summary_all_sites_clean
write_csv(
  summary_all_sites_clean,
  paste0(out.data, "tables/forest_type_structure_summary_all_sites.csv")
)
