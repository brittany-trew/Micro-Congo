#' This script conducts a PCA on ecologically significant PAI layers to map different vegetation structures within the forest.
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R")) # loads worker functions

#' [Principal Component Analysis on cumulative PAD]
df_list <- list()
X_list <- list()
keep_list <- list()
template_list <- list()
for(i in 1:length(all.samples)){
  #' Set datastream.
  sample.name <- all.samples[i]
  head.path <- paste0(mclidar.out,site.location,"/",sample.name,"/",year,"/")
  chm <- rast(paste0(head.path,"chm.tif"))
  
  all.pad <- list.files(paste0(head.path,"pad"), pattern = ".tif", full.names = T)
  pad.r <- lapply(all.pad, rast)
  pad0 <- processPAD(pad.r, chm)
  start_heights <- seq(
    from = 0,
    by = dzd,
    length.out = nlyr(pad0)
  )
  
  #' Identify slices for combining
  idx_0_2  <- which(start_heights < 2)
  idx_2_10 <- which(start_heights >= 2   & start_heights < 10)
  idx_10_Top <- which(start_heights >= 10)
  
  #' Sum slices to convert to PAI
  PAI_0_2  <- app(pad0[[idx_0_2]],  sum, na.rm = TRUE) * dzd
  PAI_2_10 <- app(pad0[[idx_2_10]], sum, na.rm = TRUE) * dzd
  PAI_10_top <- app(pad0[[idx_10_Top]], sum, na.rm = TRUE) * dzd
  
  #' Aggregate to a coarser resolution
  pad.25_0to2 <- aggregate(PAI_0_2, 5, fun = "mean", na.rm = T)
  pad.25_2to10 <- aggregate(PAI_2_10, 5, fun = "mean", na.rm = T)
  pad.25_10totop <- aggregate(PAI_10_top, 5, fun = "mean", na.rm = T)

  pai_25 <- c(pad.25_0to2, pad.25_2to10, pad.25_10totop)
  names(pai_25) <- c("PAI_0_2", "PAI_2_10", "PAI_10_Top")
  plot(pai_25)
  
  X <- terra::values(pai_25, mat = TRUE)
  keep <- rowSums(X) > 0 & complete.cases(X)
  X <- X[keep, , drop = FALSE]
  
  #' Save PAI data
  X_list[[sample.name]] <- X
  #' Save raster cell numbers for turning data back into maps
  keep_list[[sample.name]] <- keep
  #' Tamplate for maps
  template_list[[sample.name]] <- pai_25[[1]]
  #' Creates site label for all data points
  Site <- rep(sample.name, nrow(X))
  X_df <- cbind.data.frame(X, Site)
  df_list[[sample.name]] <- X_df
}

df <- bind_rows(df_list)
colnames(df) <- c("Understory", "Midstory", "Canopy", "Site")

#' run PCA on entire dataset.
pca <- prcomp(df[, c("Understory", "Midstory", "Canopy")], center = TRUE, scale. = TRUE)
pc_scores <- as.data.frame(pca$x)
pc_scores$site <- df$Site

#' Extract loadings: 
#' i.e. How much, and in which direction, each of the variables contributes to each principal component
loadings <- as.data.frame(pca$rotation)
loadings$var <- rownames(loadings)

#' Calc. how much of the structural variation each of these axes actually captures:
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_exp[1] * 100, 1), "% variance)")
pc2_lab <- paste0("PC2 (", round(var_exp[2] * 100, 1), "% variance)")
pc3_lab <- paste0("PC3 (", round(var_exp[3] * 100, 1), "% variance)")

pca_table <- loadings %>%
  mutate(
    PC1 = round(PC1, 3),
    PC2 = round(PC2, 3),
    PC3 = round(PC3, 3),
    Stratum = var
  ) %>%
  dplyr::select(Stratum, PC1, PC2, PC3)

variance_row <- data.frame(
  Stratum = "Variance explained (%)",
  PC1 = round(var_exp[1] * 100, 1),
  PC2 = round(var_exp[2] * 100, 1),
  PC3 = round(var_exp[3] * 100, 1)
)
pca_table <- bind_rows(pca_table, variance_row)
pca_table
write_csv(pca_table, paste0(out.data, "tables/pca_results.csv"))

p12 <- ggplot(pc_scores, aes(PC1, PC2)) +
  geom_point(
    alpha = 0.15,
    size = 0.4,
    colour = "#555555"
  ) +
  geom_segment(
    data = loadings,
    aes(
      x = 0,
      y = 0,
      xend = PC1 * 3,
      yend = PC2 * 3
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.4,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = loadings,
    aes(
      x = PC1 * 3.5,
      y = PC2 * 3.5,
      label = var
    ),
    inherit.aes = FALSE
  ) +
  labs(
    x = "PC1 (44.7% variance)",
    y = "PC2 (34.5% variance)"
  ) +
  theme_classic()
p12

ggsave(
  paste0(out.data, "plots/pca_pad.png"),
  plot = p12,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)


#' [Clustering.]
#' Four structural classes were retained as a parsimonious representation of the major differences 
#' in vertical vegetation structure, 
#' while avoiding subdivision into increasingly similar structural configurations.
gmm4 <- Mclust(
  pc_scores[, c("PC1", "PC2", "PC3")],
  G = 4
)
df$cluster4 <- factor(gmm4$classification)

cluster_summary <- df %>%
  dplyr::group_by(cluster4) %>%
  dplyr::summarise(
    n = n(),
    
    Understory_mean = mean(Understory),
    Understory_SD = stats::sd(Understory),
    
    Midstory_mean = mean(Midstory),
    Midstory_SD = stats::sd(Midstory),
    
    Canopy_mean = mean(Canopy),
    Canopy_SD = stats::sd(Canopy)
  )
write_csv(cluster_summary, paste0(out.data, "tables/cluster_summary.csv"))


#' Fitting the final 4-class Gaussian mixture model, 
#' assigning every PCA point to one of those four structural classes, 
#' then plotting the classes in PC1–PC2 space with ellipses and labelled centroids
tmp <- pc_scores
tmp$cluster <- factor(gmm4$classification)

centroids <- aggregate(
  cbind(PC1, PC2) ~ cluster,
  data = tmp,
  FUN = mean
)

p <- ggplot(tmp, aes(PC1, PC2, colour = cluster)) +
  geom_point(
    alpha = 0.2,
    size = 0.2
  ) +
  stat_ellipse(
    aes(group = cluster),
    linewidth = 0.5
  ) +
  geom_point(
    data = centroids,
    aes(PC1, PC2, fill = cluster),
    shape = 21,
    colour = "black",
    size = 3.5,
    stroke = 0.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = centroids,
    aes(PC1, PC2, label = cluster),
    fontface = "bold",
    colour = "black",
    nudge_y = 0.3,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = colour.p) +
  scale_fill_manual(values = colour.p) +
  coord_cartesian(
    ylim = c(-4, 4.5)
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "Structural classes in PCA space",
    x = pc1_lab,
    y = pc2_lab
  )

p
ggsave(plot = p, filename = paste0(out.data,"plots/Classes_PCA.png"))


for(i in 1:length(all.samples)){
  sample.name <- all.samples[i]
  head.path <- paste0(mclidar.out,site.location,"/",sample.name,"/", year,"/")
  
  chm <- rast(paste0(head.path, "chm.tif"))
  
  #' Identify rows belonging to this site
  site_rows <- df$Site == sample.name
  #' Get structural class for those cells
  site_clusters <- df$cluster4[site_rows]
  #' Create raster template
  cluster_r <- template_list[[sample.name]]
  #' Empty raster values
  v <- rep(NA_integer_, ncell(cluster_r))
  #' Put classifications back into original retained cells
  v[keep_list[[sample.name]]] <- as.integer(as.character(site_clusters))
  values(cluster_r) <- v
  
  #' Resample back to CHM resolution
  cluster.final <- resample(cluster_r,
    chm,
    method = "near")
  
  plot(cluster.final)
  
  out.path <- paste0(out.data,"pca_clusters/")
  dir.create(out.path,
    recursive = TRUE,
    showWarnings = FALSE)
  
  writeRaster(cluster.final,
    paste0(out.path,sample.name,"_pca_C4.tif"),
    overwrite = TRUE)
}


