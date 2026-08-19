#' This script conducts a PCA on ecologically significant PAI layers to map different vegetation structures within the forest.
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R")) # loads worker functions

df <- tibble()
X_list <- list()
keep_list <- list()
template_list <- list()
for(i in 1:length(all.samples)){
  
  #' Set datastream.
  sample.name <- all.samples[i]
  head.path <- paste0(out.data.old,site.location,"/",sample.name,"/",year,"/")
  chm <- rast(paste0(head.path,"chm.tif"))
  
  pad.processed <- processPAD(head.path, dzd, chm)
  pad.aboveG <- pad.processed[[1]]
  z0G <- pad.processed[[2]]
  
  # Identify slices for combining.
  idx_0_2  <- which(z0G >= 0.5 & z0G < 2)
  idx_2_10 <- which(z0G >= 2   & z0G < 10)
  idx_10_25 <- which(z0G >= 10   & z0G < 20)
  
  # Integrate vertically
  PAI_0_2  <- app(pad.aboveG[[idx_0_2]],  sum, na.rm = TRUE) * dzd
  PAI_2_10 <- app(pad.aboveG[[idx_2_10]], sum, na.rm = TRUE) * dzd
  PAI_10_25 <- app(pad.aboveG[[idx_10_25]], sum, na.rm = TRUE) * dzd
  
  pad.25_0to2 <- aggregate(PAI_0_2, 5, fun = "mean", na.rm = T)
  pad.25_2to10 <- aggregate(PAI_2_10, 5, fun = "mean", na.rm = T)
  pad.25_10to25 <- aggregate(PAI_10_25, 5, fun = "mean", na.rm = T)

  pai_25 <- c(pad.25_0to2, pad.25_2to10, pad.25_10to25)
  names(pai_25) <- c("PAI_0_2", "PAI_2_10", "PAI_10_25")
  plot(pai_25)
  
  X <- terra::values(pai_25, mat = TRUE)
  keep <- rowSums(X) > 0 & complete.cases(X)
  X <- X[keep, , drop = FALSE]
  
  X_list[[sample.name]] <- X
  keep_list[[sample.name]] <- keep
  template_list[[sample.name]] <- pai_25[[1]]
  
  Site <- rep(sample.name, nrow(X))
  X_df <- cbind.data.frame(X, Site)
  df <- rbind(df, X_df)
}

colnames(df) <- c("Understory", "Midstory", "Canopy", "Site Name")
#' run PCA on entire dataset.
pca <- prcomp(df[, c("Understory", "Midstory", "Canopy")], center = TRUE, scale. = TRUE)
pc_scores <- as.data.frame(pca$x[, 1:2])
pc_scores$site <- df$Site

loadings <- as.data.frame(pca$rotation[,1:2])
loadings$var <- rownames(loadings)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_exp[1] * 100, 1), "% variance)")
pc2_lab <- paste0("PC2 (", round(var_exp[2] * 100, 1), "% variance)")


pca_table <- loadings %>%
  mutate(
    PC1 = round(PC1, 3),
    PC2 = round(PC2, 3),
    Stratum = var
  ) %>%
  select(Stratum, PC1, PC2)

variance_row <- data.frame(
  Stratum = "Variance explained (%)",
  PC1 = round(var_exp[1] * 100, 1),
  PC2 = round(var_exp[2] * 100, 1)
)
pca_table <- bind_rows(pca_table, variance_row)
pca_table
write_csv(pca_table, paste0(out.data, "tables/pca_results.csv"))

p <- ggplot(pc_scores, aes(PC1, PC2, colour = site)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
    arrow = arrow(length = unit(0.2,"cm")),
    colour = "black",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = loadings,
    aes(x = PC1*4.8, y = PC2*3.4, label = var),
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("#769c43", "#5bbda7", "#db4655", "#db9346"))+
  labs(
    x=pc1_lab,
    y=pc2_lab,
    colour = "Study Site"
  ) +
  ylim(-4,4.5)+
  theme_classic()
p
ggsave(
  paste0(out.data, "plots/pca_bysite.png"),
  plot = p,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)

#' Check K clustering
library(mclust)
plot_list <- list()
for (k in 2:6) {
  gmm <- Mclust(pc_scores[, c("PC1", "PC2")], G = k)
  tmp <- pc_scores
  tmp$cluster <- factor(gmm$classification)
  
  plot_list[[as.character(k)]] <- ggplot(tmp, aes(PC1, PC2, colour = cluster)) +
    geom_point(alpha = 0.12, size = 0.35) +
    theme_classic() +
    labs(
      title = paste0("k = ", k),
      x = pc1_lab,
      y = pc2_lab
    )
}
patchwork::wrap_plots(plot_list, ncol = 3)

k_options <- c(3, 4, 5)
kmeans_list <- list()

cluster_palette <- function(k) {
  cols <- unname(site_palette)
  cols[1:k]
}

for(i in 1:length(k_options)){
  best_k <- k_options[i]
  
  cl <- Mclust(pc_scores[, c("PC1", "PC2")], G = best_k)
  kmeans_list[[paste0("k", best_k)]] <- cl$classification
  
  tmp <- pc_scores
  tmp$cluster <- factor(cl$classification)  # ← was cl$cluster
  
  centroids <- aggregate(cbind(PC1, PC2) ~ cluster, data = tmp, FUN = mean)
  
  p<- ggplot(tmp, aes(PC1, PC2, colour = cluster)) +
      geom_point(alpha = 0.2, size = 0.2) +
      stat_ellipse(aes(group = cluster), linewidth = 0.5, show.legend = FALSE) +
      geom_point(
        data = centroids,
        aes(PC1, PC2, fill = cluster),
        shape = 21, colour = "black", size = 3.5, stroke = 0.2,
        inherit.aes = FALSE
      ) +
      geom_text(
        data = centroids,
        aes(PC1, PC2, label = cluster),
        fontface = "bold", colour = "black", nudge_y = 0.3,
        show.legend = FALSE
      ) +
      scale_colour_manual(values = cluster_palette(best_k)) +
      scale_fill_manual(values = cluster_palette(best_k)) +
      ylim(-4, 4.5) +
      theme_classic() +
      labs(
        title = paste0("Structural clusters in PCA space (k = ", best_k, ")"),
        x = pc1_lab,
        y = pc2_lab,
        colour = "Structural class",
        fill = "Structural class"
      )
  plot(p)
  ggsave(plot = p, filename = paste0(out.data,"plots/pca_byStructure_",best_k,".png"))
  
}

for(i in 1:length(all.samples)){
  cluster_id <- kmeans_list[[2]]
  sample.name <- all.samples[i]
  head.path <- paste0(out.data.old,site.location,"/",sample.name,"/",year,"/")
  chm <- rast(paste0(head.path,"chm.tif"))
  
  site_rows <- df$Site == sample.name
  site_clusters <- cluster_id[site_rows]
  cluster_r <- rast(template_list[[sample.name]])
  v <- rep(NA_integer_, ncell(cluster_r))
  v[keep_list[[sample.name]]] <- site_clusters
  values(cluster_r) <- v
  cluster.final <- resample(cluster_r, chm, method = "near")
  plot(cluster.final)
  out.path <- paste0(out.data,"pca_clusters/")
  dir.create(out.path)
  writeRaster(cluster.final, paste0(out.path,sample.name,"_pca_Clust4.tif"), overwrite = T)
}


for(i in 1:length(all.samples)){
  cluster_id <- kmeans_list[[3]]  # ← remove $cluster
  sample.name <- all.samples[i]
  head.path <- paste0(out.data.old, site.location, "/", sample.name, "/", year, "/")
  chm <- rast(paste0(head.path, "chm.tif"))
  
  site_rows <- df$Site == sample.name
  site_clusters <- cluster_id[site_rows]
  cluster_r <- rast(template_list[[sample.name]])
  v <- rep(NA_integer_, ncell(cluster_r))
  v[keep_list[[sample.name]]] <- site_clusters
  values(cluster_r) <- v
  cluster.final <- resample(cluster_r, chm, method = "near")
  plot(cluster.final)
  out.path <- paste0(out.data, "pca_clusters/")
  dir.create(out.path, showWarnings = FALSE)  # ← avoids repeated warnings if folder exists
  writeRaster(cluster.final, paste0(out.path, sample.name, "_pca_Clust5.tif"),
              overwrite = TRUE)  # ← avoids error on re-runs
}

