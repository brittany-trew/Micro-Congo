#' Process Land Cover.
#' https://developers.google.com/earth-engine/datasets/catalog/ESA_WorldCover_v200#bands
lc <- rast(paste0(in.data,"ESA_LandCover/WorldCover2020_",sample.name,".tif"))
lc <- project(lc, chm_mosaic)
lc <- crop(lc, chm_mosaic)
plot(lc)

pai.path <- paste0(out.data,site.location,"/",sample.name,"/",year,"/pai/")
all.pai <- list.files(pai.path, full.names = T)
pai.stk <- rast()
for(i in 1:length(all.pai)){
  pai.lidar <- rast(all.pai[[i]])
  pai.mean <- mean(pai.lidar, na.rm = T)
  pai.stk <- c(pai.stk,pai.mean)
}

pai.stk[is.na(pai.stk)] <- 0
plot(pai.stk[[1:6]])

bare_ground_mask <- (pai.stk[[1]] < 0.5) & (pai.stk[[2]]==0)
bare_ground_numeric <- classify(bare_ground_mask, cbind(0, NA), include.lowest = TRUE)
plot(bare_ground_numeric, main = "Bare Ground Mask", col = "brown")

grass_mask <- (pai.stk[[1]] >= 0.5) & (pai.stk[[2]]==0)
grass_ground_numeric <- classify(grass_mask, cbind(0, NA), include.lowest = TRUE)
plot(grass_ground_numeric, main = "Savanna Mask", col = "lightgreen")

tree_mask <- (pai.stk[[2]] > 0)
tree_mask <- classify(tree_mask, cbind(0, NA), include.lowest = TRUE)
plot(tree_mask, main = "Tree Mask", col = "darkgreen")

lc_fine <- pai.stk[[1]]
values(lc_fine) <- 0 # Forest

lc_fine[bare_ground_mask == 1] <- 16  # Bare ground
lc_fine[tree_mask == 1]        <- 2  # Evergreen broadleaf
lc_fine[grass_mask == 1]       <- 10  # short grassland

# Plot final classification
plot(lc_fine,
     main = "Land Cover Classification")
writeRaster(lc_fine, paste0(out.data,site.location,"/",sample.name,"/",year,"/landcover.tif"), overwrite = T)
