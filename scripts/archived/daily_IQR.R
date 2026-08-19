#' Calculate IQR of Tmean and DTR.
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R"))

dzd =1
sample.name <- all.samples[1]
head.path <- paste0(out.data.old,site.location,"/",sample.name,"/",year,"/")

pai.layers <- list.files(paste0(head.path,"/pai"))
n.h <- length(pai.layers)
model.heights <- c(0.2,seq(dzd,(n.h-1)*dzd,by = dzd))

chm <- rast(paste0(head.path,"/chm.tif"))

outpath1 <- paste0(out.data,"SummaryRasters/",sample.name,"/")
dir.create(outpath1)
outpath2 <- paste0(out.data,"DTR_sliced/",sample.name,"/")
dir.create(outpath2)

for(h in 1:length(model.heights)){
  height <- sprintf("%02.1fm", model.heights[h])
  
  inpath <- paste0(volumes3,sample.name,"/dailyTemps/")
  inpath <- paste0(inpath,height,"/")
  
  ## IQR DTR (Tmax - Tmin)
  files_tmax <- list.files(inpath, pattern = "DailyMax_.*tif$", full.names = TRUE)
  if(length(files_tmax)!=21){stop("Tmax Missing years")}
  files_tmin <- list.files(inpath, pattern = "DailyMin_.*tif$", full.names = TRUE)
  if(length(files_tmin)!=21){stop("Tmin Missing years")}
  
  tmax_full <- rast(files_tmax)
  tmin_full <- rast(files_tmin)

  # tmax_mean <- mean(tmax_full)
  # plot(tmax_mean)
  # writeRaster(tmax_mean, paste0(outpath1,"Mean_Tmax_2004_2024_",height,".tif"),
  #             overwrite = TRUE)
  
  #' Max of Tmax.
  #' 95th percentile temperature:
  robust_max <- app(tmax_full, fun = function(x) quantile(x, 0.95, na.rm = TRUE))
  plot(robust_max)
  writeRaster(robust_max, paste0(outpath1,"P95_Tmax_2004_2024_",height,".tif"),
              overwrite = TRUE)
  tmax_median <- app(tmax_full, fun = median, na.rm = TRUE)
  plot(tmax_median)
  writeRaster(tmax_median, paste0(outpath1,"Median_Tmax_2004_2024_",height,".tif"),
              overwrite = TRUE)
  
  #' Daily Diurnal Temperature Range
  dtr_full <- tmax_full - tmin_full
  # Central tendency of daily temperature range:
  dtr_median <- app(dtr_full, fun = median, na.rm = TRUE)
  plot(dtr_median)
  writeRaster(dtr_median, paste0(outpath2,"DTR_Median_2004_2024_",height,".tif"),
              overwrite = TRUE)
  
  dtr_tmax_cor <- app(c(tmax_full, dtr_full), fun = function(x) {
    n <- length(x) / 2
    tx <- x[1:n]
    dt <- x[(n+1):(2*n)]
    
    # Return NA if too few complete pairs
    complete <- sum(!is.na(tx) & !is.na(dt))
    if (complete < 30) return(NA)
    cor(tx, dt, use = "complete.obs")
  })
  writeRaster(dtr_tmax_cor, paste0(outpath2,"DTR_TmaxCOR_2004_2024_",height,".tif"),
              overwrite = TRUE)
}




#' Extra code. 
#' -----------------------------------------
df <- as.data.frame(c(tmax_median, dtr_median, dtr_tmax_cor), na.rm = TRUE)
names(df) <- c("tmax", "dtr", "cor")

# Is correlation driven by tmax, DTR, or both?
ggplot(df, aes(x = tmax, y = dtr, colour = cor)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_colour_gradient2(low = "black", mid = "lightgreen", high = "pink", midpoint = median(df$cor)) +
  labs(x = "Median Tmax (°C)", y = "Median DTR (°C)", colour = "Correlation")+
  xlim(25,45)+
  ylim(5,20)




outpath <- paste0(head.path,"IQR/")
files <- list.files(outpath, pattern = "IQR_DTR_.*tif$", full.names = TRUE)
# Extract heights from filenames like "IQR_TMean_2004_2024_05.0m.tif"
heights <- as.numeric(gsub(".*_(\\d+\\.\\d+)m\\.tif$", "\\1", files))
# Load into raster list
r_list <- lapply(files, rast)
rstk <- rast(r_list)
names(rstk) <- heights

vals <- terra::values(rstk)      # matrix: rows = pixels, columns = layers
vals_df <- as.data.frame(vals)
vals_df <- vals_df %>% tidyr::drop_na()

# Convert to long format for ggplot
vals_long <- vals_df %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Height",
                      values_to = "IQR")
library(ggplot2)

tropical_cols <- c(
  "#50ADB5", "#D99000", "#2C6A73", "#D43F6E", "#99A170",
  "#4F6B32", "#2E4A1F", "#358C78", "#E06B00", "#89A35D",
  "#B54B87", "#446D5A", "#7B5CA8", "#CB6E33", "#3A5130",
  "#343E75", "#9FB567", "#8A1F3B", "#5C5A58", "#16174A",
  "black"
)

vals_long$Height <- as.numeric(as.character(vals_long$Height))

p <- vals_long %>% 
  ggplot(aes(x = factor(Height), y = IQR, fill = as.factor(Height))) +
  # add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    # adjust bandwidth
    adjust = 0.5,
    # move to the right
    justification = -0.2,
    # remove the slub interval
    .width = 0,
    point_colour = NA
  )+
  geom_boxplot(
    width = 0.3,
    # removing outliers
    alpha = 0,
  )+
  scale_fill_manual(values = tropical_cols) + # Use custom color scale for fill
  scale_color_manual(values = tropical_cols)+
  labs(
    x = "Height Above Ground (m)",
    y = "Interquantile Range of Daily Mean Temperature"
  ) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, colour = "darkred")+
  coord_flip()+
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.text = element_text(size = 10))

plot(p)
ggsave(paste0(head.path,"figures/IQR_heights_BoxPlot.png"), plot = p, width = 7, height = 8, dpi = 300)



#' Map range of daily mean and daily range for each vertical column?
## IQR TMEAN.
files <- list.files(inpath, pattern = "DailyMean_.*tif$", full.names = TRUE)
tmean_full <- rast(files)   # this will be (365*years) layers

iqr_full <- app(tmean_full, IQR, na.rm = TRUE)

plot(iqr_full, main = "IQR of Tmean")

outpath <- paste0(head.path,"IQR/")
dir.create(outpath, showWarnings = F)
writeRaster(iqr_full, paste0(outpath,"IQR_TMean_2004_2024_",height,".tif"), 
            overwrite = TRUE)
