#' Summary Statistics.
all.map1 <- list.files(out.data, 
                       pattern = "Map1_VerticalSDinDailyTmax_2004to2024.tif",
                       recursive = T,
                       full.names = T)

site.stats <- data.frame()
for(i in 1:length(all.map1)){
  
  name <- unlist(strsplit(all.map1[[i]], "/"))[[3]]
  r <- rast(all.map1[[i]])
  
  rv <- values(r, na.rm = TRUE)
  
  tmp <- data.frame(
    name   = name,
    median = round(median(rv),1),
    IQR_25 = round(quantile(rv, 0.25),1),
    IQR_75 = round(quantile(rv, 0.75),1),
    Q90    = round(quantile(rv, 0.90),1)
  )
  
  site.stats <- rbind(site.stats, tmp)
}

site.stats <- site.stats
write_csv(site.stats,paste0(out.data,"verticalsummaryMap1.csv"))


all.map4 <- list.files(out.data, 
                       pattern = "Map4_VerticalMeanofSDinDailyTmax_2004to2024.tif",
                       recursive = T,
                       full.names = T)

site.stats <- data.frame()
for(i in 1:length(all.map4)){
  
  name <- unlist(strsplit(all.map4[[i]], "/"))[[3]]
  r <- rast(all.map4[[i]])
  
  rv <- values(r, na.rm = TRUE)
  
  tmp <- data.frame(
    name   = name,
    median = round(median(rv),1),
    IQR_25 = round(quantile(rv, 0.25),1),
    IQR_75 = round(quantile(rv, 0.75),1),
    Q90    = round(quantile(rv, 0.90),1)
  )
  
  site.stats <- rbind(site.stats, tmp)
}

site.stats
write_csv(site.stats,paste0(out.data,"verticalsummaryMap4.csv"))
