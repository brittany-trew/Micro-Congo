# ── Setup ---------------------------------------------------------------------
scripts.path <- "scripts/"
source(paste0(scripts.path, "parameters.R"))

all.samples <- all.samples[-1]

#' Vertical SD versus PAI.
site.list1 <- list()
site.list2 <- list()
for(i in 1:length(all.samples)){
  
  sample.name <- all.samples[i]
  head.path <- paste0(out.data.old,site.location,"/",sample.name,"/",year,"/")
  pai <- rast(paste0(head.path,"pai/PAI_0.5m_to_canopy.tif"))
  mean.pai <- mean(pai)
  mean.pai[mean.pai==0]<- NA
  plot(mean.pai)
  
  inpath <- paste0(out.data, sample.name, "/Vertical_Summaries/")
  
  r <- rast(paste0(inpath, "Map1_VerticalSDinDailyTmax_2004to2024.tif"))
  r <- resample(r, mean.pai)
  df <- data.frame(
    values(mean.pai),
    values(r)
  )
  df <- na.omit(df)
  colnames(df) <- c("pai", "therm_het")
  site.list1[[i]] <- df
  
  r <- rast(paste0(inpath, "Map4_VerticalMeanofSDinDailyTmax_2004to2024.tif"))
  r <- resample(r, mean.pai)
  df <- data.frame(
    values(mean.pai),
    values(r)
  )
  df <- na.omit(df)
  colnames(df) <- c("pai", "therm_het")
  
  site.list2[[i]] <- df
}

df1 <- dplyr::bind_rows(site.list1)
df2 <- dplyr::bind_rows(site.list2)

#' (1) Vertical
m_total <- lm(therm_het ~ pai, data = df1)
r2 <- summary(m_total)$r.squared

p <- ggplot(df1, aes(x = pai, y = therm_het)) +
  geom_point(alpha = 0.05, size = 0.2) +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 5), linewidth = 0.7, colour = "#db4655") +
  labs(
    x = "Total plant area index",
    y = "Mean vertical SD in daily maximum temperatures"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("R² = ", round(r2, 3)),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  )+
  theme_bw()
p
ggsave(
  paste0(out.data, "plots/TotalPAI_MeanVerticalSD.png"),
  plot = p,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

#' (1) Temporal
m_total <- lm(therm_het ~ pai, data = df2)
r2 <- summary(m_total)$r.squared

p <- ggplot(df2, aes(x = pai, y = therm_het)) +
  geom_point(alpha = 0.05, size = 0.2) +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 5), linewidth = 0.7, colour = "#db4655") +
  labs(
    x = "Total plant area index",
    y = "Vertical mean of SD in daily maximum temperatures"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("R² = ", round(r2, 3)),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  )+
  theme_bw()
p
ggsave(
  paste0(out.data, "plots/TotalPAI_TimeSD.png"),
  plot = p,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)
