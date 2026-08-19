#' PAD vertical profiles
scripts.path <- "scripts/"
source(paste0(scripts.path,"parameters.R")) # loads worker functions

head.path <- paste0(out.data.old,site.location,"/",sample.name,"/",year,"/")

chm <- rast(paste0(head.path,"chm.tif"))

# ---- PAD files ----
all.pad <- list.files(paste0(head.path,"pad"), pattern = ".tif",full.names = TRUE)
pad.list <- lapply(all.pad, rast)
pad.r <- do.call(mosaic, pad.list)
plot(pad.r[[1:6]])

n <- nlyr(pad.r)
z0 <- (0:(n-1)) * dzd # bottom of each height (Z 0)
keep <- which(z0 >= dzd)
pad.aboveG <- pad.r[[keep]] # remove ground layer 0-0.5m
plot(pad.aboveG[[1]])

nG  <- nlyr(pad.aboveG)
z0G <- z0[keep]

chm_max <- as.numeric(global(chm, "max", na.rm=TRUE)[1,1])
max_start <- floor(chm_max / dzd) * dzd
start_heights <- z0G[z0G <= max_start]

pad.aboveG <- pad.aboveG[[1:length(start_heights)]]

# ---- Build long dataframe (safe) ----
df <- NULL
for(i in seq_len(nlyr(pad.aboveG))){
  v <- as.vector(pad.aboveG[[i]])    # <- THIS FIXES YOUR ERROR
  tmp <- data.frame(
    height_m = start_heights[i],
    pad = v
  )
  df <- rbind(df, tmp)
}

df <- na.omit(df)

# Remove NA rows (cleaner than inside values())
dff <- dplyr::filter(df, height_m < 25)
dff <- dplyr::filter(dff, pad > 0)

# ---- Mean profile ----
mean_profile <- aggregate(pad ~ height_m, dff, mean)

# grey envelope = spread of PAD at each height
band <- aggregate(
  pad ~ height_m,
  dff,
  function(x) {
    q <- quantile(x, c(0.01, 0.99), na.rm = TRUE)
    c(lo = q[1], hi = q[2])
  }
)

# Convert the matrix column into real columns immediately
band <- data.frame(
  height_m = band$height_m,
  lo = band$pad[,1],
  hi = band$pad[,2]
)
colnames(band) <- c("Height","lo","hi")

ggplot() +
  geom_ribbon(data = band,
              aes(y = Height, xmin = lo, xmax = hi),
              fill = "grey70", alpha = 0.6) +
  geom_path(data = mean_profile,
            aes(x = pad, y = height_m),
            linewidth = 0.5, colour = "black") +
  labs(x = "Plant Area Density", y = "Height (m)") +
  theme_classic()
