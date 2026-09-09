#' Calc. Mean Daily Maximum Temperature.
scripts.path <- "scripts/"
source(paste0(scripts.path,"imbalanga/microParameters.R")) # loads worker functions

# -------------------------------
# Cluster batch command
# -------------------------------
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
d <- array_id + offset

# -------------------------------
in.path <- paste0(head.path,"dailyTemps/")

#' PAI layers are every 0.5 m
pai.layers <- list.files(
  paste0(head.path,"/pai"),
  pattern = "_to_canopy\\.tif$"
)
n.h <- length(pai.layers)
#' Maximum available height
max.height <- n.h * 0.5
#' Model heights: 0.2 m then every 1 m
model.heights <- c(0.2, seq(1, max.height, by = 1))
n_heights <- length(model.heights)
#' One Slurm task = one height
hm <- model.heights[d]
height <- sprintf("%02.1fm", hm)


max.files <- list.files(paste0(in.path,height,"/"),
                        pattern = "DailyMax.*.tif$",
                        full.names = T)
r <- rast(max.files)
mean.r <- mean(r, na.rm = T)

pai <- rast(paste0(head.path,"pai/PAI_0.0m_to_canopy.tif"))

clean.r <- lapply(
  mean.r,
  function(r) {
    mask(
      mean.r,
      pai[[1]],
      maskvalues = NA,
      updatevalue = NA
    )
  }
)
r <- rast(clean.r)
plot(r)

out.path <- paste0(head.path,"MeanDailyTemps/")
dir.create(out.path, showWarnings = F)
writeRaster(r, paste0(out.path, "MeanDailyMax_",height,".tif"))

