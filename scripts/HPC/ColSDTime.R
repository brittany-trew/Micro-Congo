#' Calculate temporal heterogeneity in daily maximum temperature.
scripts.path <- "scripts/"
source(paste0(scripts.path,
              "imbalanga/microParameters.R"))

# -------------------------------
# Cluster batch command
# -------------------------------
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
offset <- as.numeric(Sys.getenv("ARRAY_OFFSET", "0"))
d <- array_id + offset


# -------------------------------
# Setup
# -------------------------------
chm <- rast(paste0(head.path,
                   "chm.tif"))

inpath <- paste0(head.path,
                 "dailyTemps/")

outpath <- paste0(head.path,
                  "verticalSummary/")

dir.create(outpath,
           recursive = TRUE,
           showWarnings = FALSE)

# Identify the number of modelled heights
pai.layers <- list.files(paste0(head.path,
                                "pai/"),
                         pattern = "\\.tif$",
                         full.names = TRUE)

n.h <- length(pai.layers)

model.heights <- c(0.2,
                   seq(1, n.h - 1, by = 1))

# Check that the array index represents a valid height
if (d < 1 || d > length(model.heights)) {
  stop("Array index ", d,
       " is outside the available height range 1-",
       length(model.heights), ".")
}


# -------------------------------
# Prepare working grid
# -------------------------------

# Use the first temperature raster as the working grid
template.files <- sort(list.files(
  paste0(inpath,
         "0.2m/"),
  pattern = "^DailyMax_[0-9]{4}\\.tif$",
  full.names = TRUE
))

if (length(template.files) == 0) {
  stop("No temperature rasters found in the 0.2 m directory.")
}

template <- rast(template.files[[1]])

# Resample the CHM once to the temperature grid
chm_temp <- resample(chm,
                     template)


# -------------------------------
# Process array height
# -------------------------------
hh <- model.heights[[d]]
h <- sprintf("%.1f", hh)

message("Processing height ", h, " m ...")

files_max <- sort(list.files(
  paste0(inpath,
         h, "m/"),
  pattern = "^DailyMax_[0-9]{4}\\.tif$",
  full.names = TRUE
))

if (length(files_max) == 0) {
  stop("No DailyMax files found for height ", h, " m.")
}

message("Found ", length(files_max),
        " yearly files for ", h, " m.")

# Stack every daily Tmax layer from 2004–2024
r.stk <- rast(files_max)

# Calculate temporal SD across all days
# Negative temperatures are treated as invalid model values
gridsd <- app(
  r.stk,
  fun = function(x) {
    x <- x[is.finite(x) & x >= 0]
    
    if (length(x) < 2) {
      return(NA_real_)
    }
    
    sd(x)
  }
)

# Remove heights above the local canopy
gridsd <- ifel(chm_temp <= hh,
               NA,
               gridsd)

names(gridsd) <- paste0("temporalSD_",
                        h, "m")

# Save the result for this height
outfile <- paste0(outpath,
                  "temporalSD_",
                  h,
                  "m.tif")

writeRaster(gridsd,
            outfile,
            overwrite = TRUE)

message("Output written for height ",
        h, " m: ",
        outfile)
