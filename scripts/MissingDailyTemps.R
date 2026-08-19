inpath <- paste0(volumes3, sample.name, "/dailyTemps")

yr.seq      <- seq(2004, 2024, 1)
height_dirs <- list.dirs(inpath, full.names = TRUE, recursive = FALSE)

missing <- list()

for (d in height_dirs) {
  ht <- basename(d)
  for (yr in yr.seq) {
    f <- file.path(d, paste0("DailyMax_", yr, ".tif"))
    if (!file.exists(f)) {
      missing[[length(missing) + 1]] <- data.frame(height = ht, year = yr)
    }
  }
}

if (length(missing) == 0) {
  message("All files present — ", length(height_dirs), " height bands x ", length(yr.seq), " years.")
} else {
  missing_df <- do.call(rbind, missing)
  message(nrow(missing_df), " files missing:")
  print(missing_df)
}
