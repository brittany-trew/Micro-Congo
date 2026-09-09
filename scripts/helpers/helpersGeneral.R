#' Manuscript: Vertical microclimate diversity shapes thermal refugia in Congolese rainforests.
#' Date: August 2026
#' Author: Brittany Trew (*brittany.trew@fas.harvard.edu*)
#' Description: This script contains general helper functions required for running the code.
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

#' @title Load and Install Required R Packages
#' @description
#' Loads a vector of R package names. If a package is not already installed, the function installs it from CRAN and then loads it. Useful for reproducible workflows or scripts with many dependencies.
#' @param packages Character vector. Names of the packages to load.
#' @return No return value. Side effect: packages are loaded into the session (and installed if missing).
#' @export
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
      message(paste("Installing missing package:", pkg))
      install.packages(pkg, quiet = TRUE)
    }
    suppressPackageStartupMessages(
      suppressMessages(
        suppressWarnings(
          library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        )
      )
    )
  }
}

collapse_vertical <- function(files, chm) {
  r.lst <- lapply(files, load_masked, chm = chm)
  r.stk <- rast(r.lst)
  
  r.stk[r.stk < 0] <- NA
  
  n_days <- nlyr(r.lst[[1]])
  index  <- rep(seq_len(n_days), times = length(r.lst))
  
  colmax <- tapp(
    r.stk, index,
    fun = max,
    na.rm = TRUE
  )
  
  colmin <- tapp(
    r.stk, index,
    fun = min,
    na.rm = TRUE
  )
  
  colsd <- tapp(
    r.stk, index,
    fun = sd,
    na.rm = TRUE
  )
  
  col95 <- tapp(
    r.stk,
    index,
    fun = function(x, ...) {
      quantile(
        x,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
      )
    }
  )
  
  gc()
  
  return(list(
    colmax = colmax,
    colmin = colmin,
    colsd  = colsd,
    col95  = col95
  ))
}

# Extract height (metres) from parent folder name
# e.g. ".../dailyTemps/5.0m/DailyMax_2009.tif" -> 5
get_height <- function(f) {
  as.numeric(gsub("m$", "", basename(dirname(f))))
}

# Load a single height-band raster, resample to CHM grid, then mask pixels
# where canopy height is below this height band (i.e. above canopy)
load_masked <- function(f, chm) {
  hh <- get_height(f)
  r  <- resample(rast(f), chm)
  ifel(chm <= hh, NA, r)
}
