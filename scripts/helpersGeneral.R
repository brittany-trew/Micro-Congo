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