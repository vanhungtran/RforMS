# _common.R — Global book configuration and AI4SciFig theme setup
# Sourced automatically or in setup chunks for consistent publication styling

local_renv <- file.path(getwd(), "renv/library/windows/R-4.6/x86_64-w64-mingw32")
if (dir.exists(local_renv)) {
  .libPaths(unique(c(local_renv, .libPaths())))
}

if (dir.exists("C:/Users/lucas/AppData/Local/R/win-library/4.6")) {
  .libPaths(unique(c(.libPaths(), "C:/Users/lucas/AppData/Local/R/win-library/4.6")))
}


# Source all package functions from R/
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  tryCatch(source(f), error = function(e) NULL)
}

if (requireNamespace("ggplot2", quietly = TRUE) && exists("theme_scifigure")) {
  # Apply AI4SciFig theme as the default for all ggplots in the book
  ggplot2::theme_set(theme_scifigure(base_size = 7.5, grid = "both"))

  # Define fallback aliases so existing chapters inherit AI4SciFig styling automatically
  theme_minimal <- function(base_size = 7.5, ...) theme_scifigure(base_size = base_size, grid = "both", ...)
  theme_bw      <- function(base_size = 7.5, ...) theme_scifigure(base_size = base_size, grid = "both", ...)
  theme_classic <- function(base_size = 7.5, ...) theme_scifigure(base_size = base_size, grid = "none", ...)
}

# Color palettes
if (requireNamespace("ggplot2", quietly = TRUE) && exists("scale_color_scifigure")) {
  options(
    ggplot2.discrete.colour = function(...) scale_color_scifigure(palette = "nature", ...),
    ggplot2.discrete.fill   = function(...) scale_fill_scifigure(palette = "nature", ...)
  )
}
