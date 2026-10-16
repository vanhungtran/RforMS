#' AI4SciFig Full-Box ggplot2 Theme for Scientific Publication Figures
#'
#' Implements journal-compliant visual specifications (Nature, Science, Cell):
#' - Four visible black spines (closed box, panel.border 0.5 pt)
#' - Outward-pointing ticks (direction = "out", length = 2.5 pt)
#' - 7-8 pt typography scaled for actual printed dimensions (89 mm / 180 mm)
#' - Subtle dashed reference grid (#d9d9d9, 0.3 pt)
#' - Colorblind-safe palettes (Okabe-Ito and Nature Biomedical)
#'
#' @param base_size Numeric base font size in points (default: 7 pt)
#' @param base_family Font family name (default: "sans")
#' @param grid Character indicating grid style: "both", "y", "x", or "none"
#' @return A ggplot2 theme object
#' @export
theme_scifigure <- function(base_size = 7, base_family = "sans", grid = "both") {
  thm <- ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Four visible black spines forming a closed box
      panel.border       = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line          = ggplot2::element_blank(),
      # Outward ticks on bottom and left
      axis.ticks         = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length  = grid::unit(2.5, "pt"),
      # Integrated bold panel title
      plot.title         = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.subtitle      = ggplot2::element_text(size = base_size, colour = "#444444", margin = ggplot2::margin(b = 4)),
      axis.title         = ggplot2::element_text(size = base_size + 0.5, face = "plain", colour = "black"),
      axis.text          = ggplot2::element_text(size = base_size, colour = "black"),
      legend.title       = ggplot2::element_text(size = base_size, face = "bold"),
      legend.text        = ggplot2::element_text(size = base_size - 0.5),
      legend.background  = ggplot2::element_blank(),
      legend.key         = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "#f2f2f2", colour = "black", linewidth = 0.4),
      strip.text         = ggplot2::element_text(size = base_size, face = "bold", margin = ggplot2::margin(2, 3, 2, 3)),
      plot.margin        = ggplot2::margin(4, 5, 4, 4)
    )

  if (grid == "both") {
    thm <- thm + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "#d9d9d9", linewidth = 0.3, linetype = "dashed"),
      panel.grid.minor = ggplot2::element_blank()
    )
  } else if (grid == "y") {
    thm <- thm + ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(colour = "#d9d9d9", linewidth = 0.3, linetype = "dashed"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  } else if (grid == "x") {
    thm <- thm + ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = "#d9d9d9", linewidth = 0.3, linetype = "dashed"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  } else {
    thm <- thm + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }

  thm
}

#' Colorblind-safe palettes for scientific figures
#' @keywords internal
scifig_palettes <- list(
  okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"),
  nature    = c("#0F4D92", "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"),
  science   = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6"),
  nejm      = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97")
)

#' Discrete color scale for scientific figures
#'
#' @param palette Palette name: "nature", "okabe_ito", "science", or "nejm"
#' @param ... Arguments passed to ggplot2::discrete_scale
#' @export
scale_color_scifigure <- function(palette = "nature", ...) {
  cols <- scifig_palettes[[palette]] %||% scifig_palettes$nature
  ggplot2::discrete_scale("colour", "scifig", scales::manual_pal(values = cols), ...)
}

#' Discrete fill scale for scientific figures
#'
#' @param palette Palette name: "nature", "okabe_ito", "science", or "nejm"
#' @param ... Arguments passed to ggplot2::discrete_scale
#' @export
scale_fill_scifigure <- function(palette = "nature", ...) {
  cols <- scifig_palettes[[palette]] %||% scifig_palettes$nature
  ggplot2::discrete_scale("fill", "scifig", scales::manual_pal(values = cols), ...)
}

#' Save ggplot at physical publication dimensions in millimeters
#'
#' Ensures typography and strokes match exact journal print specifications.
#'
#' @param plot A ggplot object
#' @param filename File path without extension
#' @param journal Target journal ("nature", "science", "cell", "pnas")
#' @param width_class Column width class: "single" (89 mm) or "double" (180 mm)
#' @param height_mm Figure height in millimeters (default: 65 mm)
#' @param formats Vector of formats to export: c("png", "pdf")
#' @param dpi Resolution in dots per inch (default: 300)
#' @return Invisibly returns output paths
#' @export
scifig_save <- function(plot, filename, journal = "nature", width_class = "single",
                        height_mm = 65, formats = c("png", "pdf"), dpi = 300) {
  width_mm <- if (width_class == "double") 180 else 89

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  out_files <- character()

  for (fmt in formats) {
    target <- paste0(filename, ".", fmt)
    ggplot2::ggsave(
      filename = target,
      plot = plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = dpi
    )
    out_files <- c(out_files, target)
  }

  invisible(out_files)
}

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (length(a) == 1 && is.atomic(a) && is.na(a)) return(b)
  a
}
