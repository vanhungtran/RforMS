library(pdftools)
fig_dir <- "data/figures"
pdfs <- list.files(fig_dir, pattern="\\.pdf$", full.names=TRUE)
for (f in pdfs) {
  png_path <- sub("\\.pdf$", ".png", f)
  cat(sprintf("Converting: %s ... ", basename(f)))
  bitmap <- pdf_render_page(f, dpi = 150)
  png::writePNG(bitmap, png_path)
  cat(sprintf("OK (%d x %d)\n", dim(bitmap)[2], dim(bitmap)[1]))
}
cat("Done.\n")
