#!/usr/bin/env Rscript
# Fast MTBLS38 analysis — EIC extraction + 1D peak detection on known masses
# ~2 seconds per compound instead of 3-5 minutes (100x faster than full CentWave)
suppressPackageStartupMessages({
  library(MSnbase)
  library(xcms)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

analysis_dir <- "analysis"
project_dir  <- "."
raw_dir       <- file.path(project_dir, "raw", "MTBLS38")
out_fig       <- file.path(analysis_dir, "results", "figures", "standards")
out_tbl       <- file.path(analysis_dir, "results", "tables", "standards")
dir.create(out_fig, recursive=TRUE, showWarnings=FALSE)
dir.create(out_tbl, recursive=TRUE, showWarnings=FALSE)

cat("=== MTBLS38 Fast Peak Analysis ===\n")

# Known [M+H]+ and [M-H]- masses for key compounds
known <- tibble::tribble(
  ~compound,                    ~mz_pos,    ~mz_neg,    ~rt_expected,
  "biotin",                     245.0955,   NA,         120,
  "cytidine",                   244.0928,   242.0782,   110,
  "proline",                    116.0706,   NA,          80,
  "quercetin",                  303.0499,   301.0354,   300,
  "kaempferol",                 287.0550,   285.0405,   310,
  "apigenin",                   271.0601,   NA,         320,
  "naringenin",                 273.0758,   NA,         300,
  "vanillin",                   NA,         151.0401,   260,
  "4-hydroxybenzoic_acid",      139.0390,   137.0244,   220,
  "nicotinamide",               123.0553,   NA,         100,
  "glycine_betaine",            118.0863,   NA,          80,
  "allantoin",                  159.0513,   NA,          90,
  "L-lysine",                   147.1128,   NA,          75,
  "ectoine",                    143.0815,   141.0669,    85,
  "cytosine",                   112.0505,   NA,          80,
  "guanine",                    152.0567,   NA,          90,
  "uridine",                    245.0768,   243.0623,   120,
  "jasmonic_acid",              NA,         209.1183,   350,
  "trans-zeatin",               220.1196,   NA,         200,
  "esculetin",                  179.0339,   177.0193,   280,
  "emodin",                     NA,         269.0455,   380,
  "1H-indole-1-acetic_acid",    NA,         174.0561,   250,
  "L-alanyl-L-alanine",         161.0921,   NA,          80,
  "(-)-epigallocatechin",       307.0812,   305.0667,   240,
  "D-arabinitol",               153.0758,   NA,          75,
  "citrulline",                 176.1029,   NA,          80,
  "cystathionine",              223.0747,   NA,          85,
  "adipic_acid",                NA,         145.0506,   180,
  "3-coumaric_acid",            NA,         163.0401,   270,
  "gentiobiose",                NA,         341.1089,    85
)

mzml_files <- list.files(raw_dir, pattern="\\.mzML$", full.names=TRUE)
cat(sprintf("Found %d mzML files, %d known compounds\n",
  length(mzml_files), nrow(known)))

# ── Process each known compound ─────────────────────────────────────────────
results <- list()
pb <- txtProgressBar(0, nrow(known), style=3)

for (i in seq_len(nrow(known))) {
  comp <- known$compound[i]
  mz_p <- known$mz_pos[i]
  mz_n <- known$mz_neg[i]

  # Find matching mzML file(s)
  pattern <- paste0(comp, "\\.mzML$")
  f <- grep(pattern, mzml_files, value=TRUE)[1]

  if (is.na(f)) {
    # Try negative mode
    f <- grep(paste0(comp, "_neg\\.mzML$"), mzml_files, value=TRUE)[1]
    if (is.na(f)) { results[[i]] <- NULL; setTxtProgressBar(pb,i); next }
    mz_target <- mz_n
    is_neg <- TRUE
  } else {
    mz_target <- mz_p
    is_neg <- FALSE
  }
  if (is.na(mz_target)) { results[[i]] <- NULL; setTxtProgressBar(pb,i); next }

  # Read raw data
  raw <- tryCatch(readMSData(f, mode="onDisk", msLevel.=1),
                  error=function(e) NULL)
  if (is.null(raw)) { results[[i]] <- NULL; setTxtProgressBar(pb,i); next }

  # Extract EIC around expected m/z (±25 ppm window)
  mz_win <- mz_target * c(1 - 25e-6, 1 + 25e-6)
  chr <- tryCatch(chromatogram(raw, mz=mz_win, aggregationFun="max"),
                  error=function(e) NULL)
  if (is.null(chr)) { results[[i]] <- NULL; setTxtProgressBar(pb,i); next }

  rtime_v <- rtime(chr[1,1])
  intens  <- intensity(chr[1,1])

  # Detect peaks in the 1D EIC (fast!)
  peaks <- tryCatch({
    findChromPeaks(chr, param=CentWaveParam(
      ppm=25, peakwidth=c(5, 60), snthresh=10,
      prefilter=c(3, 100), noise=100))
  }, error=function(e) NULL)

  n_peaks <- if (!is.null(peaks)) nrow(chromPeaks(peaks)) else 0

  # Find best peak (closest to expected m/z)
  best_peak <- NULL
  if (n_peaks > 0) {
    pk <- chromPeaks(peaks)
    mass_errors <- abs(pk[,"mz"] - mz_target) / mz_target * 1e6
    best_idx <- which.min(mass_errors)
    if (mass_errors[best_idx] < 100) {
      best_peak <- list(
        mz = pk[best_idx,"mz"],
        rt = pk[best_idx,"rt"],
        rtmin = pk[best_idx,"rtmin"],
        rtmax = pk[best_idx,"rtmax"],
        into = pk[best_idx,"into"],
        intb = pk[best_idx,"intb"],
        ppm_error = mass_errors[best_idx]
      )
    }
  }

  results[[i]] <- list(
    compound=comp, mz_target=mz_target, polarity=if(is_neg)"neg"else"pos",
    n_scans=length(rtime_v), tic_max=max(intens,na.rm=TRUE),
    n_peaks=n_peaks, best_peak=best_peak,
    rtime=rtime_v, intensity=intens, file=basename(f)
  )
  setTxtProgressBar(pb, i)
}
close(pb)

# ── Summary ─────────────────────────────────────────────────────────────────
cat("\n\n── Peak Detection Summary ──\n")
valid <- results[!vapply(results, is.null, logical(1))]
n_found <- sum(vapply(valid, function(r) !is.null(r$best_peak), logical(1)))
cat(sprintf("Compounds analyzed: %d / %d\n", length(valid), nrow(known)))
cat(sprintf("Peaks found at expected m/z: %d (%.0f%%)\n",
  n_found, n_found/length(valid)*100))

# Summary table
summary_rows <- lapply(valid, function(r) {
  bp <- r$best_peak
  tibble::tibble(
    compound=r$compound, mz_expected=r$mz_target, polarity=r$polarity,
    n_peaks_total=r$n_peaks,
    mz_found=if(!is.null(bp)) bp$mz else NA_real_,
    rt_found=if(!is.null(bp)) bp$rt else NA_real_,
    ppm_error=if(!is.null(bp)) bp$ppm_error else NA_real_,
    intensity=if(!is.null(bp)) bp$into else NA_real_
  )
})
summary_df <- dplyr::bind_rows(summary_rows)
write.csv(summary_df, file.path(out_tbl, "peak_validation.csv"), row.names=FALSE)
cat(sprintf("  ✓ Saved: peak_validation.csv (%d compounds)\n", nrow(summary_df)))

# ── Figures ─────────────────────────────────────────────────────────────────

# Figure 1: Gallery of 9 clean chromatograms
cat("\n── Generating figures ──\n")
top9 <- head(summary_df[order(summary_df$ppm_error), ], 9)

pdf(file.path(out_fig, "01_standards_gallery.pdf"), width=14, height=12)
plots <- list()
for (i in seq_len(min(9, nrow(top9)))) {
  comp_name <- top9$compound[i]
  r <- valid[[which(vapply(valid, function(x) x$compound==comp_name, logical(1)))[1]]]
  if (is.null(r)) next

  bp <- r$best_peak
  df <- data.frame(rt=r$rtime, int=r$intensity)

  p <- ggplot(df, aes(x=rt, y=int)) +
    geom_line(color="steelblue", linewidth=0.3) +
    labs(title=comp_name,
         subtitle=sprintf("m/z=%.4f found (%.1f ppm) | rt=%.1fs",
           bp$mz, bp$ppm_error, bp$rt),
         x="RT (s)", y="Intensity") +
    theme_minimal(base_size=9)

  # Zoom to peak region if peak found
  if (!is.null(bp)) {
    p <- p + coord_cartesian(xlim=c(max(0, bp$rtmin-30), bp$rtmax+30))
  }
  plots[[i]] <- p
}
print(wrap_plots(plots, ncol=3))
dev.off()
cat("  ✓ Saved: 01_standards_gallery.pdf\n")

# Figure 2: Peak anatomy — annotated peak for best match
best_compound <- summary_df[which.min(summary_df$ppm_error), ]
if (nrow(best_compound) > 0) {
  r <- valid[[which(vapply(valid, function(x)
    x$compound==best_compound$compound, logical(1)))[1]]]
  bp <- r$best_peak

  pdf(file.path(out_fig, "02_peak_anatomy.pdf"), width=10, height=5)
  df <- data.frame(rt=r$rtime, int=r$intensity)
  zoom <- df[df$rt >= bp$rtmin-20 & df$rt <= bp$rtmax+20, ]

  p <- ggplot(zoom, aes(x=rt, y=int)) +
    geom_area(fill="steelblue", alpha=0.15) +
    geom_line(color="steelblue", linewidth=0.6) +
    annotate("rect", xmin=bp$rtmin, xmax=bp$rtmax,
             ymin=0, ymax=max(zoom$int)*1.05,
             fill=NA, color="red", linetype="dashed", linewidth=0.8) +
    annotate("point", x=bp$rt, y=bp$into, color="red", size=3) +
    annotate("label", x=bp$rt, y=bp$into*1.1,
             label=sprintf("m/z = %.4f\nrt = %.1f s\nSNR = %.0f",
               bp$mz, bp$rt, bp$into/bp$intb),
             size=3.5, fill="white", alpha=0.8) +
    labs(title=sprintf("Peak Anatomy — %s", best_compound$compound),
         subtitle="Dashed box = CentWave peak boundaries | Red dot = apex",
         x="Retention Time (s)", y="Intensity") +
    theme_minimal(base_size=13)
  print(p)
  dev.off()
  cat("  ✓ Saved: 02_peak_anatomy.pdf\n")
}

# Figure 3: Mass accuracy distribution
pdf(file.path(out_fig, "03_mass_accuracy.pdf"), width=8, height=5)
acc <- summary_df[!is.na(summary_df$ppm_error), ]
p <- ggplot(acc, aes(x=ppm_error)) +
  geom_histogram(fill="steelblue", bins=20, alpha=0.8) +
  geom_vline(xintercept=c(-25, 25), linetype="dashed", color="red") +
  labs(title="Mass Accuracy of Detected Peaks vs Known Standards",
       subtitle=sprintf("%d compounds, median error = %.1f ppm",
         nrow(acc), median(acc$ppm_error, na.rm=TRUE)),
       x="Mass Error (ppm)", y="Count") +
  theme_minimal(base_size=12)
print(p)
dev.off()
cat("  ✓ Saved: 03_mass_accuracy.pdf\n")

cat(sprintf("\n=== Done: %s ===\n", Sys.time()))
