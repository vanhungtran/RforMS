# R for Mass Spectrometry

A comprehensive guide to analyzing mass spectrometry data using R and the R for Mass Spectrometry ecosystem.

## Book Structure

This book is organized into seven main parts:

### Part I: Foundations and Reproducible Setup
- **Preface** (`index.qmd`) - Welcome and overview
- **Chapter 1: Introduction to Mass Spectrometry Data Analysis with R** (`01-introduction.qmd`) - Workflows, concepts, and challenges
- **Chapter 2: The R and Bioconductor Ecosystem for Mass Spectrometry** (`02-r-ecosystem.qmd`) - Core packages and data structures
- **Chapter 3: Reproducible Project Setup for MS Workflows** (`03-reproducible-setup.qmd`) - Renv, targets, and project organization

### Part II: Data Import, Formats, and Quality Control
- **Chapter 4: MS File Formats and Data Conversion** (`04-file-formats.qmd`) - Vendor formats, open formats, and MSConvert
- **Chapter 5: Importing and Inspecting Raw MS Data in R** (`05-importing-data.qmd`) - Working with `Spectra` and `mzR`
- **Chapter 6: Experimental Metadata and Quality Control** (`06-metadata-qc.qmd`) - Sample annotation, batch effects, and QC reports

### Part III: Spectral Processing and Feature Quantification
- **Chapter 7: Spectral Processing** (`07-spectral-processing.qmd`) - Smoothing, baseline correction, and centroiding
- **Chapter 8: Chromatographic Peak Detection and Quantification** (`08-peak-detection.qmd`) - CentWave, EICs, and peak boundaries
- **Chapter 9: Retention Time Alignment and Feature Tables** (`09-rt-alignment.qmd`) - Correcting shifts and gap filling
- **Chapter 10: Spectral Similarity and Library Searching** (`10-spectral-library.qmd`) - Matching MS/MS spectra and building libraries

### Part IV: Visualization, Preprocessing, and Statistics
- **Chapter 11: Visualization of Mass Spectrometry Data** (`11-visualization.qmd`) - Creating publication-ready figures
- **Chapter 12: Preprocessing Intensity Matrices** (`12-preprocessing-matrices.qmd`) - Log transformation, scaling, and imputation
- **Chapter 13: Statistical Testing and Differential Abundance Analysis** (`13-statistical-testing.qmd`) - Limma, FDR, and batch correction
- **Chapter 14: Multivariate Analysis and Machine Learning** (`14-machine-learning.qmd`) - PCA, PLS-DA, Random Forests, and cross-validation

### Part V: Proteomics Workflows
- **Chapter 15: Shotgun Proteomics Data Analysis** (`15-shotgun-proteomics.qmd`) - PSMs, FDR, and protein inference
- **Chapter 16: Quantitative Proteomics** (`16-quantitative-proteomics.qmd`) - Label-free, SILAC, TMT, and iTRAQ
- **Chapter 17: Post-Translational Modification Analysis** (`17-ptm-analysis.qmd`) - Localization, motifs, and differential PTMs

### Part VI: Metabolomics Workflows
- **Chapter 18: Untargeted Metabolomics Data Processing** (`18-untargeted-metabolomics.qmd`) - Feature detection and QC with xcms
- **Chapter 19: Metabolite Annotation and Identification** (`19-metabolite-annotation.qmd`) - Accurate mass, adducts, and databases
- **Chapter 20: Targeted Metabolomics and Quantification** (`20-targeted-metabolomics.qmd`) - SRM, PRM, and calibration curves

### Part VII: Multi-Omics, Reporting, and Reproducible Delivery
- **Chapter 21: Multi-Omics Integration** (`21-multiomics-integration.qmd`) - Integrating proteomics and metabolomics
- **Chapter 22: Reproducible Reporting and Workflow Automation** (`22-reproducible-reporting.qmd`) - Quarto, Docker, and FAIR principles

### Summary and Appendices
- **Final Summary** (`summary.qmd`) - From Raw Data to Reproducible Insight
- **Appendix A: R Package Reference** (`appendix-a.qmd`)
- **Appendix B: MS File Format Reference** (`appendix-b.qmd`)
- **Appendix C: Common Adducts, Neutral Losses, and PTMs** (`appendix-c.qmd`)
- **Appendix D: Statistical and Reproducibility Quick Reference** (`appendix-d.qmd`)

Only the chapter files listed in `_quarto.yml` are included in the current book build.

## Quick Start

### Prerequisites

- R (>= 4.2.0)
- RStudio (recommended)
- Quarto (>= 1.3)

### Installation

```r
# Install core MS packages
install.packages("BiocManager")
BiocManager::install(c(
  "Spectra",
  "MsCoreUtils",
  "QFeatures",
  "xcms",
  "MetaboCoreUtils",
  "msdata",
  "MsDataHub"
))

# Install data manipulation and visualization
install.packages(c("tidyverse", "ggplot2", "plotly"))
```

### Building the Book

```bash
# In the project directory
quarto render
```

### Preview the Book

```bash
quarto preview
```

## Reading the Book

The book is available online at: https://vanhungtran.github.io/RforMS/

## Color Theme

The book uses a custom color palette featuring:
- **Headers**: Royal Blue 4 (`#27408B`)
- **Accents**: Magenta2 (`#EE00EE`) and Steel Blue (`#4682B4`)
- Clean, professional design optimized for readability

## Key R Packages Covered

- **Spectra** - Core MS data infrastructure
- **xcms** - LC-MS preprocessing and peak detection
- **QFeatures** - Quantitative proteomics workflows
- **MsCoreUtils** - MS data utilities
- **MetaboCoreUtils** - Metabolomics tools
- **PSMatch** - Peptide-spectrum matching

## Contributing

Contributions are welcome. Please feel free to:
- Report issues
- Suggest improvements
- Submit pull requests
- Request additional topics

## License

This work is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## Acknowledgments



## Contact

Author: Lucas VHH TRAN
Email: tranhungydhcm@gmail.com
Repository: https://github.com/vanhungtran/RforMS
