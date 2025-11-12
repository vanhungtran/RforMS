# R for Mass Spectrometry

A comprehensive guide to analyzing mass spectrometry data using R and the R for Mass Spectrometry ecosystem.

## 📚 Book Structure

This book is organized into five main parts:

### Part I: Foundations
- **Preface** (index.qmd) - Welcome and overview
- **Chapter 0: MS Principles** (intro.qmd) - Mass spectrometry theory and instrumentation
- **Chapter 1: Getting Started** (01-Introduction.qmd) - Hands-on introduction with R

### Part II: Core Techniques
- **Chapter 2: R Fundamentals** - Essential R programming for MS analysis
- **Chapter 3: Data Formats** - Importing and exporting MS data
- **Chapter 4: Spectral Preprocessing** - Cleaning and preparing spectra
- **Chapter 5: Peak Detection** - Finding and quantifying peaks

### Part III: Analysis & Visualization
- **Chapter 6: Data Visualization** - Creating informative plots
- **Chapter 7: Statistical Analysis** - Hypothesis testing and multivariate methods

### Part IV: Applications
- **Chapter 8: Metabolomics** - Metabolite identification and analysis
- **Chapter 9: Proteomics** - Protein identification workflows
- **Chapter 10: QFeatures** - Quantitative feature analysis

### Part V: Advanced Topics
- **Chapter 11: Advanced Methods** - Machine learning and method development

## 🚀 Quick Start

### Prerequisites

- R (≥ 4.2.0)
- RStudio (recommended)
- Quarto (≥ 1.3)

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

## 📖 Reading the Book

The book is available online at: https://vanhungtran.github.io/RforMS/

## 🎨 Color Theme

The book uses a custom color palette featuring:
- **Headers**: Royal Blue 4 (#27408B)
- **Accents**: Magenta2 (#EE00EE) and Steel Blue (#4682B4)
- Clean, professional design optimized for readability

## 📦 Key R Packages Covered

- **Spectra** - Core MS data infrastructure
- **xcms** - LC-MS preprocessing and peak detection
- **QFeatures** - Quantitative proteomics workflows
- **MsCoreUtils** - MS data utilities
- **MetaboCoreUtils** - Metabolomics tools
- **PSMatch** - Peptide-spectrum matching

## 🤝 Contributing

Contributions are welcome! Please feel free to:
- Report issues
- Suggest improvements
- Submit pull requests
- Request additional topics

## 📝 License

This work is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## 🙏 Acknowledgments

Built on the R for Mass Spectrometry ecosystem developed by:
- Laurent Gatto
- Johannes Rainer
- Sebastian Gibb
- And the entire RforMassSpectrometry community

## 📧 Contact

Author: Lucas VHH TRAN
Email: tranhungydhcm@gmail.com
Repository: https://github.com/vanhungtran/RforMS
