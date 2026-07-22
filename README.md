# Mass Spectrometry Data Analysis with R

[![Website](https://img.shields.io/badge/website-RforMS-steelblue?style=flat-square&logo=github)](https://vanhungtran.github.io/RforMS/)
[![Built with Quarto](https://img.shields.io/badge/built%20with-Quarto-397eba?style=flat-square&logo=quarto)](https://quarto.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.4-276dc2?style=flat-square&logo=r)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-%E2%89%A5%203.20-3792b5?style=flat-square&logo=r)](https://bioconductor.org/)
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/license-CC%20BY--NC--SA%204.0-lightgrey?style=flat-square)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## A Practical Guide to Reproducible Proteomics and Metabolomics with Bioconductor

By **Lucas VHH TRAN** — published at <https://vanhungtran.github.io/RforMS/>

A comprehensive, hands-on guide to analyzing mass spectrometry data with R and the
[Bioconductor](https://bioconductor.org/) ecosystem. Covers the full workflow — raw
data import, feature detection, metabolite annotation, peptide identification, protein
inference, quantification (label-free, TMT, targeted), normalization, batch correction,
missing-data handling, differential abundance, machine learning, biomarker modeling,
pathway analysis, and reproducible reporting. Built with [Quarto](https://quarto.org/).

---

## Book Structure

26 chapters in seven parts, plus summary, references, and six appendices.

Chapter 1 introduces what a mass spectrometer measures. The remaining chapters are
grouped into seven parts:

### Part I — Foundations: Reproducible MS Analysis in R (Chapters 2–6)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
|  2  | `03-r-bioconductor.qmd`             | R and Bioconductor for MS                             |
|  3  | `04-reproducible-project.qmd`       | Build a Reproducible MS Analysis Project              |
|  4  | `05-import-data.qmd`                | Import MS Data                                        |
|  5  | `06-data-objects.qmd`               | Construct Analysis-Ready Data Objects                 |
|  6  | `07-initial-qc.qmd`                 | Perform Initial Quality Control                       |

### Part II — Feature Detection and Identification (Chapters 7–11)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
|  7  | `08-detect-features.qmd`            | Detect Chromatographic Features                       |
|  8  | `09-annotate-metabolites.qmd`       | Annotate Metabolites                                  |
|  9  | `10-library-search.qmd`             | Match MS/MS Spectra to Libraries                      |
| 10  | `11-audit-psm.qmd`                  | Audit Peptide-Spectrum Matches                        |
| 11  | `12-protein-evidence.qmd`           | Protein Evidence and Inference                        |

### Part III — Quantification Workflows (Chapters 12–16)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 12  | `13-quant-proteomics.qmd`           | Quantitative Proteomics Data Object                   |
| 13  | `14-lfq-proteomics.qmd`             | Label-Free Quantification                             |
| 14  | `15-tmt-labeling.qmd`               | TMT and Isobaric Labeling                             |
| 15  | `16-targeted-quant.qmd`             | Targeted MS Quantification                            |
| 16  | `17-metabolomics-pipelines.qmd`     | Metabolomics Quantification Pipelines                 |

### Part IV — Normalization, Batch Correction, and Missing Data (Chapters 17–18)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 17  | `18-normalize-batches.qmd`          | Normalize Across Samples and Batches                  |
| 18  | `19-missing-data.qmd`               | Handle Missing Data                                   |

### Part V — Statistical Modeling and Machine Learning (Chapters 19–23)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 19  | `20-experimental-design.qmd`        | Experimental Design, Replication, and Power           |
| 20  | `21-differential-abundance.qmd`     | Differential Abundance Analysis                       |
| 21  | `22-covariates-repeated.qmd`        | Covariates and Repeated Measures                      |
| 22  | `23-machine-learning.qmd`           | Machine Learning for MS Data                          |
| 23  | `24-biomarker-modeling.qmd`         | Biomarker Modeling                                    |

### Part VI — Biological Interpretation and Reproducible Reporting (Chapters 24–25)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 24  | `25-pathway-network.qmd`            | Pathway and Network Analysis                          |
| 25  | `26-reproducible-reports.qmd`       | Build Reproducible MS Reports                         |

### Part VII — Capstone (Chapter 26)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 26  | `27-capstone-case-studies.qmd`      | Capstone: Two End-to-End Case Studies                 |

### Back Matter

| File                                     | Content                               |
|:-----------------------------------------|:--------------------------------------|
| `index.qmd`                              | Preface                               |
| `summary.qmd`                            | Summary and Future Directions         |
| `references.qmd`                         | References                            |

### Appendices

| File                                                | Content                              |
|:----------------------------------------------------|:-------------------------------------|
| `appendix-a-packages.qmd`                           | Package Reference                    |
| `appendix-b-formats.qmd`                            | MS File Formats                      |
| `appendix-c-adducts.qmd`                            | Adduct Tables                        |
| `appendix-d-statistics.qmd`                         | Statistics Reference                 |
| `appendix-e-single-cell-proteomics.qmd`             | Single-Cell Proteomics               |
| `appendix-f-ms-imaging.qmd`                         | Mass Spectrometry Imaging            |

## Quick Start

### Prerequisites

- R (>= 4.4)
- Bioconductor (>= 3.20)
- Quarto (>= 1.3)

### Building the Book

```bash
# Clone the repository
git clone https://github.com/vanhungtran/RforMS.git
cd RforMS

# Install required packages (see Appendix A for the full list)
# Then render:
quarto render
```

## Reading the Book

The book is available online at: <https://vanhungtran.github.io/RforMS/>

## Reproducibility

All runnable examples use bundled or versioned example data from Bioconductor
packages (`msdata`, `faahKO`, `MsDataHub`) so they reproduce without external
downloads. Data-access chapters that demonstrate importing from public repositories
include stable accession numbers and retrieval instructions.

## Key R Packages Covered

- **Spectra** — Core MS data infrastructure
- **xcms** — LC-MS preprocessing and peak detection
- **QFeatures** — Quantitative proteomics data structures
- **MsCoreUtils** — MS data utilities
- **MetaboCoreUtils** — Metabolomics tools
- **PSMatch** — Peptide-spectrum matching
- **limma** — Differential abundance analysis
- **DEP** — Label-free proteomics workflows

## Contributing

Contributions are welcome. Please feel free to:

- Report issues
- Suggest improvements
- Submit pull requests

## License

This work is licensed under
[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).

## Contact

- **Author:** Lucas VHH TRAN
- **Email:** <tranhungydhcm@gmail.com>
- **Repository:** <https://github.com/vanhungtran/RforMS>
