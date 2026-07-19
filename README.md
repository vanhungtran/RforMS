# Mass Spectrometry Data Analysis with R

**A Practical Guide to Reproducible Proteomics and Metabolomics with Bioconductor**

By Lucas VHH TRAN — <https://vanhungtran.github.io/RforMS/>

A comprehensive, hands-on guide to analyzing mass spectrometry data using R and the Bioconductor ecosystem. Built with Quarto.

## Book Structure

The book is organized into seven parts spanning 26 chapters, plus summary, references, and four appendices:

**Preface** (`index.qmd`)

### Part I — Foundations: Reproducible MS Analysis in R (Chapters 1–6)
| # | File | Title |
|---|------|-------|
| 1 | `02-what-is-ms.qmd` | What Is Mass Spectrometry? |
| 2 | `03-r-bioconductor.qmd` | R and Bioconductor for MS |
| 3 | `04-reproducible-project.qmd` | Build a Reproducible MS Analysis Project |
| 4 | `05-import-data.qmd` | Import MS Data |
| 5 | `06-data-objects.qmd` | Construct Analysis-Ready Data Objects |
| 6 | `07-initial-qc.qmd` | Perform Initial Quality Control |

### Part II — Feature Detection and Identification (Chapters 7–11)
| # | File | Title |
|---|------|-------|
| 7 | `08-detect-features.qmd` | Detect Chromatographic Features |
| 8 | `09-annotate-metabolites.qmd` | Annotate Metabolites |
| 9 | `10-library-search.qmd` | Match MS/MS Spectra to Libraries |
| 10 | `11-audit-psm.qmd` | Audit Peptide-Spectrum Matches |
| 11 | `12-protein-evidence.qmd` | Protein Evidence and Inference |

### Part III — Quantification Workflows (Chapters 12–16)
| # | File | Title |
|---|------|-------|
| 12 | `13-quant-proteomics.qmd` | Quantitative Proteomics Data Object |
| 13 | `14-lfq-proteomics.qmd` | Label-Free Quantification |
| 14 | `15-tmt-labeling.qmd` | TMT and Isobaric Labeling |
| 15 | `16-targeted-quant.qmd` | Targeted MS Quantification |
| 16 | `17-metabolomics-pipelines.qmd` | Metabolomics Quantification Pipelines |

### Part IV — Normalization, Batch Correction, and Missing Data (Chapters 17–18)
| # | File | Title |
|---|------|-------|
| 17 | `18-normalize-batches.qmd` | Normalize Across Samples and Batches |
| 18 | `19-missing-data.qmd` | Handle Missing Data |

### Part V — Statistical Modeling and Machine Learning (Chapters 19–23)
| # | File | Title |
|---|------|-------|
| 19 | `20-experimental-design.qmd` | Experimental Design, Replication, and Power |
| 20 | `21-differential-abundance.qmd` | Differential Abundance Analysis |
| 21 | `22-covariates-repeated.qmd` | Covariates and Repeated Measures |
| 22 | `23-machine-learning.qmd` | Machine Learning for MS Data |
| 23 | `24-biomarker-modeling.qmd` | Biomarker Modeling |

### Part VI — Biological Interpretation and Reproducible Reporting (Chapters 24–25)
| # | File | Title |
|---|------|-------|
| 24 | `25-pathway-network.qmd` | Pathway and Network Analysis |
| 25 | `26-reproducible-reports.qmd` | Build Reproducible MS Reports |

### Part VII — Capstone: End-to-End Case Studies (Chapter 26)
| # | File | Title |
|---|------|-------|
| 26 | `27-capstone-case-studies.qmd` | Capstone: Two End-to-End Case Studies |

### Back Matter
- `summary.qmd` — Summary and Future Directions
- `references.qmd` — References

### Appendices
- `appendix-a-packages.qmd` — Package Reference
- `appendix-b-formats.qmd` — MS File Formats
- `appendix-c-adducts.qmd` — Adduct Tables
- `appendix-d-statistics.qmd` — Statistics Reference

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

All runnable examples use bundled or versioned example data from Bioconductor packages (`msdata`, `faahKO`, `MsDataHub`) so they reproduce without external downloads. Data-access chapters that demonstrate importing from public repositories include stable accession numbers and retrieval instructions.

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

This work is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## Contact

Author: Lucas VHH TRAN
Email: tranhungydhcm@gmail.com
Repository: <https://github.com/vanhungtran/RforMS>
