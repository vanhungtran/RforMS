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

30 chapters in seven parts, plus summary, references, and four appendices.

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

### Part V — Statistical Modeling, Machine Learning, and Biological Interpretation (Chapters 19–24)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 19  | `20-experimental-design.qmd`        | Experimental Design, Replication, and Power           |
| 20  | `21-differential-abundance.qmd`     | Differential Abundance Analysis                       |
| 21  | `22-covariates-repeated.qmd`        | Covariates and Repeated Measures                      |
| 22  | `23-machine-learning.qmd`           | Machine Learning for MS Data                          |
| 23  | `24-biomarker-modeling.qmd`         | Biomarker Modeling                                    |
| 24  | `25-pathway-network.qmd`            | Pathway and Network Analysis                          |

### Part VI — Specialized and Emerging Applications (Chapters 25–28)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 25  | `26-single-cell-proteomics.qmd`     | Single-Cell Proteomics with `scp`                     |
| 26  | `27-ms-imaging.qmd`                 | Mass Spectrometry Imaging with `Cardinal`              |
| 27  | `28-biomarker-external-validation.qmd` | External Validation and Confounder-Robust Biomarker Panels |
| 28  | `29-isotope-tracing-networking.qmd` | Metabolite Discovery with Stable-Isotope Tracing and Molecular Networking |

### Part VII — Reproducible Reporting and Capstone (Chapters 29–30)

| Ch  | Source                              | Title                                                 |
|:---:|:------------------------------------|:------------------------------------------------------|
| 29  | `30-reproducible-reports.qmd`       | Build Reproducible MS Reports                         |
| 30  | `31-capstone-case-studies.qmd`      | Capstone: Two End-to-End Case Studies                 |

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

### Standalone Springer Nature Textbook Monograph PDFs

Each flagship chapter has been authored as an independent, publication-grade Springer Nature textbook monograph matching the rigorous pedagogy of Mary Regina Boland's *Health Data Analytics with R*. Each monograph includes formal mathematical definitions, analytical proofs, pedagogical callout boxes (*Methodological Pitfall*, *Biological Discovery*, *Takeaway Summary*, *Deep Dive*), full hands-on R laboratories with console outputs, a critical review of the *Five Deadly Sins*, multi-tier exercises, and exhaustive step-by-step worked solutions.

| Chapter | Monograph Title | Pages | Key Mathematical & Methodological Coverage | PDF Download |
|:---:|:---|:---:|:---|:---:|
| **08** | Chromatographic Feature Detection & Alignment with xcms | 15 | Continuous Wavelet Transform (Mexican Hat Ricker wavelet), ppm tolerance derivation, dynamic ROI tracking, Obiwarp dynamic time warping, Peak Density correspondence, gap-filling mechanics | [PDF](docs/latex/chapter_feature_detection_springer.pdf) |
| **14** | Quantitative Label-Free Proteomics (LFQ) | 15 | MaxLFQ Laplacian graph optimization, Top$N$ / iBAQ, Tukey's biweight loss, dynamic range compression | [PDF](docs/latex/chapter_lfq_proteomics_springer.pdf) |
| **15** | Multiplexed Isobaric Labeling (TMT / iTRAQ) | 11 | Co-isolation ratio compression proof, delta-method variance inflation, SPS-MS3/FAIMS filtering, IRS multi-batch scaling algebra | [PDF](docs/latex/chapter_tmt_labeling_springer.pdf) |
| **18** | Feature Normalization & Batch Correction | 18 | Rocke-Lorenzato two-component error model, Huber VSN $\operatorname{arsinh}$ derivation, Cyclic LOESS, PQN, ComBat empirical Bayes location-scale shrinkage, RUV | [PDF](docs/latex/chapter_normalization_batch_springer.pdf) |
| **19** | Missing Data in Mass Spectrometry | 17 | Rubin's MCAR/MAR/MNAR taxonomy, logistic diagnostic model, left-censored Gaussian Perseus/QRILC, $k$-NN, Rubin's pooling rules for multiple imputation, 4-strategy sensitivity audit | [PDF](docs/latex/chapter_missing_data_springer.pdf) |
| **21** | Differential Abundance Analysis | 19 | Smyth's empirical Bayes variance moderation, linear models, limma-trend vs. voom, Benjamini-Hochberg FDR, volcano geometry | [PDF](docs/latex/chapter_differential_abundance_springer.pdf) |
| **22** | Covariates, Longitudinal Profiling, and Repeated Measures | 17 | Linear Mixed-Effects Models (LMMs), Henderson's MME, REML derivation, Intraclass Correlation Coefficient (ICC), pseudoreplication variance inflation proof, `dream` empirical Bayes moderation, natural cubic splines | [PDF](docs/latex/chapter_covariates_repeated_springer.pdf) |
| **23** | Machine Learning in Mass Spectrometry | 18 | $p \gg n$ geometry, Cover's theorem on linear separability, Elastic Net grouping bound proof, Random Forest variance reduction, XGBoost, leakage taxonomy, Brier score | [PDF](docs/latex/chapter_machine_learning_springer.pdf) |
| **25** | Pathway Enrichment & Molecular Network Analysis | 17 | Hypergeometric ORA, background universe selection bias proof, GSEA weighted Kolmogorov-Smirnov running sums, NES, betweenness centrality, bipartite cross-omics networks | [PDF](docs/latex/chapter_pathway_network_springer.pdf) |
| **26** | Single-Cell Proteomics with `scp` and `QFeatures` | 15 | Poisson ion shot noise limit ($\text{CV} = 1/\sqrt{\lambda}$), SCoPE2 carrier ion suppression, multi-assay DAG hierarchy, pseudoreplication mixed modeling | [PDF](docs/latex/chapter_single_cell_springer.pdf) |

*Total: 10 flagship monographs, 162 publication-grade Springer textbook pages.*

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
- **scp** — Single-cell proteomics on top of `QFeatures`
- **Cardinal** — Mass spectrometry imaging

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
- **Affiliation:** CK-CARE Davos, Switzerland
- **Email:** <lucas.tran@ck-care.ch>
- **Repository:** <https://github.com/vanhungtran/RforMS>
