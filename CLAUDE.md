# RforMS — Project Context for Claude

## Book Identity

**Title:** Mass Spectrometry Data Analysis with R
**Subtitle:** A Practical Guide to Reproducible Proteomics and Metabolomics with Bioconductor
**Author:** Lucas VHH TRAN
**Published at:** https://vanhungtran.github.io/RforMS/
**Repo:** https://github.com/vanhungtran/RforMS

Built with **Quarto book** (`_quarto.yml`). Source `.qmd` files are in the project root; rendered HTML goes to `docs/`. All source files (`.qmd`, `.R`, `.bib`, `.csv`, `_quarto.yml`, `CLAUDE.md`) are tracked in git; generated artifacts (`docs/`, `*_cache/`, `*_files/`, `r_libs/`, `.quarto/`) are git-ignored.

---

## Chapter Map

> **Rendered chapter number** (col 1) = the number Quarto outputs and that all textual "Chapter N" cross-references use. File prefixes are 02–27; rendered chapters are 1–26.

| # | File | Title |
|---|------|-------|
| 1 | `02-what-is-ms.qmd` | What Is Mass Spectrometry? |
| **Part I — Foundations (Ch 2–6)** |||
| 2 | `03-r-bioconductor.qmd` | R and Bioconductor for MS |
| 3 | `04-reproducible-project.qmd` | Build a Reproducible MS Analysis Project |
| 4 | `05-import-data.qmd` | Import MS Data |
| 5 | `06-data-objects.qmd` | Construct Analysis-Ready Data Objects |
| 6 | `07-initial-qc.qmd` | Perform Initial Quality Control |
| **Part II — Feature Detection & Identification (Ch 7–11)** — *regrouped: metabolomics ID before proteomics ID* |||
| 7 | `08-detect-features.qmd` | Detect Chromatographic Features |
| 8 | `09-annotate-metabolites.qmd` | Annotate Metabolites |
| 9 | `10-library-search.qmd` | Match MS/MS Spectra to Libraries |
| 10 | `11-audit-psm.qmd` | Audit Peptide-Spectrum Matches |
| 11 | `12-protein-evidence.qmd` | Protein Evidence and Inference |
| **Part III — Quantification (Ch 12–16)** |||
| 12 | `13-quant-proteomics.qmd` | Quantitative Proteomics Data Object |
| 13 | `14-lfq-proteomics.qmd` | Label-Free Quantification |
| 14 | `15-tmt-labeling.qmd` | TMT and Isobaric Labeling |
| 15 | `16-targeted-quant.qmd` | Targeted MS Quantification |
| 16 | `17-metabolomics-pipelines.qmd` | Metabolomics Quantification Pipelines |
| **Part IV — Normalization, Batch Correction & Missing Data (Ch 17–18)** |||
| 17 | `18-normalize-batches.qmd` | Normalize Across Samples and Batches |
| 18 | `19-missing-data.qmd` | Handle Missing Data |
| **Part V — Statistical Modeling & Machine Learning (Ch 19–23)** |||
| 19 | `20-experimental-design.qmd` | Experimental Design, Replication, and Power *(new)* |
| 20 | `21-differential-abundance.qmd` | Differential Abundance Analysis |
| 21 | `22-covariates-repeated.qmd` | Covariates and Repeated Measures |
| 22 | `23-machine-learning.qmd` | Machine Learning for MS Data *(new)* |
| 23 | `24-biomarker-modeling.qmd` | Biomarker Modeling |
| **Part VI — Interpretation & Reporting (Ch 24–25)** |||
| 24 | `25-pathway-network.qmd` | Pathway and Network Analysis |
| 25 | `26-reproducible-reports.qmd` | Build Reproducible MS Reports |
| **Part VII — Capstone (Ch 26)** |||
| 26 | `27-capstone-case-studies.qmd` | Capstone: Two End-to-End Case Studies *(new)* |
| — | `summary.qmd` | Summary and Future Directions |
| — | `references.qmd` | References |
| **Appendices** |||
| A | `appendix-a-packages.qmd` | Package Reference |
| B | `appendix-b-formats.qmd` | MS File Formats |
| C | `appendix-c-adducts.qmd` | Adduct Tables |
| D | `appendix-d-statistics.qmd` | Statistics Reference |

> **Note:** `MS_basic.qmd` was archived to `MS_basic.qmd.archived` (orphaned infographic, not in the book).

---

## Chapter Descriptions

**Ch 1 — What Is Mass Spectrometry?**
Explains what a mass spectrometer measures (*m/z*, not mass), why the same molecule can produce multiple peaks (charge states, adducts, in-source fragments), and how acquisition strategy — full-scan, DDA, DIA, SRM/MRM — shapes missingness patterns before any modeling. Covers ion suppression and contamination as the two most common data artifacts in LC–ESI–MS.

**Ch 2 — R and Bioconductor for MS**
Practical orientation to R and Bioconductor for computational MS. Core data structures (vectors, data frames, lists, functions), the four key Bioconductor containers (`Spectra`, `MsExperiment`, `SummarizedExperiment`, `QFeatures`), and a complete end-to-end pipeline from raw mzML to a differential result.

**Ch 3 — Build a Reproducible MS Analysis Project**
Infrastructure for reproducible MS analysis. Project directory layout, Quarto literate programming, `renv` for package version locking, `targets` pipeline caching with dynamic branching, parameter provenance, `testthat` unit testing, and GitHub Actions CI.

**Ch 4 — Import MS Data**
History and structure of MS file formats: proprietary vendor formats (Thermo `.RAW`, Sciex `.WIFF`, etc.) vs. open standard `mzML` (HUPO-PSI, PSI-MS CV). Converting with `msconvert`. Loading mzML into R with `Spectra` + `MsBackendMzR`; choosing backends for different data sizes.

**Ch 5 — Construct Analysis-Ready Data Objects**
In-depth coverage of the four Bioconductor containers with real data: constructing and filtering `Spectra`, linking files to sample metadata in `MsExperiment`, building `SummarizedExperiment` from a feature matrix, assembling a `QFeatures` hierarchy. When to use each and how to convert between them.

**Ch 6 — Perform Initial Quality Control**
Pre-modeling QC: annotating samples (group, batch, blank, QC pool), TIC and BPC visualisation across all files, retention-time stability assessment, mass accuracy on known QC peaks, blank contamination and carryover detection, and generating a QC report.

**Ch 7 — Detect Chromatographic Features**
xcms feature detection pipeline: CentWave peak detection (`findChromPeaks`, parameters: ppm, peakwidth, snthresh), Obiwarp retention-time alignment (`adjustRtime`), feature correspondence (`groupChromPeaks`, `PeakDensityParam`), peak filling (`fillChromPeaks`), and feature matrix export as `SummarizedExperiment`. Quality metrics: detection rate, CV in QC pools.

**Ch 8 — Annotate Metabolites from Feature Tables**
MSI confidence levels (1–4) as the annotation framework. Adduct deconvolution with CAMERA, exact-mass database search with `MetaboCoreUtils`, isotope pattern scoring. Computational path to Level 2 (putative, library-matched) and Level 3 (putative class) annotations.

**Ch 9 — Match MS/MS Spectra to Libraries**
Cosine similarity spectral library search with `MetaboAnnotation`. Building local reference libraries with `CompoundDb`, querying MassBank/MoNA/HMDB, handling adducts in matching, interpreting scores with confidence thresholds.

**Ch 10 — Audit Peptide-Spectrum Matches**
PSM confidence filtering: loading PSM tables, target–decoy FDR with `PSMatch`, diagnosing shared peptides, missed cleavages, charge anomalies and contaminant hits, building a clean identification report for downstream quantification.

**Ch 11 — Protein Evidence and Inference**
The protein inference problem in bottom-up proteomics. Mapping peptides to proteins with `dplyr`, distinguishing unique/shared/razor peptides, `igraph` network analysis for protein groups, peptide-to-protein aggregation with `QFeatures`, and filtering to avoid false protein calls.

**Ch 12 — Build and Audit a Quantitative Proteomics Data Object**
`QFeatures` as the organisational backbone for quantitative proteomics. PSM → peptide → protein hierarchy, importing search-engine output, `filterFeatures()` and `filterNA()`, and how aggregation propagates uncertainty.

**Ch 13 — Label-Free Quantification**
Complete LFQ workflow with `DEP` + `QFeatures`: contaminant/reverse-hit filtering, normalization (median, quantile, cyclic-loess), median and robust peptide-to-protein summarization, IRS for multi-batch data, PCA diagnostics.

**Ch 14 — TMT and Isobaric Labeling Analysis**
Reporter-ion quantification (TMT-6/10/16, iTRAQ): loading reporter columns, sample-loading normalization, within-channel median scaling, purity correction, co-isolation interference and ratio compression.

**Ch 15 — Targeted MS Quantification**
SRM/MRM targeted quantification: loading transition data, peak integration, calibration curves, LOD/LOQ estimation, Skyline import into R.

**Ch 16 — Metabolomics Quantification Pipelines**
End-to-end untargeted metabolomics with `xcms` + `MsExperiment`: raw mzML → normalized feature matrix, bridging Ch 7 feature detection to the normalization and statistics chapters.

**Ch 17 — Normalize Across Samples and Batches**
Post-quantification matrix normalization: log-transformation, median, quantile, LOESS, TIC-based correction, ComBat batch correction. CV in QC pools and PCA to compare methods.

**Ch 18 — Handle Missing Data**
MCAR/MAR/MNAR missingness diagnosis. Imputation strategies: minimum-value, KNN, BPCA. Sensitivity analysis across imputation methods and downstream consequences on differential results.

**Ch 19 — Experimental Design, Replication, and Power**
Designing MS experiments: randomization, blocking, biological vs. technical replication, sources of variation, power analysis with `msmsTests` and simulation, sample size planning for proteomics and metabolomics studies.

**Ch 20 — Differential Abundance Analysis**
`limma` for MS: building design matrices (two-group, factorial, blocking), empirical Bayes moderation, Benjamini-Hochberg FDR, volcano plots, ranked feature lists.

**Ch 21 — Covariates and Repeated Measures**
Mixed-effects models for paired/longitudinal MS data, blocking on subject, continuous covariates, `variancePartition` for variance decomposition, paired vs. unpaired `limma` comparison.

**Ch 22 — Machine Learning for MS Data**
Supervised and unsupervised ML with `tidymodels`: regularized regression, tree ensembles, class imbalance handling, model calibration, PCA and PLS-DA diagnostics. Nested cross-validation with feature selection inside each fold.

**Ch 23 — Biomarker Modeling Without Data Leakage**
Three forms of leakage (preprocessing, feature selection, hyperparameter tuning). Nested cross-validation, random forest and regularized regression, ROC with bootstrap CIs, survival analysis, reporting to publication standards.

**Ch 24 — Pathway and Network Analysis**
Multi-omics integration in `MultiAssayExperiment`, cross-omics correlation matrices, DIABLO multi-block PLS-DA (`mixOmics`), KEGG pathway enrichment (`clusterProfiler`), circos and correlation network visualisation.

**Ch 25 — Build Reproducible MS Reports**
Parameterized Quarto reports, linking reports to `targets` pipelines, publication-ready figure export, depositing to PRIDE (proteomics) and MetaboLights (metabolomics).

---

## Key Design Decisions Made in This Session

- **Chapter 03 was split** into Ch 2 (R + Bioconductor) and Ch 3 (Reproducible Project) — original `03-first-project.qmd` was too dense.
- **Part IV** renamed from "Data Cleaning and Preprocessing" to "Matrix Normalization and Batch Correction" — the old name implied raw spectral work; this part operates on the quantified feature matrix.
- **`07-detect-features.qmd`** was rewritten — the file previously contained a duplicate of the Ch 1 MS introduction content.
- **`summary.qmd`** was fully rewritten — it described a 7-part, 17-chapter structure from an old version of the book.
- **`index.qmd` preface** was aligned to the actual book structure and title.
- **Chapter numbering:** Quarto renders `02-what-is-ms.qmd` as Chapter 1; `03-r-bioconductor.qmd` as Chapter 2; etc. The files are numbered 02–27 but Quarto outputs chapters 1–26.

---

## Coding Conventions

- All examples use `msdata`, `faahKO`, and `MsDataHub` packages for reproducibility. Bundled or versioned example data are preferred; external datasets require a stable accession, retrieval instructions, and a non-network fallback where feasible.
- Code chunks use `#| eval: false` for setup/installation blocks; `#| eval: true` for demonstration code.
- `echo: true` globally in `_quarto.yml`; `code-fold: true` for all chapters.
- Cross-references use "Chapter N" where N is the Quarto-rendered chapter number (not the file prefix).
- Citation style: `[@bibkey]` referencing `references.bib`; no bare URLs as citation text.
- **Chapter structure:** every computational chapter contains exactly one Learning Objectives, one Summary, one Exercises, and one Session Information section. The conceptual opening chapter (`02-what-is-ms.qmd`, rendered as Chapter 1) is exempt from Learning Objectives, Summary, and Session Information — it is a non-computational introduction. `summary.qmd` and `references.qmd` are also exempt as back-matter.
