# Datasets Used in "Mass Spectrometry Data Analysis with R"

## Chapter → Dataset Map

| Ch | Chapter Title | Primary Dataset(s) | Source | Data Location |
|:--:|---------------|-------------------|--------|---------------|
| 04 | Reproducible Project | `sample.mzML` (example) | Generated | `data/raw/` |
| 05 | Import Data | PXD000547 `.sepr` files | PRIDE | `raw/PXD000547/` (416 MB)† |
| 06 | Data Objects | `msdata::proteomics()` | Bioc R package | Auto-installed |
| 07 | Initial QC | Simulated QC dataset | Generated in-chapter | — |
| **08** | **Detect Features** | **MTBLS38** (71 mzML, 26 standards) | MetaboLights | `data/mtbls38/` + `data/figures/`†† |
| 09 | Audit PSM | `msdata::proteomics()` | Bioc R package | Auto-installed |
| 12 | Library Search | `msdata::metabolomics()` | Bioc R package | Auto-installed |
| 13 | Quant Proteomics | CPTAC (via MsDataHub) | Bioc R package | Auto-installed |
| **14** | **LFQ Proteomics** | PXD004886 (4,517 proteins) + `DEP::UbiLength` | PRIDE + Bioc | `data/pxd004886/` |
| 15 | TMT Labeling | PXD000001 | PRIDE / MsDataHub | Not shipped (optional) |
| 17 | Metabolomics | `faahKO` + MTBLS234 | Bioc + MetaboLights | Auto-installed |
| **19** | **Missing Data** | PXD004886 + PXD010154 + `DEP::UbiLength` | PRIDE + Bioc | `data/pxd004886/` + `data/pxd010154/` |
| **20** | **Differential Abundance** | PXD004886 (6 DE hits) + `DEP::UbiLength` | PRIDE + Bioc | `data/pxd004886/` |
| **23** | **Pathway-Network** | PXD010154 (12 organs) + simulated | PRIDE | `data/pxd010154/` + `data/figures/` |

† Large raw files — download separately with `code/download_geo_unified.R` or `code/download_mtbls.ps1`
†† 71 mzML files = ~12 GB — download separately, analysis results in `data/mtbls38/`

**Bold** = Chapters modified with real-data `.real-data` sections.

## Shipped Data Files

### `data/pxd004886/` — DIA Benchmark (PXD004886)
| File | Rows | Description |
|------|:----:|-------------|
| `DE_results.csv` | 4,517 | limma differential expression (all proteins) |
| `DE_results_annotated.csv` | 4,517 | DE results with gene symbols + significance flags |

### `data/pxd010154/` — Tissue Atlas (PXD010154)
| File | Rows | Description |
|------|:----:|-------------|
| `organ_summary.csv` | 14 | Proteins per organ, mean/median expression |
| `tissue_enriched_proteins_top5.csv` | 60 | Top 5 enriched proteins per organ (log2FC) |

### `data/mtbls38/` — Metabolite Standards (MTBLS38)
| File | Rows | Description |
|------|:----:|-------------|
| `peak_validation.csv` | 26 | CentWave-detected peaks vs. known masses |

### `data/figures/` — Pre-computed Figures
| File | From | Chapter |
|------|------|:------:|
| `01_standards_gallery.png` | MTBLS38 analysis | 08 |
| `02_peak_anatomy.png` | MTBLS38 analysis | 08 |
| `03_mass_accuracy.png` | MTBLS38 analysis | 08 |
| `01_density_after_norm.png` | PXD004886 analysis | 14 |
| `03_PCA.png` | PXD004886 analysis | 14 |
| `05_volcano.png` | PXD004886 analysis | 14 |
| `01_tissue_correlation_heatmap.png` | PXD010154 analysis | 23 |
| `02_tissue_pca.png` | PXD010154 analysis | 23 |
| `03_tissue_enriched_heatmap.png` | PXD010154 analysis | 23 |

## R Package Datasets (Auto-Available)

These datasets are bundled with their respective packages and installed automatically:

```r
BiocManager::install(c("DEP", "msdata", "faahKO", "MsDataHub", "ropls"))
```

| Dataset | Package | Load With | Chapters |
|---------|---------|-----------|:------:|
| UbiLength | DEP | `data("UbiLength")` | 14, 19, 20 |
| UbiLength_ExpDesign | DEP | `data("UbiLength_ExpDesign")` | 20 |
| faahKO | faahKO | `library(faahKO)` | 17, 19 |
| sacurine | ropls | `data("sacurine")` | 19 |
| proteomics (mzML) | msdata | `msdata::proteomics()` | 06, 09 |
| metabolomics (mzML) | msdata | `msdata::metabolomics()` | 12 |

## Downloading Large Raw Data

The full raw data files are not shipped with the book due to size (108+ GB total). Use the companion download scripts:

```bash
# Windows PowerShell — MetaboLights mzML files
./code/download_mtbls.ps1 -Datasets MTBLS38,MTBLS1455,MTBLS234

# R — Universal downloader (GEO, PRIDE, MetaboLights)
Rscript code/download_geo_unified.R
# Then: download_geo_unified(c("PXD004886","PXD010154","MTBLS404"), download=TRUE)
```

## Reproducing the Analyses

1. Download raw data using the scripts above
2. Run the numbered pipelines in `code/`:
   - `code/01_pxd004886_dia_analysis.R` → Ch 14, 19, 20
   - `code/02_pxd010154_tissue_atlas.R` → Ch 19, 23
   - `code/03_mtbls38_peak_detection.R` → Ch 08
3. Compare output to shipped `data/` files

All pipelines record session information and use fixed random seeds for reproducibility.
