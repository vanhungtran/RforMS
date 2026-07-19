# Analysis Pipelines — Real Data for RforMS Book

Each pipeline processes a publicly available dataset from raw files to publication-ready figures. The output tables and figures are integrated into the book chapters listed below. All datasets are open-access from PRIDE, MetaboLights, or MassIVE.

## Pipeline → Chapter Map

| # | Pipeline | Dataset | Book Chapters | Key Output |
|---|----------|---------|:---:|-------------|
| 01 | `01_pxd004886_dia_analysis.R` | PXD004886 (Bruderer et al., 2017) | Ch 14, 19, 20 | 4,517 proteins, 6 DE hits, volcano plot, QC figures |
| 02 | `02_pxd010154_tissue_atlas.R` | PXD010154 (Human Proteome Atlas) | Ch 19, 22, 23 | 33,812 proteins, 12 organs, tissue-enriched markers |
| 03 | `03_mtbls38_peak_detection.R` | MTBLS38 (Metabolite Standards) | Ch 08, 11 | 26/26 peaks validated, 0.57 ppm median error |
| 04 | `04_pxd000547_spectronaut.R` | PXD000547 (Clinical DIA) | Ch 05 | Spectronaut binary import, paired DE |
| 00 | `00_annotation_utils.R` | Shared module | Ch 14 | UniProt → gene symbol annotation |

## Downloaders

| Script | Repository | Accession Types |
|--------|-----------|----------------|
| `download_geo_unified.R` | GEO, PRIDE, MetaboLights | GSE*, PXD*, MTBLS* |
| `download_mtbls.ps1` | MetaboLights (Windows) | MTBLS* |
| `download_pride.py` | PRIDE | PXD* |

## Usage

All pipelines are self-contained R scripts. Run from the project root:

```bash
# PXD004886 DIA analysis (requires raw/ directory with PRIDE data)
Rscript code/01_pxd004886_dia_analysis.R

# PXD010154 tissue atlas (requires 12 organ ZIP files in raw/PXD010154/)
Rscript code/02_pxd010154_tissue_atlas.R

# MTBLS38 peak detection (requires 71 mzML files in raw/MTBLS38/)
Rscript code/03_mtbls38_peak_detection.R
```

## Data Locations

| Dataset | Local Path | Cluster Path |
|---------|-----------|-------------|
| PXD004886 | `raw/PXD004886/` | `/cluster/data/ck_care/ltran/r4ms_book/raw/PXD004886/` |
| PXD010154 | `raw/PXD010154/` | `/cluster/data/ck_care/ltran/r4ms_book/raw/PXD010154/` |
| PXD000547 | `raw/PXD000547/` | `/cluster/data/ck_care/ltran/r4ms_book/raw/PXD000547/` |
| MTBLS38 | `raw/MTBLS38/` | Local Windows only (12 GB) |

## Reproducibility

All pipelines use a fixed random seed and record session information. The exact package versions are available in the pipeline log files. To reproduce:

1. Download raw data using the downloader scripts
2. Run the numbered pipeline
3. Compare output to `data/` directory shipped with the book
