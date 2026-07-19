# Book Data Directory

## Shipped Data (Available)

| Directory | Dataset | Chapters |
|-----------|---------|:------:|
| `pxd004886/` | PXD004886 DIA benchmark | 14, 19, 20 |
| `pxd010154/` | PXD010154 Tissue Atlas | 19, 23 |
| `mtbls38/` | MTBLS38 Metabolite Standards | 08 |
| `figures/` | Pre-computed PNG figures | 08, 14, 23 |

## R Package Data (Auto-Installed)

These are available after installing the required packages:
- `DEP::UbiLength` — Ch 14, 19, 20
- `DEP::UbiLength_ExpDesign` — Ch 20
- `faahKO` — Ch 17, 19
- `ropls::sacurine` — Ch 19
- `msdata::proteomics()` — Ch 06, 09
- `msdata::metabolomics()` — Ch 12

## Optional Data (Download Required)

These files are referenced in optional code blocks (`eval: false`) and are not required for rendering:

| File | Source | Chapter |
|------|--------|:------:|
| `pxd000001_tmt_psm_reporters.csv` | PRIDE PXD000001 via MsDataHub | 15 |
| `tmt_psm_table.csv` | Generate from TMT pipeline | 15 |
| `xcms_feature_matrix.rds` | Generate from xcms pipeline | 17 |
| `raw/sample.mzML` | msdata package or your own file | 04 |

To obtain these, run the companion download scripts in `../code/` and follow the chapter's preprocessing pipeline.
