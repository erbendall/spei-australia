# Gridded SPEI for Australia 1980–2024 at 0.01° Resolution

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![License: MIT](https://img.shields.io/badge/Code_License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19447029.svg)](https://doi.org/10.5281/zenodo.19447029)

Gridded Standardised Precipitation Evapotranspiration Index (SPEI) for the Australian continent at **0.01° (~1 km) resolution**, computed at 3, 6, 12 and 24 month timescales for the period January 2010 – December 2024.

**Author:** Eli R. Bendall  
Hawkesbury Institute for the Environment, Western Sydney University, NSW, Australia  
Funding: National Environmental Science Program (NESP) Resilient Landscapes Hub  
ORCID: [0000-0002-9120-3457](https://orcid.org/0000-0002-9120-3457)  
Code developed with AI assistance from Claude (Anthropic) — https://claude.ai

---

## Dataset Overview

| Attribute | Value |
|---|---|
| Spatial resolution | 0.01° (~1 km) |
| Spatial extent | 112–154°E, 9–44°S (Australian continent) |
| Temporal coverage (output) | January 2010 – December 2024 |
| Temporal coverage (NetCDF) | January 1980 – December 2024 |
| SPEI timescales | 3, 6, 12, 24 months |
| Reference period | January 1980 – December 2009 (30 years) |
| PET method | Hargreaves (Hargreaves & Samani 1985) |
| Distribution | Log-logistic |
| Input data | ANUClimate v2.0 (Hutchinson et al. 2021) |
| CRS | GDA94 (EPSG:4283) |
| Missing value (NetCDF) | -9999 |

---

## File Structure

```
spei-australia/
├── README.md
├── LICENSE
├── CITATION.cff
│
├── SPEI3_monthly_tifs/        # SPEI-3 monthly GeoTIFFs (Jan 2010 – Dec 2024)
│   ├── SPEI3_201001.tif
│   ├── SPEI3_201002.tif
│   └── ...  (180 files)
│
├── SPEI6_monthly_tifs/        # SPEI-6 monthly GeoTIFFs (180 files)
├── SPEI12_monthly_tifs/       # SPEI-12 monthly GeoTIFFs (180 files)
├── SPEI24_monthly_tifs/       # SPEI-24 monthly GeoTIFFs (180 files)
│
├── SPEI3_1980_2024.nc         # SPEI-3 full timeseries NetCDF (540 months)
├── SPEI6_1980_2024.nc         # SPEI-6 full timeseries NetCDF (540 months)
├── SPEI12_1980_2024.nc        # SPEI-12 full timeseries NetCDF (540 months)
└── SPEI24_1980_2024.nc        # SPEI-24 full timeseries NetCDF (540 months)
```

**Monthly GeoTIFFs** cover January 2010 – December 2024 (180 files per scale).  
**NetCDF files** include the full 1980–2024 timeseries (540 months per scale), including the reference period used for SPEI calibration.

File naming convention for GeoTIFFs: `SPEI{scale}_{YYYYMM}.tif`  
(e.g. `SPEI12_202001.tif` = SPEI-12 for January 2020)

---

## Input Data: ANUClimate v2.0

Input climate data were sourced from **ANUClimate v2.0** (Hutchinson et al. 2021), a high-resolution gridded climate dataset for the Australian continent produced by the Australian National University (ANU) and hosted at the National Computational Infrastructure (NCI).

ANUClimate v2.0 is generated using thin-plate spline interpolation of observational data from the Australian Bureau of Meteorology (BoM) station network, with elevation as a covariate. The dataset provides monthly climate surfaces at 0.01° (~1 km) resolution from 1900 to present, making it one of the highest-resolution national climate datasets available for Australia.

| Variable | Description | Units |
|---|---|---|
| `rain` | Monthly total rainfall | mm month⁻¹ |
| `tmax` | Monthly mean daily maximum temperature | °C |
| `tmin` | Monthly mean daily minimum temperature | °C |
| `srad` | Monthly mean daily solar radiation | MJ m⁻² day⁻¹ |

ANUClimate v2.0 is publicly accessible via the NCI THREDDS server (no NCI account required):  
https://thredds.nci.org.au/thredds/catalog/gh70/ANUClimate/v2-0/stable/month/

**Citation:**  
Hutchinson, M., Xu, T., Kesteven, J., Marang, I., and Evans, B. (2021). ANUClimate v2.0 monthly climate data for the Australian continent, 1900 to present. NCI Australia. https://doi.org/10.25914/60a10acd183a2

---

## Methods

### 1. Potential Evapotranspiration (PET)

PET was estimated using the Hargreaves method (Hargreaves & Samani 1985):

```
PET = 0.0023 × (srad × days_in_month) × (tmean + 17.8) × (tmax − tmin)^0.5
```

where `srad` is in MJ m⁻² day⁻¹, `tmean = (tmax + tmin) / 2` in °C, and `days_in_month` converts daily radiation to a monthly total. The Hargreaves method was selected as it requires only temperature and radiation inputs, both of which are available in ANUClimate v2.0 at the required resolution and temporal coverage.

### 2. Climatic Water Balance (CWB)

Monthly CWB was computed as:

```
CWB = rainfall − PET
```

### 3. SPEI Computation

SPEI was computed from the CWB timeseries using the `SPEI` R package (Beguería & Vicente-Serrano 2023) with the following settings:

- **Distribution:** log-logistic, fitted using unbiased probability weighted moments (ub-pwm)
- **Reference period:** January 1980 – December 2009 (30 years)
- **Timescales:** 3, 6, 12, and 24 months

SPEI values were computed for all land cells across the full period (January 1980 – December 2024). Output GeoTIFFs are provided for the target period (January 2010 – December 2024) only; full-period timeseries are available in the NetCDF files.

### Note on reference period SPEI values

SPEI values are provided for the full 1980–2024 period in the NetCDF outputs, 
including the 1980–2009 reference period itself. Values within the reference 
period are self-referential — they describe how each month ranked within the 
log-logistic distribution that was fitted using all months in that same period. 
This is standard practice in SPEI computation and is consistent with the 
approach used in the global SPEIbase dataset (Beguería & Vicente-Serrano 2023). 
However, users should exercise caution when interpreting SPEI values within the 
reference period, as they do not represent anomalies relative to an independent 
baseline. The 2010–2024 monthly GeoTIFF outputs are the primary scientific 
deliverable of this dataset, as these values are genuinely out-of-sample 
relative to the reference period.

### 4. Computational Approach

Due to the continental extent and 0.01° resolution of the dataset (~7 million land cells), computation was performed on the National Computational Infrastructure (NCI) Gadi supercomputer using a spatial tiling approach. The Australian extent was divided into a 3×3 grid of 9 tiles, each processed as an independent job using `terra` (Hijmans 2024) for raster operations. Tile outputs were subsequently merged into continental-scale GeoTIFFs and NetCDF files.

---

## Known Limitations

### -Inf values in SPEI outputs

A subset of cells contain values of `-Inf` (negative infinity). These occur in hyper-arid regions and during extended drought periods where accumulated climatic water balance approaches zero over the accumulation window, causing the log-logistic distribution fitting to fail. This is a known and documented limitation of the SPEI method in arid environments, acknowledged by the original SPEI authors:

> "In some rare cases it was not possible to achieve a good fit to the log-logistic distribution... Most of the problems were in areas that were very arid, at high altitude, and at high latitude." — SPEIbase documentation (spei.csic.es/database.html)

This issue is also documented in the official SPEI R package GitHub repository (github.com/sbegueria/SPEI/issues/11), which notes that -Inf values occur when precipitation values are zero or near zero, particularly at longer accumulation scales.

**Distribution of -Inf values in this dataset:**

| Scale | Months affected | Total -Inf cells | % of all cell-months |
|---|---|---|---|
| SPEI-3 | 60 / 180 (33%) | 2,341,700 | 0.185% |
| SPEI-6 | 43 / 180 (24%) | 904,198 | 0.071% |
| SPEI-12 | 48 / 180 (27%) | 1,675,760 | 0.132% |
| SPEI-24 | 90 / 180 (50%) | 17,839,493 | 1.406% |

The spatial and temporal distribution of these values is climatologically coherent. They cluster around the 2019–2020 Black Summer drought and the 2013–2016 drought — the most extreme multi-year drought events in the Australian instrumental record — and are concentrated in hyper-arid interior regions (e.g. central Australia) where monthly rainfall totals regularly approach zero. This pattern is consistent with the global SPEIbase dataset, which reports the same artefact over the Sahara and other hyper-arid regions.

**These values are retained as-is in the GeoTIFF outputs** to preserve the raw computation and allow users to apply their own treatment. In the NetCDF outputs they are stored as the missing value flag (-9999). Users wishing to exclude or clamp extreme values for downstream analysis may apply the following in R:

```r
library(terra)
r <- rast("SPEI24_202001.tif")
r_clamped <- clamp(r, lower = -3.09, upper = 3.09)  # 0.1th / 99.9th percentile
```

---

## Reproducing This Dataset

All code used to produce this dataset is available at:  
https://github.com/erbendall/spei-australia

### Script Overview

| Script | Route | Purpose |
|---|---|---|
| `01_download_anuclimate.Rmd` | Local | Download ANUClimate v2.0 data via NCI THREDDS (no NCI account required) |
| `02_spei_tile_hpc.R` | HPC | Per-tile SPEI computation — processes one spatial tile (called by PBS array job) |
| `02_spei_merge_hpc.R` | HPC | Merge tile outputs into continental GeoTIFFs and NetCDFs |
| `02_spei_pipeline_local.R` | Local | Full end-to-end pipeline for local machines — computes CWB, SPEI tiles, and final outputs without HPC access |
| `02_spei_local_missing.R` | Local | Targeted script for recomputing specific failed or missing tile-scale combinations |
| `03_spei_merge_local.R` | Local | Merges completed local tile outputs into final GeoTIFFs and NetCDFs |
| `hpc/spei_tile_job.sh` | HPC | PBS array job script — submits 9 independent tile jobs simultaneously |
| `hpc/spei_merge_job.sh` | HPC | PBS merge job script |
| `hpc/spei_submit.sh` | HPC | Master submission script — submits array + merge job with correct dependency |

### Step-by-Step: HPC Route (NCI Gadi or compatible PBS system)

Requires an NCI account with access to project `gh70` (ANUClimate data). Recommended for reproducing at full continental scale.

1. Transfer scripts to Gadi: `scp hpc/*.sh hpc/*.R username@gadi.nci.org.au:~/spei/`
2. Update paths and project code in `02_spei_tile_hpc.R` and `hpc/spei_submit.sh`
3. Submit jobs: `bash hpc/spei_submit.sh`
4. Monitor with `qstat -u username`; logs written to `~/spei/logs/`
5. If any tile-scale combinations fail, run `02_spei_local_missing.R` locally using the CWB layers downloaded from Gadi
6. Once all tiles are present, run `03_spei_merge_local.R` to produce final outputs
7. Download outputs: `scp -r username@gadi.nci.org.au:/g/data/YOURPROJECT/spei_outputs/ .`

### Step-by-Step: Local Machine Route

No NCI account required. Requires ~50 GB disk space and at least 16 GB RAM. Estimated run time: 2–5 days on a modern desktop or laptop. This route uses the same 3×3 tiling approach as the HPC pipeline but processes tiles sequentially on a single machine.

1. Open `01_download_anuclimate.Rmd` in RStudio and follow instructions to download ANUClimate v2.0 data (~35 GB)
2. Update `base_path` and `N_CORES` in `02_spei_pipeline_local.R` to match your local setup
3. Run `02_spei_pipeline_local.R` — computes CWB layers, per-tile SPEI, and final merged outputs in sequence
4. If the pipeline is interrupted, re-run `02_spei_pipeline_local.R` — completed tiles are checkpointed and skipped automatically
5. If specific tile-scale combinations fail, use `02_spei_local_missing.R` to recompute them individually
6. If tiles are complete but final outputs are missing, run `03_spei_merge_local.R` to produce final GeoTIFFs and NetCDFs
7. Outputs are written to `{base_path}/spei/`

### R Package Dependencies

```r
install.packages(c("terra", "SPEI", "ncdf4", "lubridate", "here", "parallel"))
```

| Package | Purpose |
|---|---|
| `terra` | Raster I/O and spatial operations |
| `SPEI` | SPEI computation (log-logistic fitting) |
| `ncdf4` | NetCDF file creation and attribute handling |
| `lubridate` | Date sequence generation |
| `here` | Portable file paths |
| `parallel` | Multi-core processing (local pipeline) |

---

## Citation

If you use this dataset, please cite:

**Dataset:**
> Bendall, E.R. (2026). Gridded SPEI for Australia 1980–2024 at 0.01° resolution derived from ANUClimate v2.0. Zenodo. https://doi.org/10.5281/zenodo.19447029

**Code:**
> Bendall, E.R. (2026). SPEI Australia pipeline (v1.0). GitHub/Zenodo. https://doi.org/10.5281/zenodo.19629231

**Input data:**
> Hutchinson, M., Xu, T., Kesteven, J., Marang, I., and Evans, B. (2021). ANUClimate v2.0 monthly climate data for the Australian continent, 1900 to present. NCI Australia. https://doi.org/10.25914/60a10acd183a2

**SPEI method:**
> Vicente-Serrano, S.M., Beguería, S., López-Moreno, J.I. (2010). A multiscalar drought index sensitive to global warming: the Standardized Precipitation Evapotranspiration Index. *Journal of Climate*, 23(7), 1696–1718. https://doi.org/10.1175/2009JCLI2909.1

**SPEI R package:**
> Beguería, S. & Vicente-Serrano, S.M. (2023). SPEI: Calculation of the Standardised Precipitation-Evapotranspiration Index. R package version 1.8.1. https://cran.r-project.org/package=SPEI

**PET method:**
> Hargreaves, G.H. & Samani, Z.A. (1985). Reference crop evapotranspiration from temperature. *Applied Engineering in Agriculture*, 1(2), 96–99. https://doi.org/10.13031/2013.26773

---

## Acknowledgements

This dataset was produced using resources from the National Computational Infrastructure (NCI), which is supported by the Australian Government through the Western Sydney University partner share. This work was supported by Western Sydney University and with funding from the Australian Government under the National Environmental Science Program's Resilient Landscapes Hub. 

Code developed with AI assistance from Claude (Anthropic). https://claude.ai

---

## Licence

**Dataset:** Creative Commons Attribution 4.0 International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/))  
**Code:** MIT Licence — see [LICENSE](LICENSE)

---

## Contact

Eli R. Bendall  
Hawkesbury Institute for the Environment  
Western Sydney University  
ORCID: [0000-0002-9120-3457](https://orcid.org/0000-0002-9120-3457)  
GitHub: [github.com/erbendall](https://github.com/erbendall)
