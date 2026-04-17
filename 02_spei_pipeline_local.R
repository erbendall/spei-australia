# =============================================================================
# SPEI Pipeline — Full Local Machine Version
# =============================================================================
# Author:      Eli R. Bendall
#              Hawkesbury Institute for the Environment,
#              Western Sydney University, NSW, Australia
#
# AI Partner:  Claude (Anthropic, claude.ai) — assisted in code development,
#              debugging, and pipeline design throughout this project.
#              All scientific decisions and final code review remain the
#              author's responsibility.
#
# Date:        April 2026
# Version:     1.0
#
# Description:
#   Full end-to-end SPEI pipeline for users without HPC access. Computes
#   gridded SPEI at 3, 6, 12, and 24 month timescales for Australia using
#   ANUClimate v2.0 monthly climate data downloaded via 01_download_anuclimate.Rmd.
#
#   The pipeline uses a spatial tiling approach (3x3 grid = 9 tiles) to manage
#   memory, processing one tile at a time using mclapply for parallelisation.
#   Intermediate CWB layers and tile outputs are written to disk so the pipeline
#   can be safely interrupted and resumed without recomputing completed tiles.
#
# Prerequisites:
#   - ANUClimate v2.0 data downloaded via 01_download_anuclimate.Rmd
#   - R packages: terra, SPEI, ncdf4, lubridate, parallel
#
# Input data structure (set base_path below):
#   base_path/
#   ├── rain/   ANUClimate_v2-0_rain_monthly_YYYYMM.nc  (1980-2024)
#   ├── tmax/   ANUClimate_v2-0_tmax_monthly_YYYYMM.nc
#   ├── tmin/   ANUClimate_v2-0_tmin_monthly_YYYYMM.nc
#   └── srad/   ANUClimate_v2-0_srad_monthly_YYYYMM.nc
#
# Outputs (written to base_path/spei/):
#   SPEI3_1980_2024.nc, SPEI6_1980_2024.nc, SPEI12_1980_2024.nc, SPEI24_1980_2024.nc
#   SPEI3_monthly_tifs/, SPEI6_monthly_tifs/, SPEI12_monthly_tifs/, SPEI24_monthly_tifs/
#
# SPEI parameters:
#   PET method:       Hargreaves (Hargreaves & Samani 1985)
#   Distribution:     Log-logistic, ub-pwm fitting
#   Reference period: January 1980 - December 2009
#   Output period:    January 2010 - December 2024
#   Timescales:       3, 6, 12, 24 months
#   CRS:              GDA94 (EPSG:4283), 0.01 degree resolution
#
# Memory and time notes:
#   Estimated peak RAM: 16-24 GB (varies by N_CORES and tile)
#   Estimated run time: 2-5 days on a modern desktop or laptop
#   If RAM is limited (<16 GB), set N_CORES <- 1
#   For interrupted runs: resume by re-running — completed tiles are skipped
#   For failed individual tile-scale combinations, use 02_spei_local_missing.R
#
# Known limitations:
#   SPEI values of -Inf occur in hyper-arid regions and during extreme drought
#   periods where accumulated CWB approaches zero. This is a known limitation
#   of the SPEI method documented by Vicente-Serrano et al. and in the SPEIbase
#   global dataset. See README.md for full details and the -Inf cell count table.
#
# Citation:
#   Bendall, E.R. (2026). Gridded SPEI for Australia 1980-2024 at 0.01 degree
#   resolution derived from ANUClimate v2.0. Zenodo.
#   https://doi.org/10.5281/zenodo.XXXXXXX
#   Code: https://github.com/erbendall/spei-australia
# =============================================================================

library(terra)
library(SPEI)
library(ncdf4)
library(lubridate)
library(parallel)

# ---- USER SETTINGS ----------------------------------------------------------

# Root directory containing downloaded ANUClimate data subfolders
# (rain/, tmax/, tmin/, srad/)
base_path <- "/path/to/anuclimate_data"  # Update to your local data directory

# Number of cores for parallel tile processing
# Reduce to 1 or 2 if RAM is limited (<16 GB)
N_CORES <- 4

# -----------------------------------------------------------------------------

spei_dir  <- file.path(base_path, "spei")
temp_dir  <- file.path(base_path, "terra_temp")
cwb_dir   <- file.path(temp_dir, "cwb_layers")
tiles_dir <- file.path(temp_dir, "tiles")

for (d in c(spei_dir, temp_dir, cwb_dir, tiles_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

terraOptions(tempdir = temp_dir, memfrac = 0.6, progress = 0)

# ---- Date sequence ----------------------------------------------------------

dates      <- seq(as.Date("1980-01-01"), as.Date("2024-12-01"), by = "month")
n_layers   <- length(dates)
out_idx    <- which(dates >= as.Date("2010-01-01"))
out_dates  <- dates[out_idx]
days_in_mo <- days_in_month(dates)

ref_start  <- as.Date("1980-01-01")
ref_end    <- as.Date("2009-12-01")

nc_path <- function(var, d) {
  file.path(base_path, var,
            paste0("ANUClimate_v2-0_", var, "_monthly_",
                   format(d, "%Y%m"), ".nc"))
}

# ---- Define 3x3 tile grid ---------------------------------------------------

# Australian extent: 112-154E, 9-44S
lon_breaks <- seq(112, 154, length.out = 4)
lat_breaks <- seq(-44, -9, length.out = 4)

tiles <- vector("list", 9)
idx   <- 1
for (row in 1:3) {
  for (col in 1:3) {
    tiles[[idx]] <- ext(
      lon_breaks[col], lon_breaks[col + 1],
      lat_breaks[row], lat_breaks[row + 1]
    )
    idx <- idx + 1
  }
}

# ---- Step 1: Compute and write CWB layers -----------------------------------

message("\n=== Step 1: Computing CWB layers ===")

for (i in seq_along(dates)) {
  cwb_file <- file.path(cwb_dir, paste0("cwb_", format(dates[i], "%Y%m"), ".tif"))
  if (file.exists(cwb_file)) next

  r_rain <- rast(nc_path("rain", dates[i]))
  r_tmax <- rast(nc_path("tmax", dates[i]))
  r_tmin <- rast(nc_path("tmin", dates[i]))
  r_srad <- rast(nc_path("srad", dates[i]))

  r_tmean <- (r_tmax + r_tmin) / 2
  r_pet   <- 0.0023 * (r_srad * days_in_mo[i]) * (r_tmean + 17.8) *
             (r_tmax - r_tmin)^0.5
  r_cwb   <- r_rain - r_pet

  writeRaster(r_cwb, cwb_file, datatype = "FLT4S", overwrite = TRUE)

  if (i %% 12 == 0) message("  CWB computed: ", format(dates[i], "%Y-%m"))
  rm(r_rain, r_tmax, r_tmin, r_srad, r_tmean, r_pet, r_cwb)
}
message("CWB layers complete.")

# ---- Step 2: Compute SPEI per tile per scale --------------------------------

message("\n=== Step 2: Computing SPEI tiles ===")

scales <- c(3, 6, 12, 24)

process_tile <- function(tile_idx) {
  tile_ext <- tiles[[tile_idx]]

  # Load all CWB layers for this tile
  cwb_files <- file.path(cwb_dir,
                         paste0("cwb_", format(dates, "%Y%m"), ".tif"))
  r_cwb <- rast(cwb_files)
  r_cwb <- crop(r_cwb, tile_ext, snap = "near")

  # Check for ocean-only tile
  vals <- values(r_cwb[[1]])
  if (all(is.na(vals))) {
    message("  Tile ", tile_idx, ": ocean-only, skipping")
    writeLines("ocean", file.path(tiles_dir,
               paste0("tile_", tile_idx, "_ocean.txt")))
    return(invisible(NULL))
  }

  for (sc in scales) {
    out_file <- file.path(tiles_dir,
                          paste0("tile_", tile_idx, "_SPEI", sc, ".tif"))
    if (file.exists(out_file)) {
      message("  Tile ", tile_idx, " SPEI-", sc, ": already done, skipping")
      next
    }

    message("  Tile ", tile_idx, " SPEI-", sc, ": computing...")

    cwb_mat <- values(r_cwb)  # cells x layers matrix

    spei_mat <- matrix(NA_real_, nrow = nrow(cwb_mat), ncol = n_layers)

    for (cell in seq_len(nrow(cwb_mat))) {
      x <- cwb_mat[cell, ]
      if (all(is.na(x))) next
      spei_mat[cell, ] <- tryCatch(
        as.numeric(spei(
          ts(x, start = c(1980, 1), frequency = 12),
          scale = sc,
          ref.start = c(1980, 1),
          ref.end   = c(2009, 12),
          na.rm     = TRUE
        )$fitted),
        error = function(e) rep(NA_real_, n_layers)
      )
    }

    # Write output TIF (all layers, full period)
    r_spei <- setValues(r_cwb, spei_mat)
    writeRaster(r_spei, out_file, datatype = "FLT4S", overwrite = TRUE)
    message("  Tile ", tile_idx, " SPEI-", sc, ": done")
    rm(spei_mat, r_spei)
    gc()
  }

  rm(r_cwb, cwb_mat)
  gc()
  invisible(NULL)
}

mclapply(seq_along(tiles), process_tile, mc.cores = N_CORES)

message("All tiles complete.")

# ---- Step 3: Merge tiles and write final outputs ----------------------------

message("\n=== Step 3: Merging tiles and writing outputs ===")

for (sc in scales) {
  message("  Merging SPEI-", sc, "...")

  # Collect tile files (skip ocean tiles)
  tile_files <- file.path(tiles_dir,
                           paste0("tile_", seq_along(tiles), "_SPEI", sc, ".tif"))
  tile_files <- tile_files[file.exists(tile_files)]

  if (length(tile_files) == 0) {
    warning("No tile files found for SPEI-", sc, " — skipping")
    next
  }

  # Write monthly TIFs for output period
  tif_dir <- file.path(spei_dir, paste0("SPEI", sc, "_monthly_tifs"))
  dir.create(tif_dir, showWarnings = FALSE)

  for (i in seq_along(out_idx)) {
    layer_idx <- out_idx[i]
    out_tif   <- file.path(tif_dir,
                            paste0("SPEI", sc, "_", format(out_dates[i], "%Y%m"), ".tif"))
    if (file.exists(out_tif)) next

    tile_layers <- lapply(tile_files, function(f) rast(f)[[layer_idx]])
    merged      <- do.call(merge, tile_layers)
    writeRaster(merged, out_tif, datatype = "FLT4S", overwrite = TRUE)
    rm(tile_layers, merged)
  }

  message("    Monthly TIFs written: SPEI-", sc)

  # Write full-period NetCDF (1980-2024)
  nc_out <- file.path(spei_dir, paste0("SPEI", sc, "_1980_2024.nc"))
  if (!file.exists(nc_out)) {
    # Merge all layers month by month and build NetCDF
    r_first    <- rast(tile_files[1])
    lon_vals   <- xFromCol(r_first, 1:ncol(r_first))
    lat_vals   <- yFromRow(r_first, 1:nrow(r_first))

    lon_dim  <- ncdim_def("lon",  "degrees_east",  lon_vals)
    lat_dim  <- ncdim_def("lat",  "degrees_north", lat_vals)
    time_dim <- ncdim_def("time", "months since 1980-01-01",
                           as.numeric(dates - as.Date("1980-01-01")) / 30.4375,
                           unlim = TRUE)

    spei_var <- ncvar_def(
      name     = paste0("SPEI", sc),
      units    = "unitless",
      dim      = list(lon_dim, lat_dim, time_dim),
      missval  = -9999,
      longname = paste0("Standardised Precipitation Evapotranspiration Index — ",
                        sc, "-month timescale")
    )

    nc <- nc_create(nc_out, spei_var)
    ncatt_put(nc, 0, "title",
              paste0("Gridded SPEI-", sc, " for Australia 1980-2024"))
    ncatt_put(nc, 0, "institution",
              "Hawkesbury Institute for the Environment, Western Sydney University")
    ncatt_put(nc, 0, "source", "ANUClimate v2.0 (doi:10.25914/60a10acd183a2)")
    ncatt_put(nc, 0, "reference_period", "1980-01-01 to 2009-12-31")
    ncatt_put(nc, 0, "PET_method", "Hargreaves (Hargreaves & Samani 1985)")
    ncatt_put(nc, 0, "missing_value", -9999)

    for (i in seq_len(n_layers)) {
      tile_layers <- lapply(tile_files, function(f) rast(f)[[i]])
      merged      <- do.call(merge, tile_layers)
      layer_vals  <- values(merged)[, 1]
      layer_vals[is.infinite(layer_vals)] <- -9999
      layer_vals[is.nan(layer_vals)]      <- -9999
      ncvar_put(nc, spei_var, layer_vals,
                start = c(1, 1, i), count = c(-1, -1, 1))
      rm(tile_layers, merged)
    }

    nc_close(nc)
    message("    NetCDF written: ", basename(nc_out))
  }
}

message("\n=== Pipeline complete ===")
message("Outputs written to: ", spei_dir)
