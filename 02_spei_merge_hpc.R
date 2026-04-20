# =============================================================================
# SPEI Pipeline — HPC Merge Script (Generalised)
# =============================================================================
# Author:      Eli R. Bendall
#              Hawkesbury Institute for the Environment,
#              Western Sydney University, NSW, Australia
#
# Funding:     This work was supported by Western Sydney University and with
#              funding from the Australian Government under the National 
#              Environmental Science Program's Resilient Landscapes Hub.
#
# AI Partner:  Claude (Anthropic, claude.ai) — assisted in code development,
#              debugging, and pipeline design throughout this project.
#              This script was developed collaboratively between the author
#              and Claude over multiple working sessions in March-April 2026.
#              Claude's contributions included: pipeline architecture, memory
#              optimisation strategies, PBS job scripting, error diagnosis,
#              and iterative debugging. All scientific decisions, parameter
#              choices and final code review remain the author's responsibility.
#
# Date:        April 2026
# Version:     1.0
#
# Description:
#   Merges per-tile SPEI TIF outputs from the HPC tile pipeline into
#   full-Australia rasters and writes final NetCDF and GeoTIFF outputs.
#
#   This script is called by spei_merge_job.sh after all tile jobs complete.
#   It can also be run manually after verifying all tiles are complete.
#
# Input:
#   Per-tile SPEI TIF files written by spei_tile.R
#   Expected location: OUTPUT_DIR/tiles/SPEI{scale}/SPEI{scale}_{YYYYMM}_tile{NN}.tif
#
# Output:
#   - Monthly GeoTIFFs: OUTPUT_DIR/SPEI{scale}_monthly_tifs/SPEI{scale}_{YYYYMM}.tif
#   - Full NetCDF:      OUTPUT_DIR/SPEI{scale}_1980_2024.nc
#
# Known limitations:
#   SPEI values of -Inf occur in hyper-arid regions during extreme drought
#   periods where log-logistic distribution fitting fails. These are retained
#   as -Inf in GeoTIFF outputs and stored as missing values (-9999) in NetCDF
#   outputs. See README.md for full documentation and references.
#
# Citation:
#   Bendall, E.R. (2026). Gridded SPEI for Australia 1980-2024 at 0.01 degree
#   resolution derived from ANUClimate v2.0. Zenodo.
#   https://doi.org/10.5281/zenodo.XXXXXXX
#   Code: https://github.com/erbendall/spei-australia
#
# Known limitations:
#   SPEI values of -Inf occur in hyper-arid regions and during extreme drought
#   periods where log-logistic distribution fitting fails. These are retained
#   as -Inf in GeoTIFF outputs and stored as missing values (-9999) in NetCDF
#   outputs. See README.md for full documentation and references.
#
#   SPEI values within the reference period (January 1980 - December 2009)
#   are self-referential and should be interpreted with caution. The 2010-2024
#   monthly GeoTIFF outputs are the primary scientific deliverable.
# =============================================================================

library(terra)
library(ncdf4)
library(lubridate)

# ---- Paths — set via environment variables from PBS job script --------------
# These are passed from spei_merge_job.sh via -v flag
# Can also be set manually for standalone use

output_dir <- Sys.getenv("SPEI_OUTPUT_DIR",
                          unset = "/scratch/YOUR_PROJECT/spei_outputs")
scratch    <- Sys.getenv("SPEI_SCRATCH_DIR",
                          unset = "/scratch/YOUR_PROJECT/spei_temp")

tile_dir   <- file.path(output_dir, "tiles")
final_dir  <- output_dir
merged_dir <- file.path(scratch, "merge_temp")

jobfs <- Sys.getenv("PBS_JOBFS", unset = scratch)
terraOptions(tempdir = jobfs, memfrac = 0.6, progress = 0)

message("========================================")
message("SPEI Merge Script")
message("Start: ", format(Sys.time()))
message("Output: ", final_dir)
message("========================================")

# ---- Setup ------------------------------------------------------------------

dates     <- seq(as.Date("1980-01-01"), as.Date("2024-12-01"), by = "month")
out_idx   <- which(dates >= as.Date("2010-01-01") & dates <= as.Date("2024-12-01"))
out_dates <- dates[out_idx]
n_layers  <- length(dates)
n_tiles   <- 9L

# ---- Check tile completion --------------------------------------------------

message("\n---- Checking tile completion ----")
incomplete <- integer(0)
ocean      <- integer(0)
land_tiles <- integer(0)

for (i in seq_len(n_tiles)) {
  complete_file <- file.path(tile_dir, sprintf("tile%02d_complete.txt", i))
  ocean_file    <- file.path(tile_dir, sprintf("tile%02d_ocean.txt", i))

  if (file.exists(complete_file)) {
    message("  Tile ", i, ": complete")
    land_tiles <- c(land_tiles, i)
  } else if (file.exists(ocean_file)) {
    message("  Tile ", i, ": ocean — skipping")
    ocean <- c(ocean, i)
  } else {
    # Allow manual sentinel files for tiles killed before writing sentinel
    tif_check <- file.path(tile_dir, "SPEI3",
                           sprintf("SPEI3_198001_tile%02d.tif", i))
    if (file.exists(tif_check)) {
      message("  Tile ", i, ": complete (manual verification)")
      land_tiles <- c(land_tiles, i)
    } else {
      message("  Tile ", i, ": MISSING")
      incomplete <- c(incomplete, i)
    }
  }
}

if (length(incomplete) > 0) {
  stop("Tiles not yet complete: ", paste(incomplete, collapse = ", "),
       "\nResubmit failed tiles before running merge.")
}
message("Land tiles: ", paste(land_tiles, collapse = ", "))

# ---- Verify tile TIF counts -------------------------------------------------

message("\n---- Verifying tile TIF counts ----")
for (sc in c(3, 6, 12, 24)) {
  for (i in land_tiles) {
    sc_dir <- file.path(tile_dir, paste0("SPEI", sc))
    n_tifs <- length(list.files(sc_dir,
                                pattern = sprintf("SPEI%d_.*_tile%02d\\.tif$",
                                                  sc, i)))
    if (n_tifs != n_layers) {
      stop("Tile ", i, " SPEI-", sc, ": expected ", n_layers,
           " TIFs, found ", n_tifs)
    }
  }
  message("  SPEI-", sc, ": all tile TIFs verified")
}

# ---- Get spatial reference --------------------------------------------------

ref_tif <- rast(file.path(tile_dir, "SPEI3",
                           sprintf("SPEI3_198001_tile%02d.tif",
                                   land_tiles[1])))
# Load first merged month to get full extent
message("Building spatial reference from tiles...")
test_tiles <- lapply(file.path(tile_dir, "SPEI3",
                                sprintf("SPEI3_198001_tile%02d.tif",
                                        land_tiles)), rast)
r_ref  <- do.call(terra::merge, test_tiles)
nrows  <- nrow(r_ref)
ncols  <- ncol(r_ref)
lons   <- xFromCol(r_ref, 1:ncols)
lats   <- yFromRow(r_ref, 1:nrows)
rm(test_tiles, r_ref, ref_tif); gc()
message("  Grid: ", nrows, " rows x ", ncols, " cols")

# ---- Merge function ---------------------------------------------------------

merge_scale <- function(scale) {

  message("\n======== Merging SPEI-", scale, " ========")
  t0 <- Sys.time()

  sc_tile_dir  <- file.path(tile_dir,   paste0("SPEI", scale))
  sc_merge_dir <- file.path(merged_dir, paste0("SPEI", scale))
  sc_final_dir <- file.path(final_dir,  paste0("SPEI", scale, "_monthly_tifs"))
  nc_out       <- file.path(final_dir,  paste0("SPEI", scale, "_1980_2024.nc"))

  dir.create(sc_merge_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(sc_final_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Step 1: Merge tiles month by month -----------------------------------
  message("  Step 1: Merging ", n_layers, " months...")

  for (j in seq_along(dates)) {
    date_str    <- format(dates[j], "%Y%m")
    merged_file <- file.path(sc_merge_dir,
                             paste0("SPEI", scale, "_", date_str, ".tif"))
    if (file.exists(merged_file)) next

    tile_files <- file.path(sc_tile_dir,
                             sprintf("SPEI%d_%s_tile%02d.tif",
                                     scale, date_str, land_tiles))
    missing <- !file.exists(tile_files)
    if (any(missing)) {
      stop("Missing tile TIFs for SPEI-", scale, " month ", date_str,
           ": tiles ", paste(land_tiles[missing], collapse = ", "))
    }

    tile_rasts <- lapply(tile_files, rast)
    r_merged   <- if (length(tile_rasts) == 1) tile_rasts[[1]] else
                  do.call(terra::merge, tile_rasts)
    rm(tile_rasts)
    writeRaster(r_merged, merged_file, overwrite = TRUE)
    rm(r_merged)

    if (j %% 60 == 0) {
      gc()
      message("  Merged ", j, "/", n_layers, " months")
    }
  }

  n_merged <- length(list.files(sc_merge_dir, pattern = "\\.tif$"))
  message("  Merged TIFs: ", n_merged, "/", n_layers)

  # ---- Step 2: Write output TIFs for 2010-2024 ------------------------------
  message("  Step 2: Writing output TIFs (2010-2024)...")

  for (i in seq_along(out_idx)) {
    date_str    <- format(out_dates[i], "%Y%m")
    merged_file <- file.path(sc_merge_dir,
                             paste0("SPEI", scale, "_", date_str, ".tif"))
    out_file    <- file.path(sc_final_dir,
                             paste0("SPEI", scale, "_", date_str, ".tif"))
    if (!file.exists(out_file)) file.copy(merged_file, out_file)
  }

  n_out <- length(list.files(sc_final_dir, pattern = "\\.tif$"))
  message("  Output TIFs: ", n_out, "/", length(out_idx))

  # ---- Step 3: Write NetCDF using ncdf4 -------------------------------------
  # Uses ncdf4 directly to avoid writeCDF scale/offset issues.
  # -Inf values are replaced with missing value (-9999) in NetCDF output.
  # Original -Inf values are preserved in the GeoTIFF outputs.
  message("  Step 3: Writing NetCDF (", n_layers, " time steps)...")

  dim_lon  <- ncdim_def("lon", "degrees_east",  lons)
  dim_lat  <- ncdim_def("lat", "degrees_north", lats)
  dim_time <- ncdim_def("time", "months since 1980-01-01",
                        0:(n_layers - 1L), unlim = TRUE)

  var_spei <- ncvar_def(
    name     = paste0("SPEI", scale),
    units    = "unitless",
    dim      = list(dim_lon, dim_lat, dim_time),
    missval  = -9999,
    longname = paste0("SPEI ", scale, "-month timescale, ",
                      "reference period 1980-2009"),
    prec     = "float"
  )

  nc <- nc_create(nc_out, list(var_spei), force_v4 = TRUE)
  ncatt_put(nc, 0, "title",
            paste0("SPEI-", scale, " for Australia 1980-2024"))
  ncatt_put(nc, 0, "source",
            "ANUClimate v2.0 (doi:10.25914/60a10acd183a2)")
  ncatt_put(nc, 0, "reference_period", "January 1980 - December 2009")
  ncatt_put(nc, 0, "PET_method", "Hargreaves")
  ncatt_put(nc, 0, "distribution", "log-logistic")
  ncatt_put(nc, 0, "missing_value_note",
            paste("-9999 indicates missing or failed distribution fit.",
                  "-Inf values in GeoTIFF outputs indicate extreme drought",
                  "where log-logistic fitting failed (see README.md)."))
  ncatt_put(nc, 0, "author",
            "Eli R. Bendall, Western Sydney University")
  ncatt_put(nc, 0, "orcid", "0000-0002-9120-3457")
  ncatt_put(nc, 0, "code",
            "https://github.com/erbendall/spei-australia")

  for (j in seq_along(dates)) {
    date_str    <- format(dates[j], "%Y%m")
    merged_file <- file.path(sc_merge_dir,
                             paste0("SPEI", scale, "_", date_str, ".tif"))
    r_j  <- rast(merged_file)
    mat  <- values(r_j)
    mat[is.nan(mat)]      <- -9999
    mat[is.na(mat)]       <- -9999
    mat[is.infinite(mat)] <- -9999
    mat  <- matrix(mat, nrow = nrows, ncol = ncols)
    ncvar_put(nc, var_spei, mat,
              start = c(1, 1, j), count = c(ncols, nrows, 1))
    rm(r_j, mat)
    if (j %% 60 == 0) { gc(); message("  Written ", j, "/", n_layers) }
  }

  nc_close(nc)
  nc_size <- round(file.size(nc_out) / 1e9, 2)
  message("  NetCDF: ", nc_size, " GB")

  # ---- Validate output ------------------------------------------------------
  r_check  <- rast(nc_out)
  test_pts <- data.frame(
    x = c(150.9, 135.0, 144.9, 116.9),
    y = c(-33.9, -25.0, -37.8, -31.9)
  )
  vals    <- terra::extract(r_check[[out_idx[1]]], test_pts)
  n_valid <- sum(!is.na(vals[, 2]) & vals[, 2] != -9999)
  rng     <- range(vals[, 2][vals[, 2] != -9999], na.rm = TRUE)
  message(sprintf("  Validation: %d/4 valid, range %.2f to %.2f",
                  n_valid, rng[1], rng[2]))
  if (n_valid == 0) warning("SPEI-", scale, ": all NA — check tile outputs!")
  rm(r_check); gc()

  # Clean up temp merged TIFs
  unlink(sc_merge_dir, recursive = TRUE)

  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  message("  SPEI-", scale, " complete in ", elapsed, " min.")
}

# ---- Run merge for all scales -----------------------------------------------

for (sc in c(3, 6, 12, 24)) {
  nc_out  <- file.path(final_dir, paste0("SPEI", sc, "_1980_2024.nc"))
  tif_dir <- file.path(final_dir, paste0("SPEI", sc, "_monthly_tifs"))
  n_tifs  <- length(list.files(tif_dir, pattern = "\\.tif$"))

  if (file.exists(nc_out) && n_tifs == length(out_idx)) {
    message("SPEI-", sc, " already complete — skipping.")
    next
  }

  merge_scale(sc)
}

# ---- Final summary ----------------------------------------------------------

message("\n========================================")
message("MERGE COMPLETE: ", format(Sys.time()))
message("Output: ", final_dir)
message("")
for (sc in c(3, 6, 12, 24)) {
  nc_out  <- file.path(final_dir, paste0("SPEI", sc, "_1980_2024.nc"))
  tif_dir <- file.path(final_dir, paste0("SPEI", sc, "_monthly_tifs"))
  n_tifs  <- length(list.files(tif_dir, pattern = "\\.tif$"))
  nc_size <- ifelse(file.exists(nc_out),
                    paste0(round(file.size(nc_out)/1e9, 1), " GB"),
                    "MISSING")
  message(sprintf("  SPEI-%2d: NC=%s  TIFs=%d/%d",
                  sc, nc_size, n_tifs, length(out_idx)))
}
message("========================================")
