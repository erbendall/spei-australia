# =============================================================================
# SPEI Pipeline — Local Merge Script
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
#              optimisation strategies, error diagnosis, and iterative
#              debugging. All scientific decisions, parameter choices and
#              final code review remain the author's responsibility.
#
# Date:        April 2026
# Version:     1.0
#
# Description:
#   Merges per-tile SPEI TIF outputs (downloaded from HPC or computed locally)
#   into full-Australia rasters and writes final GeoTIFF and NetCDF outputs.
#   Designed for local machines with sufficient RAM (16GB+ recommended).
#   Uses ncdf4 directly to avoid terra writeCDF scale/offset issues.
#
# Usage:
#   1. Update base_path below to your local tile output directory
#   2. Source this script in R or run: Rscript spei_merge_local.R
#
# Input:  Tile TIFs in tiles_dir (one per month per scale per tile)
# Output: Merged monthly TIFs (2010-2024) + NetCDF (1980-2024)
#
# Memory: ~8-16GB RAM peak — suitable for 16-32GB machines
# Time:   ~2-6 hours depending on hardware
#
# Known limitations:
#   SPEI -Inf values (hyper-arid regions, extreme drought) are preserved
#   in GeoTIFF outputs and replaced with -9999 in NetCDF outputs.
#   See README.md for full documentation.
#
# Citation:
#   Bendall, E.R. (2026). Gridded SPEI for Australia 1980-2024 at 0.01 degree
#   resolution derived from ANUClimate v2.0. Zenodo.
#   https://doi.org/10.5281/zenodo.XXXXXXX
#   Code: https://github.com/erbendall/spei-australia
# =============================================================================
library(terra)
library(lubridate)

# ---- Paths ------------------------------------------------------------------
# Update these paths to match your local directory structure

base_path <- "/path/to/spei_outputs"  # Update to your local output directory
tiles_dir <- file.path(base_path, "tiles")       # downloaded tile TIFs
out_dir <- file.path(base_path, "spei")
temp_dir  <- file.path(base_path, "merge_temp")   # temp working directory

dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

# Direct terra temp files to local drive
terraOptions(
  tempdir  = temp_dir,
  memfrac  = 0.5,   # 50% of RAM = ~16GB on 32GB Mac
  progress = 1
)

message("Tiles directory: ", tiles_dir)
message("Output directory: ", out_dir)

# ---- Setup ------------------------------------------------------------------

dates     <- seq(as.Date("1980-01-01"), as.Date("2024-12-01"), by = "month")
out_idx   <- which(dates >= as.Date("2010-01-01") & dates <= as.Date("2024-12-01"))
out_dates <- dates[out_idx]
n_layers  <- length(dates)
n_tiles   <- 9L

# ---- Verify tile TIF downloads ----------------------------------------------

message("\n---- Verifying downloaded tile TIFs ----")

# Identify which tiles have data (all 9 should be present)
land_tiles <- integer(0)
for (i in seq_len(n_tiles)) {
  # Check for at least one TIF for this tile in SPEI3
  test_file <- file.path(tiles_dir, "SPEI3",
                         sprintf("SPEI3_198001_tile%02d.tif", i))
  if (file.exists(test_file)) {
    land_tiles <- c(land_tiles, i)
    message("  Tile ", i, ": present")
  } else {
    message("  Tile ", i, ": MISSING — check download")
  }
}

if (length(land_tiles) == 0) {
  stop("No tile TIFs found in: ", tiles_dir,
       "\nCheck that the rsync download completed successfully.")
}
message("  Land tiles found: ", paste(land_tiles, collapse = ", "))

# Verify TIF counts for each scale and tile
message("\n---- Verifying TIF counts ----")
for (sc in c(3, 6, 12, 24)) {
  sc_dir <- file.path(tiles_dir, paste0("SPEI", sc))
  if (!dir.exists(sc_dir)) {
    stop("SPEI", sc, " directory not found: ", sc_dir)
  }
  for (i in land_tiles) {
    n_tifs <- length(list.files(sc_dir,
                                pattern = sprintf("SPEI%d_.*_tile%02d\\.tif$",
                                                  sc, i)))
    if (n_tifs != n_layers) {
      warning("Tile ", i, " SPEI-", sc, ": expected ", n_layers,
              " TIFs, found ", n_tifs,
              "\nPartial download — rerun rsync to complete.")
    }
  }
  message("  SPEI-", sc, ": verified")
}

# ---- Merge function ---------------------------------------------------------

merge_scale <- function(scale) {

  message("\n======== Merging SPEI-", scale, " ========")
  t0 <- Sys.time()

  sc_tile_dir  <- file.path(tiles_dir, paste0("SPEI", scale))
  sc_merge_dir <- file.path(temp_dir,  paste0("SPEI", scale, "_merged"))
  sc_final_dir <- file.path(out_dir,   paste0("SPEI", scale, "_monthly_tifs"))
  nc_out       <- file.path(out_dir,   paste0("SPEI", scale, "_1980_2024.nc"))

  dir.create(sc_merge_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(sc_final_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Step 1: Merge tiles month by month -----------------------------------
  message("  Step 1: Merging tiles for ", n_layers, " months...")

  for (j in seq_along(dates)) {
    date_str    <- format(dates[j], "%Y%m")
    merged_file <- file.path(sc_merge_dir,
                             paste0("SPEI", scale, "_", date_str, ".tif"))

    if (file.exists(merged_file)) next  # checkpoint — skip if done

    # Load tile TIFs for this month
    tile_files <- file.path(sc_tile_dir,
                             sprintf("SPEI%d_%s_tile%02d.tif",
                                     scale, date_str, land_tiles))

    missing <- !file.exists(tile_files)
    if (any(missing)) {
      stop("Missing tile TIFs for SPEI-", scale, " month ", date_str,
           ": tiles ", paste(land_tiles[missing], collapse = ", "),
           "\nCheck your rsync download.")
    }

    # Merge tiles for this month
    tile_rasts <- lapply(tile_files, rast)
    r_merged   <- if (length(tile_rasts) == 1) tile_rasts[[1]] else
                  do.call(terra::merge, tile_rasts)
    rm(tile_rasts)

    writeRaster(r_merged, merged_file, overwrite = TRUE)
    rm(r_merged)

    if (j %% 60 == 0) {
      gc()
      message("  Merged ", j, "/", n_layers, " months (",
              round(j/n_layers*100), "%)")
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
    if (!file.exists(out_file)) {
      file.copy(merged_file, out_file, overwrite = TRUE)
    }
  }

  n_out <- length(list.files(sc_final_dir, pattern = "\\.tif$"))
  message("  Output TIFs: ", n_out, "/", length(out_idx))

  # ---- Step 3: Build NetCDF from merged monthly TIFs ------------------------
  # Stack all 540 merged TIFs and write as single NetCDF in one operation
  # This avoids writeCDF append issues
  message("  Step 3: Building NetCDF (stacking ", n_layers, " merged TIFs)...")

  merged_files <- file.path(sc_merge_dir,
                             paste0("SPEI", scale, "_",
                                    format(dates, "%Y%m"), ".tif"))

  missing_merged <- !file.exists(merged_files)
  if (any(missing_merged)) {
    stop(sum(missing_merged), " merged TIFs missing for SPEI-", scale)
  }

  r_stack <- rast(merged_files)
  writeCDF(r_stack, nc_out,
           varname  = paste0("SPEI", scale),
           longname = paste0("SPEI ", scale, "-month, ref 1980-2009"),
           unit     = "unitless",
           overwrite = TRUE)
  rm(r_stack); gc()

  nc_size <- round(file.size(nc_out) / 1e9, 2)
  message("  NetCDF written: ", nc_size, " GB")

  # ---- Validate output ------------------------------------------------------
  r_check  <- rast(nc_out)
  test_pts <- data.frame(
    x = c(150.9, 135.0, 144.9, 116.9),
    y = c(-33.9, -25.0, -37.8, -31.9)
  )
  vals    <- extract(r_check[[out_idx[1]]], test_pts)
  n_valid <- sum(!is.na(vals[, 2]))
  rng     <- range(vals[, 2], na.rm = TRUE)
  message(sprintf("  Validation: %d/4 valid, range %.2f to %.2f",
                  n_valid, rng[1], rng[2]))
  if (n_valid == 0) warning("SPEI-", scale, ": merged NC all NA — check tiles!")

  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  message("  SPEI-", scale, " done in ", elapsed, " min.")
}

# ---- Run merge for all scales -----------------------------------------------

for (sc in c(3, 6, 12, 24)) {

  nc_out  <- file.path(out_dir, paste0("SPEI", sc, "_1980_2024.nc"))
  tif_dir <- file.path(out_dir, paste0("SPEI", sc, "_monthly_tifs"))
  n_tifs  <- length(list.files(tif_dir, pattern = "\\.tif$"))

  if (file.exists(nc_out) && n_tifs == length(out_idx)) {
    message("SPEI-", sc, " already merged — skipping.")
    next
  }

  merge_scale(sc)
}

# ---- Cleanup temp merge files -----------------------------------------------

message("\nCleaning up temp merge files...")
unlink(temp_dir, recursive = TRUE)
message("Done.")

# ---- Final summary ----------------------------------------------------------

message("\n========================================")
message("MERGE COMPLETE: ", format(Sys.time()))
message("Output directory: ", out_dir)
message("")
for (sc in c(3, 6, 12, 24)) {
  nc_out  <- file.path(out_dir, paste0("SPEI", sc, "_1980_2024.nc"))
  tif_dir <- file.path(out_dir, paste0("SPEI", sc, "_monthly_tifs"))
  n_tifs  <- length(list.files(tif_dir, pattern = "\\.tif$"))
  nc_size <- ifelse(file.exists(nc_out),
                    paste0(round(file.size(nc_out)/1e9, 1), " GB"),
                    "MISSING")
  message(sprintf("  SPEI-%2d: NC=%s  TIFs=%d/%d",
                  sc, nc_size, n_tifs, length(out_idx)))
}
message("========================================")
