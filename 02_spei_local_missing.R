# =============================================================================
# SPEI Pipeline — Local Missing Tile-Scale Script
# =============================================================================
# Author:      Eli R. Bendall
#              Hawkesbury Institute for the Environment,
#              Western Sydney University, NSW, Australia
#
# AI Partner:  Claude (Anthropic, claude.ai) — assisted in code development,
#              debugging, and pipeline design throughout this project.
#              This script was developed collaboratively between the author
#              and Claude over multiple working sessions in March-April 2026.
#              Claude's contributions included: pipeline architecture, memory
#              optimisation strategies, parallel processing design, error
#              diagnosis, and iterative debugging. All scientific decisions,
#              parameter choices and final code review remain the author's
#              responsibility.
#
# Date:        April 2026
# Version:     1.0
#
# Description:
#   Computes SPEI for specific missing tile-scale combinations on a local
#   machine using existing CWB layers. Designed for cases where HPC jobs
#   were killed before completing all scales for a tile.
#
#   Processes ONE tile and ONE scale at a time to manage memory.
#   Supports parallel computation via mclapply (2-4 cores recommended).
#
# Usage:
#   1. Set TILE_INDEX, SCALE and N_CORES below
#   2. Source in R or run: Rscript spei_local_missing.R
#   3. Monitor RAM — if >28GB with 4 cores, drop to N_CORES <- 2
#   4. Re-source after any interruption — crop step is checkpointed
#
# Missing combinations to run in order (example):
#   SPEI-6:  tiles 4, 5, 6
#   SPEI-12: tiles 2, 4, 5, 6
#   SPEI-24: tiles 2, 4, 5, 6, 9
#
# Memory: ~8-20GB peak (varies by N_CORES and tile size)
# Time:   ~4-15 hours per combination depending on hardware and N_CORES
#
# Citation:
#   Bendall, E.R. (2026). Gridded SPEI for Australia 1980-2024 at 0.01 degree
#   resolution derived from ANUClimate v2.0. Zenodo.
#   https://doi.org/10.5281/zenodo.XXXXXXX
#   Code: https://github.com/erbendall/spei-australia
# =============================================================================
library(terra)
library(SPEI)
library(lubridate)
library(parallel)

# ---- USER SETTINGS ----------------------------------------------------------
# Change these for each run

TILE_INDEX <- 4    # tile to process (1-9)
SCALE      <- 6    # SPEI scale (6, 12, or 24)
N_CORES    <- 4    # cores for parallel SPEI (try 4 first, drop to 2 if RAM > 28GB)

# -----------------------------------------------------------------------------

# ---- Paths ------------------------------------------------------------------

base_path <- "/path/to/spei_outputs"      # Update to your local output directory
cwb_path  <- file.path(base_path, "terra_temp", "cwb_layers")
tiles_dir <- file.path(base_path, "tiles")
temp_dir  <- file.path(base_path, "merge_temp")

sc_dir <- file.path(tiles_dir, paste0("SPEI", SCALE))
dir.create(sc_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

terraOptions(tempdir = temp_dir, memfrac = 0.4, progress = 0)

message("========================================")
message("SPEI Local — Tile ", TILE_INDEX, " Scale ", SCALE)
message("Cores: ", N_CORES)
message("Start: ", format(Sys.time()))
message("========================================")

# ---- Check if already complete ----------------------------------------------

test_file <- file.path(sc_dir,
                       sprintf("SPEI%d_198001_tile%02d.tif", SCALE, TILE_INDEX))
if (file.exists(test_file)) {
  message("Tile ", TILE_INDEX, " SPEI-", SCALE, " already complete — skipping.")
  message("Change TILE_INDEX or SCALE and rerun.")
  stop("Already complete.")
}

# ---- Setup ------------------------------------------------------------------

dates    <- seq(as.Date("1980-01-01"), as.Date("2024-12-01"), by = "month")
n_layers <- length(dates)
n_time   <- n_layers

# ---- Define tile extent -----------------------------------------------------

lon_breaks <- seq(112, 154, length.out = 4)
lat_breaks <- seq(-44,  -9, length.out = 4)

n_cols  <- 3L
col_idx <- ((TILE_INDEX - 1L) %% n_cols) + 1L
row_idx <- ((TILE_INDEX - 1L) %/% n_cols) + 1L

xmin <- lon_breaks[col_idx]
xmax <- lon_breaks[col_idx + 1L]
ymin <- lat_breaks[4L - row_idx]
ymax <- lat_breaks[4L - row_idx + 1L]

tile_ext <- ext(xmin, xmax, ymin, ymax)
message(sprintf("Tile %d extent: lon %.2f-%.2f, lat %.2f-%.2f",
                TILE_INDEX, xmin, xmax, ymin, ymax))

# ---- Verify CWB layers ------------------------------------------------------

cwb_files <- file.path(cwb_path, paste0("cwb_", format(dates, "%Y%m"), ".tif"))
n_missing <- sum(!file.exists(cwb_files))
if (n_missing > 0) {
  stop(n_missing, " CWB layer files missing from:\n  ", cwb_path)
}
message("CWB layers: ", n_layers, " files verified")

# ---- Crop CWB to tile — one layer at a time ---------------------------------
# This step is skipped on restart if tile layers already exist

tile_cwb_dir   <- file.path(temp_dir, sprintf("tile%02d_cwb", TILE_INDEX))
tile_cwb_files <- file.path(tile_cwb_dir, paste0("cwb_", format(dates, "%Y%m"), ".tif"))
n_existing     <- sum(file.exists(tile_cwb_files))

dir.create(tile_cwb_dir, showWarnings = FALSE, recursive = TRUE)

if (n_existing == n_layers) {
  message("All ", n_layers, " tile CWB layers already exist — skipping crop.")
} else {
  message("Cropping CWB to tile (", n_existing, "/", n_layers, " exist)...")

  r_ref            <- rast(cwb_files[1])
  tile_ext_snapped <- align(tile_ext, r_ref)
  rm(r_ref); gc()

  for (i in seq_len(n_layers)) {
    if (file.exists(tile_cwb_files[i])) next
    r_i <- crop(rast(cwb_files[i]), tile_ext_snapped)
    writeRaster(r_i, tile_cwb_files[i], overwrite = FALSE)
    rm(r_i)
    if (i %% 60 == 0) {
      gc()
      message("  Cropped ", i, "/", n_layers)
    }
  }
  message("All tile CWB layers cropped.")
}

# ---- Load tile as matrix ----------------------------------------------------

message("Loading tile CWB into matrix...")
r_cwb_tile <- rast(tile_cwb_files)
n_land     <- global(!is.na(r_cwb_tile[[1]]), "sum")[1, 1]
message("Land cells: ", format(round(n_land), big.mark = ","))

cwb_matrix <- values(r_cwb_tile)
r_template <- rast(tile_cwb_files[1])
rm(r_cwb_tile); gc()

message("Matrix: ", nrow(cwb_matrix), " x ", ncol(cwb_matrix),
        " (", round(object.size(cwb_matrix) / 1e6, 1), " MB)")

land_idx <- which(!is.na(cwb_matrix[, 1]))
n_cells  <- nrow(cwb_matrix)
message("Processing ", format(length(land_idx), big.mark = ","),
        " land cells using ", N_CORES, " cores")

# ---- Compute SPEI — parallel mclapply approach ------------------------------
# Splits land cells into N_CORES chunks and processes in parallel.
# Each worker gets a copy of the data — RAM scales with N_CORES.
# If RAM exceeds 28GB, terminate and rerun with N_CORES <- 2

message("\nComputing SPEI-", SCALE, " with ", N_CORES, " cores...")
message("Monitor RAM in Activity Monitor — if > 28GB, stop and set N_CORES <- 2")
t0 <- Sys.time()

# Split land cells evenly across cores
chunks <- split(land_idx,
                cut(seq_along(land_idx), N_CORES, labels = FALSE))

# Each chunk processes its assigned land cells independently
chunk_results <- mclapply(chunks, function(idx) {

  chunk_mat <- matrix(NA_real_, nrow = length(idx), ncol = n_time)

  for (k in seq_along(idx)) {
    cell <- idx[k]
    x    <- cwb_matrix[cell, ]

    chunk_mat[k, ] <- tryCatch({
      as.numeric(SPEI::spei(
        ts(x, start = c(1980, 1), frequency = 12),
        scale     = SCALE,
        ref.start = c(1980, 1),
        ref.end   = c(2009, 12)
      )$fitted)
    }, error = function(e) rep(NA_real_, n_time))
  }

  list(idx = idx, result = chunk_mat)

}, mc.cores = N_CORES)

# Reassemble results into full matrix
message("Reassembling results...")
result <- matrix(NA_real_, nrow = n_cells, ncol = n_time)
for (chunk in chunk_results) {
  result[chunk$idx, ] <- chunk$result
}
rm(chunk_results); gc()

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
message("SPEI computation done in ", elapsed, " min.")

# ---- Validate ---------------------------------------------------------------

samp    <- result[sample(land_idx, min(10, length(land_idx))), 361]
n_valid <- sum(!is.na(samp))
rng     <- range(samp, na.rm = TRUE)
message(sprintf("Validation: %d/10 valid, range %.2f to %.2f",
                n_valid, rng[1], rng[2]))
if (n_valid == 0) stop("Validation failed — all NA. Check CWB layers.")

# ---- Write TIFs -------------------------------------------------------------

message("Writing ", n_time, " monthly TIFs...")

for (j in seq_len(n_time)) {
  date_str <- format(dates[j], "%Y%m")
  out_file <- file.path(sc_dir,
                        sprintf("SPEI%d_%s_tile%02d.tif",
                                SCALE, date_str, TILE_INDEX))
  writeRaster(setValues(r_template, result[, j]), out_file, overwrite = TRUE)
  if (j %% 60 == 0) message("  Written ", j, "/", n_time)
}

# ---- Cleanup ----------------------------------------------------------------

rm(cwb_matrix, result); gc()

message("\n========================================")
message("DONE: Tile ", TILE_INDEX, " SPEI-", SCALE)
message("Time: ", elapsed, " min")
message("Finished: ", format(Sys.time()))
message("========================================")
message("\nNext steps:")
message("  Change TILE_INDEX and/or SCALE and rerun")
message("  Crop step will be skipped for same tile")

