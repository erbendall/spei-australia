# =============================================================================
# SPEI Pipeline — HPC Tile Processing Script
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
#              optimisation strategies, PBS job scripting, error diagnosis,
#              and iterative debugging. All scientific decisions, parameter
#              choices and final code review remain the author's responsibility.
#
# Date:        April 2026
# Version:     1.0
#
# Description:
#   Processes one spatial tile of the SPEI tiled pipeline using a matrix
#   approach for memory-efficient SPEI computation on PBS-based HPC systems.
#   Called by spei_tile_job.sh with tile index passed as argument.
#
# Memory design:
#   - CWB layers cropped one at a time (~500MB peak during crop)
#   - Cropped stack read as a numeric matrix via values()
#   - SPEI applied cell-by-cell in a single R session
#   - No parallel workers — predictable, low memory usage
#   - Peak RAM per tile: ~20GB (set PBS mem to 96GB for safety)
#
# Tile grid (3 cols x 3 rows):
#   1  2  3   (lat: -9   to -20.7)
#   4  5  6   (lat: -20.7 to -32.3)
#   7  8  9   (lat: -32.3 to -44)
#
# Usage: Rscript 02_spei_tile.R <tile_index>
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

# ---- Parse tile index -------------------------------------------------------

args       <- commandArgs(trailingOnly = TRUE)
tile_index <- suppressWarnings(as.integer(args[1]))

if (length(args) == 0 || is.na(tile_index) ||
    tile_index < 1L || tile_index > 9L) {
  stop("Invalid tile index. Must be 1-9. Got: '",
       ifelse(length(args) == 0, "none", args[1]), "'")
}

message("========================================")
message("SPEI Tile ", tile_index, " of 9")
message("Start: ", format(Sys.time()))
message("========================================")

# ---- Paths ------------------------------------------------------------------

# ---- Paths — set via environment variables from PBS job script --------------
# These are passed from spei_tile_job.sh via -v flag.
# Can also be set manually for standalone use.

scratch  <- Sys.getenv("SPEI_SCRATCH_DIR",
                        unset = "/scratch/YOUR_PROJECT/spei_temp")
out_dir  <- Sys.getenv("SPEI_OUTPUT_DIR",
                        unset = "/scratch/YOUR_PROJECT/spei_outputs")

cwb_dir      <- file.path(scratch, "cwb_layers")
tile_dir     <- file.path(out_dir, "tiles")
tile_cwb_dir <- file.path(scratch, sprintf("tile%02d_cwb", tile_index))

dir.create(tile_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(tile_cwb_dir, showWarnings = FALSE, recursive = TRUE)

# Force temp files to project scratch — avoids HPC module scratch inode limits
tmpdir <- file.path(scratch, paste0("tmp_", Sys.getpid()))
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
jobfs  <- Sys.getenv("PBS_JOBFS", unset = tmpdir)
terraOptions(tempdir = jobfs, memfrac = 0.4, progress = 0)
message("Terra temp: ", jobfs)

# ---- Setup ------------------------------------------------------------------

dates     <- seq(as.Date("1980-01-01"), as.Date("2024-12-01"), by = "month")
out_idx   <- which(dates >= as.Date("2010-01-01") & dates <= as.Date("2024-12-01"))
out_dates <- dates[out_idx]
n_layers  <- length(dates)

# ---- Define tile extent -----------------------------------------------------

lon_breaks <- seq(112, 154, length.out = 4)
lat_breaks <- seq(-44,  -9, length.out = 4)

n_cols  <- 3L
col_idx <- ((tile_index - 1L) %% n_cols) + 1L
row_idx <- ((tile_index - 1L) %/% n_cols) + 1L

xmin <- lon_breaks[col_idx]
xmax <- lon_breaks[col_idx + 1L]
ymin <- lat_breaks[4L - row_idx]
ymax <- lat_breaks[4L - row_idx + 1L]

tile_ext <- ext(xmin, xmax, ymin, ymax)
message(sprintf("Tile %d extent: lon %.2f-%.2f, lat %.2f-%.2f",
                tile_index, xmin, xmax, ymin, ymax))

# ---- Verify CWB layers ------------------------------------------------------

cwb_files <- file.path(cwb_dir, sprintf("cwb_%04d.tif", seq_len(n_layers)))
n_missing <- sum(!file.exists(cwb_files))
if (n_missing > 0) {
  stop(n_missing, " CWB layer files missing from: ", cwb_dir)
}
message("CWB layers: ", n_layers, " files verified")

# ---- Crop CWB to tile — one layer at a time ---------------------------------
# Peak RAM during this step: ~500MB (one layer at a time)

tile_cwb_files <- file.path(tile_cwb_dir,
                             sprintf("cwb_%04d.tif", seq_len(n_layers)))
n_existing     <- sum(file.exists(tile_cwb_files))

if (n_existing == n_layers) {
  message("All ", n_layers, " tile CWB layers already exist — skipping crop.")
} else {
  message("Cropping CWB to tile (", n_existing, "/", n_layers, " exist)...")

  # Snap extent to grid using first layer
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

  n_written <- sum(file.exists(tile_cwb_files))
  if (n_written != n_layers) {
    stop("Expected ", n_layers, " tile layers, found ", n_written)
  }
  message("All ", n_written, " tile CWB layers cropped.")
}

# ---- Load tile as virtual stack then read to matrix -------------------------
# Virtual stack of small cropped TIFs — much smaller than full Australia stack
# values() reads entire stack into RAM as a matrix in one efficient call

message("Loading tile CWB virtual stack...")
r_cwb_tile <- rast(tile_cwb_files)
message(sprintf("Tile raster: %d rows x %d cols x %d layers",
                nrow(r_cwb_tile), ncol(r_cwb_tile), nlyr(r_cwb_tile)))

# Check for land cells
n_land <- global(!is.na(r_cwb_tile[[1]]), "sum")[1, 1]
if (n_land == 0) {
  message("Tile ", tile_index, ": no land cells — skipping.")
  writeLines(paste0("OCEAN_", tile_index),
             file.path(tile_dir, sprintf("tile%02d_ocean.txt", tile_index)))
  quit(save = "no", status = 0)
}
message("Land cells: ", format(round(n_land), big.mark = ","))

# Read entire tile into matrix — ONE read operation, no parallel overhead
# This is the key step: rows = cells, cols = 540 months
message("Reading tile into matrix...")
cwb_matrix <- values(r_cwb_tile)
rm(r_cwb_tile); gc()

message("Matrix: ", nrow(cwb_matrix), " x ", ncol(cwb_matrix),
        " (", round(object.size(cwb_matrix)/1e6, 1), " MB)")

# Template raster for output — single layer with correct extent/CRS
r_template <- rast(tile_cwb_files[1])

# Land cell indices — only process non-NA cells
land_idx <- which(!is.na(cwb_matrix[, 1]))
message("Processing ", format(length(land_idx), big.mark=","), " land cells")

# ---- SPEI computation — matrix approach -------------------------------------

compute_spei_matrix <- function(cwb_mat, land_idx, scale) {

  n_cells <- nrow(cwb_mat)
  n_time  <- ncol(cwb_mat)
  result  <- matrix(NA_real_, nrow = n_cells, ncol = n_time)

  message("  Computing SPEI-", scale, " for ",
          format(length(land_idx), big.mark=","), " cells...")

  for (k in seq_along(land_idx)) {
    cell <- land_idx[k]
    x    <- cwb_mat[cell, ]

    result[cell, ] <- tryCatch({
      as.numeric(SPEI::spei(
        ts(x, start = c(1980, 1), frequency = 12),
        scale     = scale,
        ref.start = c(1980, 1),
        ref.end   = c(2009, 12)
      )$fitted)
    }, error = function(e) rep(NA_real_, n_time))

    if (k %% 50000 == 0) {
      message("    ", format(k, big.mark=","), "/",
              format(length(land_idx), big.mark=","),
              " (", round(k/length(land_idx)*100), "%)")
    }
  }
  result
}

# ---- Write outputs ----------------------------------------------------------

write_tile_outputs <- function(result_mat, scale, tile_index) {

  # NetCDF — full timeseries 1980-2024
  nc_file <- file.path(tile_dir,
                       sprintf("SPEI%d_1980_2024_tile%02d.nc",
                               scale, tile_index))

  message("  Building output raster stack...")
  r_out <- rast(lapply(seq_len(ncol(result_mat)), function(j) {
    setValues(r_template, result_mat[, j])
  }))

  writeCDF(r_out, nc_file,
           varname  = paste0("SPEI", scale),
           longname = paste0("SPEI ", scale, "-month, ref 1980-2009"),
           unit     = "unitless", overwrite = TRUE)
  message("  NC: ", round(file.size(nc_file)/1e6, 1), " MB")

  # Monthly TIFs — 2010-2024
  sc_dir <- file.path(tile_dir, paste0("SPEI", scale))
  dir.create(sc_dir, showWarnings = FALSE, recursive = TRUE)
  for (i in seq_along(out_idx)) {
    writeRaster(
      setValues(r_template, result_mat[, out_idx[i]]),
      file.path(sc_dir, sprintf("SPEI%d_%s_tile%02d.tif",
                                scale, format(out_dates[i], "%Y%m"),
                                tile_index)),
      overwrite = TRUE
    )
  }
  message("  TIFs: ", length(out_idx), " files")
  rm(r_out); gc()
}

validate_result <- function(result_mat, land_idx, scale, tile_index) {
  samp    <- result_mat[sample(land_idx, min(10, length(land_idx))), out_idx[1]]
  n_valid <- sum(!is.na(samp))
  if (n_valid == 0) {
    warning("Tile ", tile_index, " SPEI-", scale, ": all NA")
    return(FALSE)
  }
  rng <- range(samp, na.rm = TRUE)
  message(sprintf("  Validation: %d/10 valid, range %.2f to %.2f",
                  n_valid, rng[1], rng[2]))
  TRUE
}

# ---- Run SPEI for all scales ------------------------------------------------

for (sc in c(3, 6, 12, 24)) {

  nc_file   <- file.path(tile_dir,
                         sprintf("SPEI%d_1980_2024_tile%02d.nc",
                                 sc, tile_index))
  sc_dir    <- file.path(tile_dir, paste0("SPEI", sc))
  tif_count <- length(list.files(sc_dir,
                                 pattern = sprintf("_tile%02d\\.tif$",
                                                   tile_index)))

  if (file.exists(nc_file) && tif_count == length(out_idx)) {
    message("SPEI-", sc, " tile ", tile_index, " done — skipping.")
    next
  }

  message("\n--- SPEI-", sc, " tile ", tile_index, " ---")
  t0         <- Sys.time()
  result_mat <- compute_spei_matrix(cwb_matrix, land_idx, sc)

  validate_result(result_mat, land_idx, sc, tile_index)
  write_tile_outputs(result_mat, sc, tile_index)

  rm(result_mat); gc()
  message(sprintf("  Done in %.1f min.",
                  as.numeric(difftime(Sys.time(), t0, units="mins"))))
}

# ---- Cleanup ----------------------------------------------------------------

rm(cwb_matrix); gc()
unlink(tile_cwb_dir, recursive = TRUE)
message("Tile CWB temp files cleaned up.")

writeLines(paste0("COMPLETE_", tile_index, "_", format(Sys.time())),
           file.path(tile_dir, sprintf("tile%02d_complete.txt", tile_index)))

message("\n========================================")
message("Tile ", tile_index, " complete: ", format(Sys.time()))
message("========================================")
