#!/usr/bin/env bash
# =============================================================================
# SPEI Pipeline — Master PBS Submission Script
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
#              All scientific decisions and final code review remain the
#              author's responsibility.
#
# Date:        April 2026
# Version:     1.0
#
# Description:
#   Master submission script for the SPEI tiled pipeline on any PBS-based
#   HPC system. Submits the tile array job and merge job with correct
#   dependency in a single command.
#
# HPC Compatibility:
#   This script was developed and tested on NCI Gadi (Australia) but is
#   designed to work on ANY PBS-based HPC system. To adapt for your system:
#   1. Update PROJECT and path variables in USER SETTINGS
#   2. Update QUEUE to match your system queue name
#   3. Update MODULE_CMD to load R on your system
#   4. Remove or modify the STORAGE variable (NCI Gadi-specific)
#
#   Tested on:    NCI Gadi (PBS Pro, normal queue)
#   Should work:  Pawsey Setonix, MASSIVE M3, university PBS/TORQUE clusters
#
# Usage:
#   bash spei_submit.sh
#
# Requirements:
#   - PBS job scheduler (qsub)
#   - spei_tile_job.sh and spei_merge_job.sh in the same directory
#   - R environment with terra, SPEI, lubridate packages
#   - ANUClimate v2.0 data accessible at DATA_PATH
#
# Before running:
#   1. Update all USER SETTINGS below
#   2. Ensure spei_tile.R and spei_merge_hpc.R are in SCRIPT_DIR
#   3. Run: bash spei_submit.sh
# =============================================================================

# ---- USER SETTINGS — update these before running ---------------------------

PROJECT="YOUR_PROJECT_CODE"          # HPC project/allocation code
SCRIPT_DIR="$HOME/spei"              # directory containing all R and sh scripts
OUTPUT_DIR="/scratch/${PROJECT}/spei_outputs"  # where outputs are written
SCRATCH_DIR="/scratch/${PROJECT}/spei_temp"    # scratch space for temp files
DATA_PATH="/g/data/gh70/ANUClimate/v2-0/stable/month"  # ANUClimate data path

# PBS settings
QUEUE="normal"          # PBS queue name
TILE_NCPUS=2            # CPUs per tile job
TILE_MEM="96GB"         # memory per tile job
TILE_WALLTIME="12:00:00" # walltime per tile job
MERGE_NCPUS=4           # CPUs for merge job
MERGE_MEM="64GB"        # memory for merge job
MERGE_WALLTIME="04:00:00" # walltime for merge job

# Storage flags (PBS -l storage) — adjust for your HPC system
# Format: filesystem/project+filesystem/project
# Example for NCI Gadi: gdata/gh70+gdata/dk92+scratch/YOUR_PROJECT
STORAGE="scratch/${PROJECT}"

# -----------------------------------------------------------------------------

echo "============================================"
echo "SPEI Pipeline Submission"
echo "Project:    ${PROJECT}"
echo "Scripts:    ${SCRIPT_DIR}"
echo "Outputs:    ${OUTPUT_DIR}"
echo "Data:       ${DATA_PATH}"
echo "============================================"

# Verify scripts exist
for f in 02_spei_tile_hpc.R spei_tile_job.sh 02_spei_merge_hpc.R spei_merge_job.sh; do
    if [ ! -f "${SCRIPT_DIR}/${f}" ]; then
        echo "ERROR: ${SCRIPT_DIR}/${f} not found"
        exit 1
    fi
done
echo "All scripts found."

# Create output directories
mkdir -p "${OUTPUT_DIR}/tiles"
mkdir -p "${SCRATCH_DIR}/cwb_layers"
echo "Output directories created."

# Export settings for job scripts to use
export SPEI_PROJECT="${PROJECT}"
export SPEI_SCRIPT_DIR="${SCRIPT_DIR}"
export SPEI_OUTPUT_DIR="${OUTPUT_DIR}"
export SPEI_SCRATCH_DIR="${SCRATCH_DIR}"
export SPEI_DATA_PATH="${DATA_PATH}"

# Submit tile array job
echo ""
echo "Submitting tile array job (9 tiles)..."
ARRAY_ID=$(qsub \
    -P "${PROJECT}" \
    -q "${QUEUE}" \
    -l ncpus="${TILE_NCPUS}" \
    -l mem="${TILE_MEM}" \
    -l walltime="${TILE_WALLTIME}" \
    -l storage="${STORAGE}" \
    -l wd \
    -j oe \
    -J 1-9 \
    -o "${SCRIPT_DIR}/" \
    -v SPEI_PROJECT,SPEI_SCRIPT_DIR,SPEI_OUTPUT_DIR,SPEI_SCRATCH_DIR,SPEI_DATA_PATH \
    "${SCRIPT_DIR}/spei_tile_job.sh")

if [ -z "${ARRAY_ID}" ]; then
    echo "ERROR: Tile array job submission failed"
    exit 1
fi
echo "Tile array job submitted: ${ARRAY_ID}"

# Submit merge job with dependency on all tiles completing
echo ""
echo "Submitting merge job (depends on all tiles)..."
MERGE_ID=$(qsub \
    -P "${PROJECT}" \
    -q "${QUEUE}" \
    -l ncpus="${MERGE_NCPUS}" \
    -l mem="${MERGE_MEM}" \
    -l walltime="${MERGE_WALLTIME}" \
    -l storage="${STORAGE}" \
    -l wd \
    -j oe \
    -o "${SCRIPT_DIR}/spei_merge.log" \
    -W depend=afterok:"${ARRAY_ID}" \
    -v SPEI_PROJECT,SPEI_SCRIPT_DIR,SPEI_OUTPUT_DIR,SPEI_SCRATCH_DIR \
    "${SCRIPT_DIR}/spei_merge_job.sh")

if [ -z "${MERGE_ID}" ]; then
    echo "ERROR: Merge job submission failed"
    exit 1
fi
echo "Merge job submitted: ${MERGE_ID}"

echo ""
echo "============================================"
echo "Submission complete!"
echo "Tile array job: ${ARRAY_ID}"
echo "Merge job:      ${MERGE_ID}"
echo ""
echo "Monitor with:  qstat -u \$USER"
echo "Tile logs:     ${SCRIPT_DIR}/spei_tile.o*.log"
echo "Merge log:     ${SCRIPT_DIR}/spei_merge.log"
echo "============================================"
