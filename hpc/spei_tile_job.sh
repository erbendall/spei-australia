#!/usr/bin/env bash
# =============================================================================
# SPEI Pipeline — PBS Array Job Script (Tile Processing)
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
#              All scientific decisions and final code review remain the
#              author's responsibility.
#
# Date:        April 2026
# Version:     1.0
#
# Description:
#   PBS array job script that submits 9 independent tile jobs simultaneously.
#   Each job processes one spatial tile of Australia using spei_tile.R.
#
# HPC Compatibility:
#   This script was developed and tested on NCI Gadi (Australia) but is
#   designed to work on ANY PBS-based HPC system. To adapt for your system:
#   1. Update the #PBS directives (queue name, memory, walltime as needed)
#   2. Update MODULE_CMD to load R on your system (e.g. "module load R/4.3.1")
#   3. Update R_LIB to your personal R library path
#   4. Remove or modify the #PBS -l storage line (NCI Gadi-specific)
#   5. Update USER SETTINGS paths to your project scratch directories
#
#   Tested on:    NCI Gadi (PBS Pro, normal queue)
#   Should work:  Pawsey Setonix, MASSIVE M3, university PBS/TORQUE clusters
#
# Usage:
#   Submitted automatically by spei_submit.sh, or manually:
#   qsub -P YOUR_PROJECT spei_tile_job.sh
#
# After all tiles complete, submit the merge job:
#   qsub -W depend=afterok:<ARRAY_JOB_ID> spei_merge_job.sh
#
# Customisation:
#   Update PBS directives below to match your HPC system requirements.
#   Update USER SETTINGS section with your project paths.
#
# Citation:
#   Bendall, E.R. (2026). Gridded SPEI for Australia 1980-2024 at 0.01 degree
#   resolution derived from ANUClimate v2.0. Zenodo.
#   https://doi.org/10.5281/zenodo.XXXXXXX
#   Code: https://github.com/erbendall/spei-australia
# =============================================================================

# ---- PBS directives — update for your HPC system ----------------------------

#PBS -N spei_tile
#PBS -q normal
#PBS -l ncpus=2
#PBS -l mem=96GB
#PBS -l walltime=12:00:00
#PBS -l wd
#PBS -j oe
#PBS -J 1-9

# Storage directive — update for your HPC system and data paths
# NCI Gadi example: #PBS -l storage=gdata/gh70+gdata/dk92+scratch/YOUR_PROJECT
# Remove or modify this line for non-NCI systems
#PBS -l storage=gdata/gh70+scratch/YOUR_PROJECT

# ---- USER SETTINGS — update these ------------------------------------------

SCRIPT_DIR="${SPEI_SCRIPT_DIR:-$HOME/spei}"
OUTPUT_DIR="${SPEI_OUTPUT_DIR:-/scratch/YOUR_PROJECT/spei_outputs}"
SCRATCH_DIR="${SPEI_SCRATCH_DIR:-/scratch/YOUR_PROJECT/spei_temp}"

# R library path — update to your R library location
R_LIB="${HOME}/R/x86_64-conda-linux-gnu-library/4.3.1"

# Module load command — update for your HPC system
# NCI Gadi:    module use /g/data/dk92/apps/Modules/modulefiles && module load NCI-data-analysis/2024.01
# Pawsey:      module load r/4.3.1
# Generic:     module load R/4.3.1
MODULE_CMD="module use /g/data/dk92/apps/Modules/modulefiles && module load NCI-data-analysis/2024.01"

# -----------------------------------------------------------------------------

echo "=========================================="
echo "SPEI Tile Job"
echo "Job ID:      $PBS_JOBID"
echo "Array index: $PBS_ARRAY_INDEX of 9"
echo "Node:        $(hostname)"
echo "Start:       $(date)"
echo "=========================================="

# Load R environment
eval "${MODULE_CMD}"

# Force TMPDIR to project scratch — avoids module scratch inode limits
export TMPDIR="${SCRATCH_DIR}/tmp_$$"
mkdir -p "${TMPDIR}"

# Set R library
export R_LIBS_USER="${R_LIB}"

# Export paths for R script
export SPEI_OUTPUT_DIR="${OUTPUT_DIR}"
export SPEI_SCRATCH_DIR="${SCRATCH_DIR}"

echo "Script dir:  ${SCRIPT_DIR}"
echo "Output dir:  ${OUTPUT_DIR}"
echo "Scratch dir: ${SCRATCH_DIR}"
echo "TMPDIR:      ${TMPDIR}"
echo ""

# Run tile script
echo "Starting tile ${PBS_ARRAY_INDEX}..."
Rscript "${SCRIPT_DIR}/02_spei_tile_hpc.R" "${PBS_ARRAY_INDEX}"

# Report exit status
EXIT_CODE=$?
echo ""
echo "=========================================="
echo "Tile ${PBS_ARRAY_INDEX} finished: $(date)"
echo "Exit code: ${EXIT_CODE}"
if [ "${EXIT_CODE}" -eq 0 ]; then
    echo "Status: SUCCESS"
else
    echo "Status: FAILED — check log"
fi

# Clean up tmp directory
rm -rf "${TMPDIR}"

echo "=========================================="
exit "${EXIT_CODE}"
