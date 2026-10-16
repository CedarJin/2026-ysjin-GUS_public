#!/bin/bash
#SBATCH --job-name=gus_country_n30
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-country-n30-%j.out
#SBATCH --error=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-country-n30-%j.err
source /etc/profile.d/modules.sh
module load conda
set -euo pipefail
PROJECT=/quobyte/angelazgrp/2026-ysjin-GUS_public
conda activate "$PROJECT/software/envs/investigut"
export PYTHONNOUSERSITE=1
export PYTHONPATH="$PROJECT/scripts"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
OUT="$PROJECT/results/investigut_gus/country_n30_figures_${SLURM_JOB_ID}_r${SLURM_RESTART_COUNT:-0}"
mkdir "$OUT"
cp "$PROJECT/scripts/plot_gus_country_n30.py" "$OUT/"
python "$OUT/plot_gus_country_n30.py" --out "$OUT"
