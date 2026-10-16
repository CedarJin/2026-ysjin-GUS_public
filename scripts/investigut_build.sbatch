#!/bin/bash
#SBATCH --job-name=gus_ig_build
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-build-%j.out
#SBATCH --error=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-build-%j.err

source /etc/profile.d/modules.sh
module load conda
set -euo pipefail
PROJECT=/quobyte/angelazgrp/2026-ysjin-GUS_public
export PYTHONNOUSERSITE=1
conda activate "$PROJECT/software/envs/investigut"
test -f "$PROJECT/software/ENVIRONMENT_READY"
cd "$PROJECT/software/InvestiGUT/data"
test -f DOWNLOAD_VERIFIED
if [[ -f investigut.dmnd ]]; then
    diamond dbinfo --db investigut.dmnd
    touch DATABASE_READY
    exit 0
fi
# Each restart writes a fresh intermediate. Only a successful build is promoted.
BUILD_DB="investigut-building-${SLURM_JOB_ID}-${SLURM_RESTART_COUNT:-0}"
cat Investigut_Proteins.fa.xz_0 Investigut_Proteins.fa.xz_1 | \
    xz --decompress --stdout | \
    diamond makedb --in - --db "$BUILD_DB" --threads "$SLURM_CPUS_PER_TASK"
diamond dbinfo --db "$BUILD_DB.dmnd"
mv "$BUILD_DB.dmnd" investigut.dmnd
touch DATABASE_READY
