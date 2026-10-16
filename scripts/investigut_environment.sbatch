#!/bin/bash
#SBATCH --job-name=gus_ig_env
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-env-%j.out
#SBATCH --error=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-env-%j.err

source /etc/profile.d/modules.sh
module load conda
set -euo pipefail
PROJECT=/quobyte/angelazgrp/2026-ysjin-GUS_public
export CONDA_PKGS_DIRS="$PROJECT/software/conda-pkgs"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$PROJECT/software/envs" "$CONDA_PKGS_DIRS"
ENV_PREFIX="$PROJECT/software/envs/investigut"
# Requeue resumes an existing environment through install instead of replacing it.
ACTION=create
[[ ! -f "$ENV_PREFIX/conda-meta/history" ]] || ACTION=install
mamba "$ACTION" -y --no-channel-priority --override-channels -c conda-forge -c bioconda \
    -p "$ENV_PREFIX" python=3.11 numpy=1.24.3 scipy=1.10.1 matplotlib-base \
    seaborn=0.13.0 libsqlite=3.48.0 pandas=1.5.3 statsmodels=0.13.5 \
    ete3=3.1.2 openpyxl=3.0.10 diamond=2.1.11
conda activate "$ENV_PREFIX"
python -c 'import numpy, scipy, matplotlib, seaborn, pandas, statsmodels, ete3, openpyxl; print("All required imports passed")'
diamond version
conda list --explicit > "$PROJECT/software/investigut-environment.explicit.txt"
touch "$PROJECT/software/ENVIRONMENT_READY"
