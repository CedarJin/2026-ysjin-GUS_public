#!/bin/bash
#SBATCH --job-name=gus_ig_audit
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-audit-%j.out
#SBATCH --error=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-audit-%j.err

source /etc/profile.d/modules.sh
module load conda
set -euo pipefail
PROJECT=/quobyte/angelazgrp/2026-ysjin-GUS_public
conda activate "$PROJECT/software/envs/investigut"
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
python "$PROJECT/scripts/audit_investigut_inputs.py"
