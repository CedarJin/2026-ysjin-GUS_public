#!/bin/bash
#SBATCH --job-name=gus_ig_analysis
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-analysis-%j.out
#SBATCH --error=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-analysis-%j.err

source /etc/profile.d/modules.sh
module load conda
set -euo pipefail
PROJECT=/quobyte/angelazgrp/2026-ysjin-GUS_public
conda activate "$PROJECT/software/envs/investigut"
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MPLBACKEND=Agg
SEARCH_JOB=${1:?Pass the completed search job ID}
shopt -s nullglob
matches=("$PROJECT/results/investigut_gus/"id95_cov90_evalue1e-5_"${SEARCH_JOB}"_r*/SEARCH_COMPLETED)
[[ ${#matches[@]} == 1 ]] || { echo 'Expected exactly one completed search attempt' >&2; exit 1; }
RUN_DIR=$(dirname "${matches[0]}")
ATTEMPT="$RUN_DIR/postprocess_${SLURM_JOB_ID}_r${SLURM_RESTART_COUNT:-0}"
mkdir "$ATTEMPT"
cp "$PROJECT/scripts/analyze_gus.py" "$PROJECT/scripts/filter_gus_hits.py" "$ATTEMPT/"
sha256sum "$ATTEMPT/"*.py > "$ATTEMPT/code.sha256"
python "$ATTEMPT/analyze_gus.py" --run "$RUN_DIR" --root "$PROJECT" \
    --subject-map "$PROJECT/results/investigut_gus/input_audit/subject_map.tsv" \
    --out "$ATTEMPT/analysis"
touch "$ATTEMPT/TABLES_AND_FIGURES_GENERATED"
# Scientific review and independent output verification are still required.
