#!/bin/bash
#SBATCH --job-name=gus_ig_search
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-search-%j.out
#SBATCH --error=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-search-%j.err

# Initial budget, not a benchmark. Single multithreaded CPU task; no GPU.
# Unique logs and result directories prevent multiple writers on Quobyte.
source /etc/profile.d/modules.sh
module load conda
set -euo pipefail

PROJECT=/quobyte/angelazgrp/2026-ysjin-GUS_public
REPO="$PROJECT/software/InvestiGUT"
conda activate "$PROJECT/software/envs/investigut"
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

MIN_ID=95
MIN_COVER=90
SENSITIVE="${GUS_SENSITIVE:-0}"
RUN_DIR="$PROJECT/results/investigut_gus/id${MIN_ID}_cov${MIN_COVER}_evalue1e-5_${SLURM_JOB_ID}_r${SLURM_RESTART_COUNT:-0}"
mkdir -p "$PROJECT/results/investigut_gus"
mkdir "$RUN_DIR"
cd "$RUN_DIR"

test -s "$REPO/data/investigut.dmnd"
test -s "$REPO/data/metagenome_files_set_d019ca77b2cac7064e174cfe73032533.pkl"
test -s "$REPO/data/mag_files_set_d019ca77b2cac7064e174cfe73032533.pkl"
cp "$PROJECT/reference/active_gus.faa" "$RUN_DIR/active_gus.faa"
printf "%s  %s\n" a12926d1b7a4739e8708cacd2e3825553bb5f17be7dfc11e28f63dc4a52cec2a active_gus.faa | sha256sum --check -
sha256sum active_gus.faa > input.sha256
diamond version > diamond-version.txt
diamond dbinfo --db "$REPO/data/investigut.dmnd" > database-info.txt
stat --printf="%n\t%s\t%y\n" "$REPO/data/investigut.dmnd" > database-stat.tsv
git -C "$REPO" rev-parse HEAD > investigut-commit.txt
git -C "$REPO" diff > investigut-local.patch
printf 'identity=%s\nquery_cover=%s\nsubject_cover=%s\nevalue_strict_lt=1e-5\nsensitive=%s\nthreads=%s\n' \
    "$MIN_ID" "$MIN_COVER" "$MIN_COVER" "$SENSITIVE" "$SLURM_CPUS_PER_TASK" > parameters.txt

# Only temporary alignment data go to local scratch; final outputs stay on group storage.
TASK_SCRATCH="${TMPDIR:-/tmp}/investigut-gus-${SLURM_JOB_ID}"
mkdir -p "$TASK_SCRATCH"
extra_args=()
if [[ "$SENSITIVE" == 1 ]]; then
    extra_args+=(--sensitive)
fi

diamond blastp \
    --db "$REPO/data/investigut.dmnd" \
    --query "$RUN_DIR/active_gus.faa" \
    --out "$RUN_DIR/active_gus_hits.raw.tsv" \
    --id "$MIN_ID" --query-cover "$MIN_COVER" --subject-cover "$MIN_COVER" \
    --evalue 1e-5 --max-target-seqs 0 \
    --threads "$SLURM_CPUS_PER_TASK" --tmpdir "$TASK_SCRATCH" \
    "${extra_args[@]}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen \
        qstart qend sstart send evalue bitscore qcovhsp scovhsp full_sseq

# Apply the strict boundary and retain zero-hit queries in the search QC.
python "$PROJECT/scripts/filter_gus_hits.py" \
    --query "$RUN_DIR/active_gus.faa" \
    --input "$RUN_DIR/active_gus_hits.raw.tsv" \
    --output "$RUN_DIR/active_gus_hits.tsv" \
    --qc "$RUN_DIR/search_qc.tsv"
sha256sum active_gus_hits.raw.tsv active_gus_hits.tsv > hits.sha256
touch "$RUN_DIR/SEARCH_COMPLETED"
# This marks search only. Subject-level prevalence and FPF are separate downstream steps.
