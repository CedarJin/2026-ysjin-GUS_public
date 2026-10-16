#!/bin/bash
#SBATCH --job-name=gus_ig_download
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-download-%j.out
#SBATCH --error=/quobyte/angelazgrp/2026-ysjin-GUS_public/logs/ig-download-%j.err

set -euo pipefail
PROJECT=/quobyte/angelazgrp/2026-ysjin-GUS_public
cd "$PROJECT/software/InvestiGUT/data"
# Sequential transfers; only this job writes the two parts. wget resumes on requeue.
records=(8398776 8404062)
hashes=(5500330faac0ffd6d725065ef26dd018 e0d2cab65d4a5df41669462b13318b8d)
for i in 0 1; do
    part="Investigut_Proteins.fa.xz_$i"
    if [[ ! -f "$part" ]] || [[ $(stat -c %s "$part") != 39394241094 ]]; then
        wget --continue --tries=20 --timeout=120 --waitretry=30 --progress=dot:giga \
            "https://zenodo.org/records/${records[$i]}/files/$part?download=1" -O "$part"
    fi
    printf '%s  %s\n' "${hashes[$i]}" "$part" | md5sum --check -
done
cat Investigut_Proteins.fa.xz_0 Investigut_Proteins.fa.xz_1 | xz --test
touch DOWNLOAD_VERIFIED
