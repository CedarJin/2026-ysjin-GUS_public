# GUS prevalence and count distribution

Code for country-level prevalence of 11 reference gut microbial
beta-glucuronidases (GUS), and the distribution of the number of distinct
reference GUS detected per healthy individual.

The current scripts retain countries with **at least 30 people (n >= 30)**.
Prevalence is the percentage of included people with a detection in any eligible
sample. The count distribution pools healthy individuals from eligible countries
and counts distinct reference GUS (0–11) per person.

## Workflow

1. Prepare the environment and InvestiGUT protein database:
   `investigut_environment.sbatch`, `investigut_download.sbatch`, and
   `investigut_build.sbatch`.
2. Search the reference sequences with `investigut_search.sbatch` and
   `filter_gus_hits.py` (identity >= 95%, query/subject coverage >= 90%,
   E-value < 1e-5).
3. Audit metadata with `audit_investigut_inputs.py`, resolve duplicate subjects
   with `build_subject_map.py`, and generate sample/person tables with
   `analyze_gus.py` (`investigut_analysis.sbatch SEARCH_JOB_ID`).
4. Generate all-cohort country prevalence with `plot_gus_country_n30.py`.
   Run `analyze_healthy_gus.py` to select healthy individuals and generate
   healthy-cohort country prevalence, then
   `plot_healthy_gus_count_country_n30.py` for the matching count distribution.

Scripts are in `scripts/`; `.sbatch` files are the corresponding HPC launchers. 

## Requirements and inputs

Python dependencies and DIAMOND are specified in `investigut_environment.sbatch`.
Provide these inputs locally:

- `reference/active_gus.faa`: the 11 reference sequences with clade/species headers.
- `software/InvestiGUT/`: upstream software and its protein database, sample
  caches, metadata, taxonomy tables, and `Leviatan_species_abundance.xlsx`.

## Reference
Schmitz, M.A., Dimonaco, N.J., Clavel, T. et al. Lineage-specific microbial protein prediction enables large-scale exploration of protein ecology within the human gut. Nat Commun 16, 3204 (2025). https://doi.org/10.1038/s41467-025-58442-w

