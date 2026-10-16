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
`plot_gus_characteristics.py` and the label helpers are retained because the
healthy-cohort analysis calls them. Existing upstream scripts also produce
additional summaries as part of their processing.

## Requirements and inputs

Python dependencies and DIAMOND are specified in `investigut_environment.sbatch`.
Provide these inputs locally (they are not included in this code-only repository):

- `reference/active_gus.faa`: the 11 reference sequences with clade/species headers.
- `software/InvestiGUT/`: upstream software and its protein database, sample
  caches, metadata, taxonomy tables, and `Leviatan_species_abundance.xlsx`.
- `reference/cohort_audit/`: `duplicates_table.tsv`,
  `duplicates_table_source.txt`, and `ena_li_nielsen_overlap.tsv` for subject mapping.

The scripts preserve the original analysis paths, job IDs, and cohort-size
assertions. Update them for a new location or run, and create `logs/` and
`results/investigut_gus/` before starting. Run analysis on a compute node;
the supplied UC Davis Hive launchers use the preemptible `low` partition.

Only the selected workflow code and this README are intended for Git.
Inputs, generated figures/tables, logs, installed software, and unrelated scripts
are excluded by the `.gitignore` allowlist.
