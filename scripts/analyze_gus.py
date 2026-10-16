"""Joint GUS analysis; run on a compute node after SEARCH_COMPLETED.

FPF uses the published sparse workbook convention, with a separate strict
missing-value sensitivity output. Subject merges require an explicit evidence map.
"""
import argparse
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from filter_gus_hits import qualifies


def save(frame, path):
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def query_groups(fasta):
    clades = {}
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            clades[line[1:].split()[0]] = re.search(r"clade=(\S+)", line)[1]
    groups = {q: [q] for q in clades}
    groups["any_GUS"] = list(clades)
    for clade in dict.fromkeys(clades.values()):
        groups[f"clade_{clade}"] = [q for q, c in clades.items() if c == clade]
    return clades, groups


def choose_sample(rows):
    """Choose without consulting hits; use time only if complete within subject."""
    for column in ("days_from_first_collection", "visit_number", "infant_age", "age"):
        if column in rows:
            time = pd.to_numeric(rows[column], errors="coerce")
            if time.notna().all() and time.nunique() > 1:
                chosen = rows.assign(_time=time).sort_values(["_time", "study_name", "sample_id"]).iloc[0]
                return chosen.name, f"minimum_{column}"
    chosen = rows.sort_values(["study_name", "sample_id"]).iloc[0]
    return chosen.name, "single_sample" if len(rows) == 1 else "lexical_fallback_time_unresolved"


def prevalence_table(matrix, groups, rule):
    return pd.DataFrame([dict(query_id=q, positive_people=int(matrix[q].sum()),
                              included_people=len(matrix),
                              prevalence_percent=100 * matrix[q].mean(), rule=rule)
                         for q in groups])


def analyze(run, root, subject_map=None, output=None):
    if not (run / "SEARCH_COMPLETED").exists():
        raise ValueError("Search has not completed")
    out = output or run / "analysis"
    out.mkdir(exist_ok=False)
    data = root / "software/InvestiGUT/data"
    audit = root / "results/investigut_gus/input_audit"
    clades, groups = query_groups(run / "active_gus.faa")
    metadata = pd.read_csv(audit / "included_metadata.tsv", sep="\t", dtype=str)
    assert not metadata.duplicated(["study_name", "sample_id"]).any()
    assert metadata.subject_id.notna().all()
    metadata["study_subject_key"] = metadata.study_name + "::" + metadata.subject_id
    metadata["person_id"] = metadata.study_subject_key
    if subject_map:
        mapping = pd.read_csv(subject_map, sep="\t", dtype=str)
        assert not mapping.study_subject_key.duplicated().any()
        assert mapping.evidence.notna().all()
        m = mapping.set_index("study_subject_key").person_id
        metadata["person_id"] = metadata.study_subject_key.map(m).fillna(metadata.study_subject_key)
        save(mapping, out / "subject_merge_evidence.tsv")
    sample_index = {tuple(v): i for i, v in enumerate(metadata[["study_name", "sample_id"]].values)}
    sample_presence = pd.DataFrame(0, index=metadata.index, columns=list(clades), dtype=int)
    genomes = {q: set() for q in clades}
    positive_rows, source_counts = [], defaultdict(int)
    with open(run / "active_gus_hits.tsv") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            assert qualifies(fields)
            q, sid = fields[:2]
            assert q in clades
            parts = sid.split("|")
            assert len(parts) == 5, sid
            study, sample = parts[1].split("__")
            source_counts[study] += 1
            if study == "LeviatanS_2022":
                genomes[q].add(sample)
                positive_rows.append(dict(query_id=q, clade=clades[q], Rep_ID=sample,
                                          protein_id=sid, identity=float(fields[2]),
                                          query_coverage=float(fields[12]), subject_coverage=float(fields[13]),
                                          evalue=float(fields[10]), bitscore=float(fields[11])))
            elif (study, sample) in sample_index:
                sample_presence.loc[sample_index[(study, sample)], q] = 1
            else:
                raise ValueError(f"Hit source absent from inclusion list: {study} {sample}")
    for name, members in groups.items():
        sample_presence[name] = sample_presence[members].max(axis=1)
    sample = pd.concat([metadata, sample_presence], axis=1)
    save(sample, out / "gus_sample_presence.tsv")
    people = sample.groupby("person_id", sort=True)[list(groups)].max()
    person_info = sample.groupby('person_id').agg(
        n_samples=('sample_id','size'),
        study_subject_keys=('study_subject_key',lambda x: ';'.join(sorted(set(x)))),
        studies=('study_name',lambda x: ';'.join(sorted(set(x)))))
    person_info['aggregation_rule']='any_eligible_sample_positive'
    save(person_info.join(people).reset_index(), out / "gus_subject_presence.tsv")
    save(prevalence_table(people, groups, "any_eligible_sample_positive"), out / "gus_prevalence.tsv")
    strict_map_path = audit / 'subject_map_strict.tsv'
    if strict_map_path.exists():
        strict_map = pd.read_csv(strict_map_path, sep='\t').set_index('study_subject_key').person_id
        strict_keys = sample.study_subject_key.map(strict_map)
        assert strict_keys.notna().all()
        strict_people = sample.groupby(strict_keys)[list(groups)].max()
        save(prevalence_table(strict_people, groups, 'documented_duplicate_samples_only'),
             out / 'gus_prevalence_identity_sensitivity.tsv')
    selections = [(*choose_sample(rows), person) for person, rows in metadata.groupby("person_id")]
    selected = pd.DataFrame(selections, columns=["row_index", "selection_rule", "person_id"])
    baseline = sample.loc[selected.row_index].copy()
    baseline["selection_rule"] = selected.selection_rule.values
    save(baseline, out / "baseline_sample_selection.tsv")
    save(prevalence_table(baseline, groups, "one_prespecified_sample"), out / "gus_prevalence_single_sample.tsv")
    strata = []
    for column in ('study_name', 'country'):
        for value, rows in sample.groupby(column, dropna=False):
            table = prevalence_table(rows.groupby('person_id')[list(groups)].max(), groups,
                                     'any_eligible_sample_within_stratum')
            table['stratum'] = column
            table['stratum_value'] = value
            strata.append(table)
    for column in ("study_name", "country", "disease"):
        for value, rows in baseline.groupby(column, dropna=False):
            table = prevalence_table(rows, groups, "one_prespecified_sample")
            table["stratum"] = column
            table["stratum_value"] = value
            strata.append(table)
    save(pd.concat(strata), out / "gus_prevalence_stratified.tsv")
    sample_check = pd.DataFrame([dict(query_id=q, positive_samples=int(sample[q].sum()),
                                     included_samples=len(sample),
                                     sample_detection_percent=100*sample[q].mean()) for q in groups])
    save(sample_check, out / 'gus_sample_detection_crosscheck.tsv')
    candidates = metadata.groupby("subject_id").study_name.nunique()
    save(metadata[metadata.subject_id.isin(candidates[candidates > 1].index)], out / "cross_study_id_candidates.tsv")
    save(pd.DataFrame([dict(included_samples=len(sample), study_scoped_subjects=metadata.study_subject_key.nunique(),
                            included_people=len(people), repeated_people=int((sample.groupby('person_id').size()>1).sum()),
                            excluded_samples=0, exclusion_reason='none_all_cached_samples_have_subject_metadata',
                            merged_subject_keys=metadata.study_subject_key.nunique()-len(people),
                            unresolved_time_selections=int(baseline.selection_rule.str.startswith('lexical').sum()))]),
         out / "prevalence_cohort_qc.tsv")

    abundance = pd.read_pickle(audit / "abundance.pkl")
    assert not abundance.index.duplicated().any() and not abundance.columns.duplicated().any()
    assert np.isfinite(abundance.fillna(0).values).all()
    assert (abundance.fillna(0).values >= 0).all()
    assert np.allclose(abundance.sum(axis=1), 1, atol=1e-10)
    # Specific to this audited workbook: blank cells follow the author's sparse
    # absence convention. Retain strict-NA sensitivity and the blank cell counts.
    cohort = pd.Series(["IL" if str(i).startswith("IL") else "NLD" if str(i).startswith("NLD") else "UNKNOWN"
                        for i in abundance.index], index=abundance.index, name="cohort")
    assert not cohort.eq("UNKNOWN").any()
    expanded = {name: set().union(*(genomes[q] for q in members)) for name, members in groups.items()}
    tax = pd.concat([pd.read_csv(data / "metadata" / name, sep="\t", usecols=["user_genome", "classification"])
                     for name in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv")])
    positive = pd.DataFrame(positive_rows, columns=["query_id", "clade", "Rep_ID", "protein_id", "identity", "query_coverage", "subject_coverage", "evalue", "bitscore"])
    positive = positive.drop_duplicates().merge(tax, how="left", left_on="Rep_ID", right_on="user_genome", validate="many_to_one")
    save(positive, out / "gus_positive_genomes.tsv")
    membership = pd.DataFrame({q: [int(g in expanded[q]) for g in abundance.columns] for q in groups}, index=abundance.columns)
    membership.index.name = "Rep_ID"
    save(membership.reset_index(), out / "gus_genome_membership.tsv")
    fractions = pd.DataFrame(index=abundance.index)
    long, strict, checks = [], [], []
    for q, members in expanded.items():
        missing_cols = sorted(members - set(abundance.columns))
        cols = sorted(members & set(abundance.columns))
        status = "missing_genome_column" if missing_cols else "no_positive_reference_genome" if not members else "sparse_blank_as_absent"
        values = abundance[cols].sum(axis=1)
        strict_values = abundance[cols].sum(axis=1, skipna=False)
        if missing_cols:
            values[:] = np.nan
            strict_values[:] = np.nan
        fractions[q] = values
        blank_counts = abundance[cols].isna().sum(axis=1)
        long.append(pd.DataFrame(dict(individual_id=abundance.index, cohort=cohort.values, query_id=q,
                                     fpf_fraction=values.values, fpf_percent=100*values.values,
                                     blank_positive_genome_cells=blank_counts.values, status=status)))
        strict.append(pd.DataFrame(dict(individual_id=abundance.index, query_id=q,
                                       fpf_fraction_if_blanks_unknown=strict_values.values)))
        # Independent scalar check on a deterministic small subset; repeats cannot
        # change the set of genome columns used for the sum.
        if not missing_cols:
            for i in [0, len(abundance)//2, len(abundance)-1]:
                scalar = sum(float(abundance.iloc[i][g]) for g in cols if pd.notna(abundance.iloc[i][g]))
                assert np.isclose(scalar, values.iloc[i])
        checks.append(dict(query_id=q, positive_reference_genomes=len(members), missing_columns=';'.join(missing_cols),
                           max_difference_from_upstream=float((values-abundance[cols].sum(axis=1)).abs().max()) if not missing_cols else np.nan))
    long = pd.concat(long, ignore_index=True)
    save(long, out / "gus_fpf_individual.tsv")
    save(pd.concat(strict), out / "gus_fpf_missingness_sensitivity.tsv")
    wide = fractions.mul(100).rename_axis("individual_id").reset_index()
    save(wide[["individual_id"]+list(clades)], out / "gus_fpf_matrix.tsv")
    save(wide[["individual_id"]+[q for q in groups if q not in clades]], out / "gus_fpf_union.tsv")
    summaries = []
    for name, frame in [("ALL", long), *(long.groupby("cohort"))]:
        for q, rows in frame.groupby("query_id", sort=False):
            x = rows.fpf_percent
            summaries.append(dict(cohort=name, query_id=q, included=len(rows), missing=int(x.isna().sum()),
                                  zero=int(x.eq(0).sum()), mean=x.mean(), median=x.median(),
                                  q25=x.quantile(.25), q75=x.quantile(.75)))
    save(pd.DataFrame(summaries), out / "gus_fpf_summary.tsv")
    save(pd.DataFrame(checks), out / "fpf_upstream_crosscheck.tsv")
    save(pd.DataFrame([dict(individuals=len(abundance), genomes=len(abundance.columns),
                            IL=int(cohort.eq('IL').sum()), NLD=int(cohort.eq('NLD').sum()),
                            blank_cells=int(abundance.isna().sum().sum()),
                            unit='fraction', blank_policy='author_sparse_absence_with_strict_NA_sensitivity',
                            min_row_sum=abundance.sum(axis=1).min(), max_row_sum=abundance.sum(axis=1).max())]),
         out / "abundance_matrix_qc.tsv")
    assert (people.any_GUS >= people[list(clades)].max(axis=1)).all()
    assert (fractions.any_GUS.dropna() <= 1+1e-10).all()
    assert (fractions.any_GUS + 1e-10 >= fractions[list(clades)].max(axis=1)).all()
    assert (fractions.filter(like='clade_').le(fractions.any_GUS+1e-10, axis=0)).all().all()
    manifest=[]
    manifest_paths=[run/'active_gus.faa', run/'active_gus_hits.tsv', audit/'included_metadata.tsv',
                    audit/'abundance.pkl', data/'Leviatan_species_abundance.xlsx',
                    data/'metadata/gtdbtk.bac120.summary.tsv', data/'metadata/gtdbtk.ar53.summary.tsv']
    if subject_map:manifest_paths.append(subject_map)
    if strict_map_path.exists():manifest_paths.append(strict_map_path)
    for path in manifest_paths:
        digest=hashlib.sha256()
        with open(path,'rb') as handle:
            for chunk in iter(lambda: handle.read(1024*1024), b''): digest.update(chunk)
        manifest.append(dict(path=str(path), sha256=digest.hexdigest(),
                             identity_min=95, query_coverage_min=90, subject_coverage_min=90,
                             evalue_strict_lt=1e-5,
                             investigut_commit=(run/'investigut-commit.txt').read_text().strip() if (run/'investigut-commit.txt').exists() else 'test_fixture',
                             diamond_version=(run/'diamond-version.txt').read_text().strip() if (run/'diamond-version.txt').exists() else 'test_fixture'))
    save(pd.DataFrame(manifest), out/'input_manifest.tsv')
    make_figures(out, groups, clades)
    (out/'analysis_checks.json').write_text(json.dumps(dict(thresholds='95/90/90/<1e-5',
        query_count=len(clades), group_count=len(groups), sources=source_counts,
        tables_computed=True, final_review_required=True), indent=2))


def make_figures(out, groups, clades):
    from gus_plot_labels import gus_label, FPF_NOTE
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figures = out / 'figures'
    figures.mkdir(exist_ok=True)
    prevalence = pd.read_csv(out/'gus_prevalence.tsv', sep='\t')
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.barh([gus_label(q) for q in prevalence.query_id], prevalence.prevalence_percent)
    for i, row in prevalence.iterrows():
        ax.text(row.prevalence_percent+.5, i, f'{row.positive_people}/{row.included_people}', va='center', fontsize=8)
    ax.set(xlabel='People positive (%)', title='GUS prevalence: any eligible sample positive', xlim=(0,105))
    fig.tight_layout(); fig.savefig(figures/'prevalence.png', dpi=180); plt.close(fig)
    stratified=pd.read_csv(out/'gus_prevalence_stratified.tsv',sep='\t')
    for field in ('study_name','country'):
        selected=stratified[(stratified.stratum==field) & (stratified.query_id=='any_GUS') &
                            (stratified.rule=='any_eligible_sample_within_stratum')].copy()
        selected=selected.sort_values('prevalence_percent').reset_index(drop=True)
        selected['stratum_value']=selected.stratum_value.fillna('Unknown')
        fig,ax=plt.subplots(figsize=(11,max(5,len(selected)*.28)))
        ax.barh(selected.stratum_value,selected.prevalence_percent)
        for i,row in selected.iterrows():
            ax.text(row.prevalence_percent+.5,i,f'{row.positive_people}/{row.included_people}',va='center',fontsize=7)
        ax.set(xlabel='People positive (%)',title=f'Any reference GUS by {field}',xlim=(0,115))
        fig.tight_layout();fig.savefig(figures/f'prevalence_any_by_{field}.png',dpi=180);plt.close(fig)
    fpf = pd.read_csv(out/'gus_fpf_individual.tsv', sep='\t')
    fig, axes = plt.subplots(6, 3, figsize=(19, 24))
    for ax, q in zip(axes.flat, groups):
        for cohort, rows in fpf[fpf.query_id.eq(q)].groupby('cohort'):
            x = np.sort(rows.fpf_percent.dropna())
            ax.step(x, np.arange(1,len(x)+1)/len(x), where='post', label=f'{cohort} n={len(x)}')
        ax.set(title=gus_label(q,fpf=True), xlabel='FPF (%)', ylabel='Cumulative proportion', xlim=(0,100), ylim=(0,1))
        ax.legend(fontsize=7)
    fig.tight_layout(rect=(0,.025,1,1)); fig.text(.02,.008,FPF_NOTE,fontsize=10); fig.savefig(figures/'fpf_ecdf.png', dpi=150); plt.close(fig)
    matrix = pd.read_csv(out/'gus_fpf_matrix.tsv', sep='\t', index_col=0)
    fig, ax = plt.subplots(figsize=(15,12))
    im = ax.imshow(matrix.values, aspect='auto', vmin=0, vmax=100, interpolation='nearest')
    ax.set_xticks(range(len(clades)), [gus_label(q,fpf=True) for q in clades], rotation=65, ha='right')
    ax.set(ylabel='Individuals in original workbook order', title='Individual FPF; all individuals retained')
    fig.colorbar(im, ax=ax, label='FPF (%)'); fig.tight_layout(rect=(0,.05,1,1)); fig.text(.02,.015,FPF_NOTE,fontsize=9)
    fig.savefig(figures/'fpf_individual_heatmap.png', dpi=180); plt.close(fig)
    membership = pd.read_csv(out/'gus_genome_membership.tsv', sep='\t')[list(clades)]
    overlap = membership.T.dot(membership)
    fig, ax = plt.subplots(figsize=(16,13))
    im = ax.imshow(overlap.values)
    ax.set_xticks(range(len(clades)), [gus_label(q) for q in clades], rotation=65, ha='right')
    ax.set_yticks(range(len(clades)), [gus_label(q) for q in clades])
    ax.set_title('Shared positive representative genomes')
    fig.colorbar(im, ax=ax, label='Genomes'); fig.tight_layout()
    fig.savefig(figures/'gus_genome_overlap.png', dpi=180); plt.close(fig)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--root',type=Path,default=Path('/quobyte/angelazgrp/2026-ysjin-GUS_public'))
    parser.add_argument('--subject-map',type=Path)
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    analyze(args.run,args.root,args.subject_map,args.out)
