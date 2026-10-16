"""Merge study-scoped subjects only through documented sample duplication."""
import csv
import re
from collections import defaultdict
from pathlib import Path

root = Path('/quobyte/angelazgrp/2026-ysjin-GUS_public')
audit = root / 'results/investigut_gus/input_audit'
reference = root / 'reference/cohort_audit'
with open(audit / 'included_metadata.tsv') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
samples = {(r['study_name'], r['sample_id']): r['study_name']+'::'+r['subject_id'] for r in rows}
parent = {key: key for key in samples.values()}
edges = []


def find(key):
    while parent[key] != key:
        parent[key] = parent[parent[key]]
        key = parent[key]
    return key


def join(a, b, evidence):
    if a == b:
        return
    edges.append((a, b, evidence))
    a, b = find(a), find(b)
    parent[max(a, b)] = min(a, b)


source = (reference / 'duplicates_table_source.txt').read_text().strip()
with open(reference / 'duplicates_table.tsv') as handle:
    for number, fields in enumerate(csv.reader(handle, delimiter='\t'), 1):
        assert len(fields) == 4
        a, b = tuple(fields[:2]), tuple(fields[2:])
        if a in samples and b in samples:
            join(samples[a], samples[b], f'{source}#L{number}')
runs = defaultdict(set)
for row in rows:
    for accession in re.findall(r'[SED]RR\d+', row.get('NCBI_accession', '')):
        runs[accession].add(samples[(row['study_name'], row['sample_id'])])
for accession, subjects in runs.items():
    subjects = sorted(subjects)
    for other in subjects[1:]:
        join(subjects[0], other, f'shared_run_accession:{accession}')
# Preserve the strictly documented duplicate-sample mapping as a sensitivity case.
strict_parent = {key: find(key) for key in parent}
with open(reference / 'ena_li_nielsen_overlap.tsv') as handle:
    ena = {r['run_accession']: r for r in csv.DictReader(handle, delimiter='\t')}
by_subject = defaultdict(list)
for row in rows:
    if row['study_name'] in ('LiJ_2014', 'NielsenHB_2014'):
        by_subject[row['subject_id']].append(row)
for subject, records in by_subject.items():
    if len({r['study_name'] for r in records}) != 2:
        continue
    selected = [min((r for r in records if r['study_name'] == study), key=lambda r: r['sample_id'])
                for study in ('LiJ_2014', 'NielsenHB_2014')]
    accessions = [r['NCBI_accession'].split(';')[0] for r in selected]
    aliases = [ena[a]['sample_alias'] for a in accessions]
    assert all(re.sub(r'-\d+$', '', alias).replace('.', '_') == subject for alias in aliases)
    for field in ('country', 'gender', 'age'):
        known = {r[field] for r in selected if r.get(field)}
        assert len(known) <= 1, (subject, field, known)
    join(*(samples[(r['study_name'], r['sample_id'])] for r in selected),
         'inferred_MetaHIT_longitudinal_identity:matching_curated_subject_and_ENA_alias_stem:' + ','.join(accessions))
evidence = defaultdict(set)
for a, b, source in edges:
    evidence[find(a)].add(source)
with open(audit / 'subject_map.tsv', 'w') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['study_subject_key', 'person_id', 'evidence'])
    for key in sorted(parent):
        canonical = find(key)
        writer.writerow([key, canonical, ';'.join(sorted(evidence[canonical])) or 'study_scoped_subject_id'])
with open(audit / 'subject_map_strict.tsv', 'w') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['study_subject_key', 'person_id'])
    writer.writerows(sorted(strict_parent.items()))
with open(audit / 'subject_merge_edges.tsv', 'w') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['subject_a', 'subject_b', 'evidence'])
    writer.writerows(sorted(set(edges)))
print(f'Study-scoped subjects: {len(parent)}; merged people: {len({find(k) for k in parent})}; evidence edges: {len(edges)}')
