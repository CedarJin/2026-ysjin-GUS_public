"""Audit abundance workbook and cohort metadata on a compute node."""
from pathlib import Path
import json
import pickle
import pandas as pd

root = Path('/quobyte/angelazgrp/2026-ysjin-GUS_public')
data = root / 'software/InvestiGUT/data'
out = root / 'results/investigut_gus/input_audit'
out.mkdir(exist_ok=True)
df = pd.read_excel(data / 'Leviatan_species_abundance.xlsx', index_col=0)
numeric = df.apply(pd.to_numeric, errors='raise')
audit = dict(shape=df.shape, first_ids=list(map(str, df.index[:10])),
             last_ids=list(map(str, df.index[-10:])), first_columns=list(df.columns[:10]),
             duplicate_ids=int(df.index.duplicated().sum()),
             duplicate_columns=int(df.columns.duplicated().sum()),
             missing_cells=int(df.isna().sum().sum()),
             minimum=float(numeric.min().min()), maximum=float(numeric.max().max()),
             row_sum_quantiles=numeric.sum(axis=1).quantile([0,.25,.5,.75,1]).to_dict())
numeric.to_pickle(out / 'abundance.pkl')
with open(data / 'metagenome_files_set_d019ca77b2cac7064e174cfe73032533.pkl','rb') as handle:
    samples = pickle.load(handle)
frames=[]
for path in sorted((data / 'metadata').glob('*/*.tsv*')):
    frame=pd.read_csv(path, sep='\t', dtype=str)
    frame['metadata_file']=str(path.relative_to(data))
    frames.append(frame)
meta=pd.concat(frames, ignore_index=True)
meta=meta.loc[[tuple(x) in samples for x in meta[['study_name','sample_id']].values]].copy()
meta.to_csv(out / 'included_metadata.tsv', sep='\t', index=False)
audit['sample_cache_count']=len(samples)
audit['metadata_rows']=len(meta)
audit['metadata_columns']=list(meta.columns)
audit['duplicate_sample_keys']=int(meta.duplicated(['study_name','sample_id']).sum())
audit['subject_missing']=int(meta.subject_id.isna().sum())
audit['subjects']=len(meta[['study_name','subject_id']].drop_duplicates())
with open(out / 'audit.json','w') as handle:
    json.dump(audit,handle,indent=2)
print(json.dumps(audit,indent=2))
