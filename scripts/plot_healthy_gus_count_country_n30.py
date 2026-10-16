from pathlib import Path
import argparse,json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);o=p.parse_args().out;o.mkdir(parents=True,exist_ok=True)
r=Path('/quobyte/angelazgrp/2026-ysjin-GUS_public');a=r/'results/investigut_gus/healthy_only_24761698_r0'
qs=[l[1:].split()[0] for l in (r/'reference/active_gus.faa').read_text().splitlines() if l.startswith('>')]
d=pd.read_csv(a/'gus_subject_presence.tsv',sep='\t');assert len(qs)==11 and len(d)==4244 and d.person_id.is_unique
assert d[qs].isin([0,1]).all().all()
country_n=pd.read_csv(a/'country_n30/country_denominators.tsv',sep='\t').set_index('stratum_value').included_people
samples=pd.read_csv(a/'gus_sample_presence.tsv',sep='\t',usecols=['person_id','country']+qs)
assert samples.groupby('person_id').country.nunique().eq(1).all()
selected=samples[samples.country.isin(country_n.index)]
actual=selected.groupby('country').person_id.nunique().reindex(country_n.index)
assert actual.equals(country_n)
d=d[d.person_id.isin(selected.person_id)].copy()
assert len(d)==4218 and len(country_n)==22 and country_n.ge(30).all()
country_n.to_csv(o/'country_denominators.tsv',sep='\t')
n=d[qs].sum(axis=1);counts=n.value_counts().reindex(range(12),fill_value=0).sort_index()
assert counts.sum()==4218
# Cross-check per-person counts from original eligible healthy sample detections.
s=pd.read_csv(a/'gus_sample_presence.tsv',sep='\t',usecols=['person_id']+qs)
ind=s.groupby('person_id')[qs].max().sum(axis=1)
assert np.array_equal(ind.reindex(d.person_id).to_numpy(),n.to_numpy())
t=pd.DataFrame({'distinct_reference_GUS':range(12),'people':counts.values,'percent':counts.values/len(d)*100});t.to_csv(o/'active_gus_count_distribution.tsv',sep='\t',index=False)
d[['person_id']].assign(distinct_reference_GUS=n).to_csv(o/'individual_active_gus_counts.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(11,7));fig.subplots_adjust(left=.10,right=.98,top=.97,bottom=.14)
ax.bar(range(12),counts.values,width=.72,color='#31738F')
for x,y in counts.items():ax.text(x,y+max(counts)*.012,f'{y:,}',ha='center',va='bottom',fontsize=11)
ax.set_xticks(range(12));ax.set_xlim(-.6,11.6);ax.set_ylim(0,max(counts)*1.15);ax.set_xlabel('Number of distinct reference GUS detected',fontsize=14);ax.set_ylabel('Number of people',fontsize=14);ax.spines[['top','right']].set_visible(False)
for ext in ['png','pdf','svg']:fig.savefig(o/f'01_active_gus_count_distribution.{ext}',dpi=300,bbox_inches='tight')
(o/'validation.json').write_text(json.dumps({'status':'passed','people':len(d),'samples':len(selected),'countries':len(country_n),'minimum_country_n':30,'heatmap_country_denominators_match':True,'any_positive':int(n.gt(0).sum()),'zero_detected':int(counts.loc[0]),'sample_matrix_crosscheck':'passed','mean_detected':float(n.mean()),'median_detected':float(n.median())},indent=2))
print(t.to_string(index=False))
