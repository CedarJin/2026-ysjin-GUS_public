"""Healthy-only individual prevalence, metadata strata and plots."""
from pathlib import Path
import argparse,json,subprocess,sys,hashlib
import numpy as np
import pandas as pd
from analyze_gus import choose_sample,query_groups,prevalence_table
from gus_plot_labels import gus_label
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--policy',choices=['strict','sample'],default='strict');args=p.parse_args();out=args.out;out.mkdir(exist_ok=True,parents=True)
root=Path('/quobyte/angelazgrp/2026-ysjin-GUS_public');src=root/'results/investigut_gus/id95_cov90_evalue1e-5_24754089_r0/postprocess_24754212_r0/analysis'
s=pd.read_csv(src/'gus_sample_presence.tsv',sep='\t',low_memory=False);clades,groups=query_groups(root/'reference/active_gus.faa');qs=list(groups)
h=s.disease.eq('healthy')
# Explicit conflicting diagnoses or case status exclude a person from the strict population.
case_condition=s.study_condition.notna()&~s.study_condition.isin(['control','cephalosporins'])
known_disease=s.disease.notna()&~s.disease.eq('healthy')
case_subtype=s.disease_subtype.notna()&~s.disease_subtype.isin(['healthy','T1D_nonconverter','T1D_seroconverter'])
bad=known_disease|case_condition|case_subtype
badpeople=set(s.loc[bad,'person_id'])
sample_eligible=h&~case_condition&~case_subtype
eligible=sample_eligible & (~s.person_id.isin(badpeople) if args.policy=='strict' else True)
good=s.loc[eligible].copy();assert len(good)>0
person=good.groupby('person_id',sort=True)[qs].max();assert person.index.is_unique
selection=[]
for pid,g in good.groupby('person_id',sort=True):
 idx,rule=choose_sample(g);selection.append((idx,rule))
b=good.loc[[i for i,_ in selection]].copy();b['selection_rule']=[r for _,r in selection];assert b.person_id.nunique()==len(person)==len(b)
# Retain consistent input names for reusable plotting scripts.
good.to_csv(out/'gus_sample_presence.tsv',sep='\t',index=False);b.to_csv(out/'baseline_sample_selection.tsv',sep='\t',index=False);person.to_csv(out/'gus_subject_presence.tsv',sep='\t')
prev=prevalence_table(person,groups,'any_eligible_healthy_sample_positive');prev.to_csv(out/'gus_prevalence.tsv',sep='\t',index=False)
prevalence_table(b,groups,'one_selected_healthy_sample').to_csv(out/'gus_prevalence_single_sample.tsv',sep='\t',index=False)
strata=[]
for field in ['country','study_name']:
 for level,g in good.groupby(field,dropna=False):
  x=prevalence_table(g.groupby('person_id')[qs].max(),groups,'any_eligible_sample_within_stratum');x['stratum']=field;x['stratum_value']=level;strata.append(x)
st=pd.concat(strata);st.to_csv(out/'gus_prevalence_stratified.tsv',sep='\t',index=False)
# Explicit per-sample reasons; preserve existing identity mapping.
audit=s[['person_id','study_name','sample_id','disease','study_condition','disease_subtype']].copy();audit['explicit_healthy']=h;audit['conflicting_case_record']=bad;audit['person_has_case_record']=s.person_id.isin(badpeople);audit['included_sample']=eligible;audit.to_csv(out/'healthy_inclusion_audit.tsv',sep='\t',index=False)
# Show how the alternative inclusion policy changes results.
pol=[]
for name,mask in [('healthy_sample_only',sample_eligible),('strict_healthy_people',sample_eligible&~s.person_id.isin(badpeople))]:
 mat=s.loc[mask].groupby('person_id')[qs].max();t=prevalence_table(mat,groups,name);t['samples']=int(mask.sum());pol.append(t)
pd.concat(pol).to_csv(out/'healthy_definition_sensitivity.tsv',sep='\t',index=False)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(13,12));ax.barh(range(len(qs)),prev.prevalence_percent,color='#167D9A');ax.set_yticks(range(len(qs)),[gus_label(q) for q in qs]);ax.invert_yaxis();ax.set_xlim(0,110)
for i,row in prev.iterrows():ax.text(row.prevalence_percent+.7,i,f'{row.prevalence_percent:.2f}% ({row.positive_people}/{row.included_people})',va='center',fontsize=8)
ax.set_xlabel('Prevalence (%)');ax.set_title('Healthy individuals: reference GUS prevalence\nAny eligible healthy sample positive');fig.tight_layout()
for ext in ['png','pdf','svg']:fig.savefig(out/f'healthy_GUS_prevalence.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
# Adapt the established descriptive plotter with an explicit healthy-only source.
code=(root/'scripts/plot_gus_characteristics.py').read_text();code=code.replace("src=root/'results/investigut_gus/id95_cov90_evalue1e-5_24754089_r0/postprocess_24754212_r0/analysis'",'src=Path('+repr(str(out))+')')
code=code.replace('assert len(b)==6559 and b.person_id.nunique()==6559','assert len(b)>0 and b.person_id.nunique()==len(b)').replace('assert len(check)==6559','assert len(check)==len(b)').replace('assert int(b.any_GUS.sum())==4904','assert set(b.disease)=={"healthy"}')
code=code.replace("rows=[]\nfor key", "features.pop('04_disease')  # Disease is fixed by healthy-only inclusion.\nrows=[]\nfor key")
code=code.replace("fig.suptitle('Any-GUS prevalence across individual characteristics'","axs.flat[-1].axis('off')\nfig.suptitle('Healthy individuals: any-GUS prevalence across characteristics'")
code=code.replace('One selected sample per person; 4,904/6,559 positive overall (74.77%).',f'Healthy only; one sample per person; {int(b.any_GUS.sum()):,}/{len(b):,} positive ({100*b.any_GUS.mean():.2f}%).')
code=code.replace("title+'\\nOne selected sample", "'Healthy individuals: '+title+'\\nOne selected sample").replace('Disease groups overlap; country/host differences may reflect study composition.','Healthy-only cohort; country/host differences may reflect study composition.').replace('Disease groups overlap. Missing-data groups retained.','Missing-data groups retained.').replace("'figures':9","'figures':8")
script=out/'plot_healthy_characteristics.py';script.write_text(code);subprocess.run([sys.executable,str(script),'--out',str(out/'characteristics')],check=True)
# Main country figures use any eligible healthy sample, preserving the main prevalence estimand.
for threshold in [0,30]:
 code=(root/'scripts/plot_gus_country_n30.py').read_text().replace("a=r/'results/investigut_gus/id95_cov90_evalue1e-5_24754089_r0/postprocess_24754212_r0/analysis'",'a=Path('+repr(str(out))+')')
 code=code.replace('n=n[n>=30]',f'n=n[n>={threshold}]').replace('v.shape==(23,12)','v.shape==(len(countries),12)').replace('and n.sum()==6524','and (n>=1).all()')
 code=code.replace("'6,524 people | 23 countries with n ≥ 30 | Ordered by included people'",f'f"Healthy only | {{int(n.sum()):,}} people | {{len(countries)}} countries | n ≥ {threshold}"')
 code=code.replace('Reference GUS homolog prevalence by country','Healthy individuals: GUS prevalence by country').replace('Countries with n < 30 excluded from this view.', 'Countries with n < 30 excluded from this view.' if threshold else 'All countries retained; † denotes n < 30.')
 code=code.replace("'minimum_country_n':30",f"'minimum_country_n':{threshold}").replace('GUS_country_prevalence_n_ge_30','healthy_GUS_country_'+('n_ge_30' if threshold else 'all'))
 script=out/f'plot_healthy_country_{threshold}.py';script.write_text(code);subprocess.run([sys.executable,str(script),'--out',str(out/f'country_n{threshold}')],check=True)
# Verify unions and denominators from explicit included rows.
assert np.array_equal(good[list(clades)].max(axis=1),good.any_GUS)
for q in qs:assert int(person[q].sum())==good.loc[good[q].eq(1),'person_id'].nunique()
summary=dict(policy=args.policy,original_samples=len(s),original_people=s.person_id.nunique(),healthy_samples=len(good),healthy_people=len(b),studies=good.study_name.nunique(),countries=good.country.nunique(),positive_any=int(person.any_GUS.sum()),prevalence_any=100*person.any_GUS.mean(),single_sample_positive=int(b.any_GUS.sum()),single_sample_prevalence=100*b.any_GUS.mean(),time_unresolved=int(b.selection_rule.eq('lexical_fallback_time_unresolved').sum()),healthy_labeled_people_excluded_for_case_records=len(set(s.loc[h,'person_id'])&badpeople),validation='passed',input_sha256=hashlib.sha256((src/'gus_sample_presence.tsv').read_bytes()).hexdigest())
(out/'healthy_analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
