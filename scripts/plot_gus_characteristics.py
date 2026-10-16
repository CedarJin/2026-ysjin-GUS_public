"""Descriptive individual-level GUS strata; matched sample metadata and detection."""
from pathlib import Path
from country_plot_labels import COUNTRY_NAMES
from gus_plot_labels import gus_label
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);out=p.parse_args().out;out.mkdir(parents=True,exist_ok=True)
root=Path('/quobyte/angelazgrp/2026-ysjin-GUS_public');src=root/'results/investigut_gus/id95_cov90_evalue1e-5_24754089_r0/postprocess_24754212_r0/analysis'
b=pd.read_csv(src/'baseline_sample_selection.tsv',sep='\t',low_memory=False);assert len(b)==6559 and b.person_id.nunique()==6559
clades=['L1','FMN','mL1','hybrid','NL','L2'];colors=['#0072B2','#E69F00','#009E73','#CC79A7','#D55E00','#6B65A8'];meta=[]
for l in (root/'reference/active_gus.faa').read_text().splitlines():
 if l.startswith('>'):
  q,c,s=l[1:].split();meta.append((q,c.split('=')[1],s.split('=')[1]))
meta.sort(key=lambda x:clades.index(x[1]));queries=[x[0] for x in meta]+['any_GUS'];allq=pd.read_csv(src/'gus_prevalence.tsv',sep='\t').query_id.tolist()
# Independently match each selected observation to the complete sample matrix.
s=pd.read_csv(src/'gus_sample_presence.tsv',sep='\t',low_memory=False);check=b.merge(s[['study_name','sample_id']+allq],on=['study_name','sample_id'],validate='one_to_one',suffixes=('_selected','_sample'))
assert len(check)==6559
for q in allq:assert np.array_equal(check[q+'_selected'],check[q+'_sample'])
assert int(b.any_GUS.sum())==4904
age=pd.to_numeric(b.age,errors='coerce');bmi=pd.to_numeric(b.BMI,errors='coerce')
bad=bmi.notna()&((bmi<10)|(bmi>70));b.loc[bad,['person_id','study_name','sample_id','age','BMI']].to_csv(out/'BMI_flagged_records.tsv',sep='\t',index=False)
# Broad plausibility screen is an analysis choice, not a correction to source records.
bmi=bmi.mask(bad)
agegroup=pd.cut(age,[-.001,3,18,40,60,np.inf],right=False,labels=['0–<3 years','3–<18 years','18–<40 years','40–<60 years','60+ years']).astype(object).fillna('Missing age')
bmig=pd.cut(bmi,[10,18.5,25,30,70.001],right=False,labels=['BMI <18.5','BMI 18.5–<25','BMI 25–<30','BMI 30+']).astype(object)
bmig.loc[age.isna()]='Age missing (not classified)';bmig.loc[age<18]='Under 18 (not classified)';bmig.loc[(age>=18)&bmi.isna()]='Adult BMI missing / flagged'
assert bmig.notna().all()
# Disease tokens may overlap. Healthy is exact source label; no conflation with control.
tokens=b.disease.fillna('').str.split(';').apply(set);sub=b.disease_subtype.fillna('').str.split(';').apply(set)
disease={}
for lab,tok in [('Healthy label','healthy'),('CRC','CRC'),('Adenoma','adenoma'),('IBD (all)','IBD'),('T2D','T2D'),('T1D','T1D'),('Rheumatoid arthritis','RA'),('Fatty liver','fatty_liver'),('Cirrhosis','cirrhosis'),('Hypertension','hypertension'),('Impaired glucose tolerance','IGT'),('STEC infection','STEC'),('C. difficile infection','CDI')]:
 disease[lab]=(b.disease.eq('healthy') if tok=='healthy' else tokens.apply(lambda x:tok in x))
as_mask=b.study_condition.eq('AS')|sub.apply(lambda x:'AS' in x)
b.loc[as_mask & disease['Rheumatoid arthritis'],['person_id','study_name','sample_id','study_condition','disease','disease_subtype']].to_csv(out/'disease_label_conflicts_AS_RA.tsv',sep='\t',index=False)
disease['Rheumatoid arthritis']=disease['Rheumatoid arthritis']&~as_mask
disease['Ankylosing spondylitis (AS)']=as_mask
disease['Crohn disease (IBD subtype)']=disease['IBD (all)']&sub.apply(lambda x:'CD' in x)
disease['Ulcerative colitis (IBD subtype)']=disease['IBD (all)']&sub.apply(lambda x:'UC' in x)
covered=pd.concat(disease,axis=1).any(axis=1);disease['Other recorded conditions']=~covered&b.disease.notna();disease['Missing disease']=b.disease.isna()
def masks(series,order=None):
 if order is None:order=list(series.drop_duplicates())
 return {str(x):series.eq(x) for x in order}
features={
 '01_age':('Age',masks(agegroup,['0–<3 years','3–<18 years','18–<40 years','40–<60 years','60+ years','Missing age'])),
 '02_sex':('Recorded sex',masks(b.gender.map({'female':'Female','male':'Male'}).fillna('Missing sex'),['Female','Male','Missing sex'])),
 '03_BMI':('Adult BMI',masks(bmig,['BMI <18.5','BMI 18.5–<25','BMI 25–<30','BMI 30+','Adult BMI missing / flagged','Under 18 (not classified)','Age missing (not classified)'])),
 '04_disease':('Recorded disease (groups overlap)',disease),
 '05_antibiotics':('Current antibiotic use',masks(b.antibiotics_current_use.map({'yes':'Yes','no':'No'}).fillna('Missing'),['Yes','No','Missing'])),
 '06_smoking':('Recorded smoking',masks(b.smoker.map({'yes':'Yes','no':'No'}).fillna('Missing'),['Yes','No','Missing'])),
 '07_country':('Country',masks(b.country,b.country.value_counts().index)),
 '08_lifestyle':('Non-westernized label',masks(b.non_westernized.map({'yes':'Non-westernized: yes','no':'Non-westernized: no'}))) }
def wilson(k,n):
 z=1.959963984540054;p=k/n;den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den;return 100*(mid-half),100*(mid+half)
rows=[]
for key,(title,groups) in features.items():
 if key!='04_disease':assert np.all(pd.concat(groups,axis=1).sum(axis=1)==1)
 for label,mask in groups.items():
  g=b.loc[mask];n=len(g)
  if not n:continue
  for q in allq:
   k=int(g[q].sum());lo,hi=wilson(k,n);rows.append(dict(feature=key,group=label,query_id=q,positive_people=k,included_people=n,prevalence_percent=100*k/n,ci95_lower=lo,ci95_upper=hi,studies=g.study_name.nunique(),rule='one_selected_sample_with_matched_metadata'))
t=pd.DataFrame(rows);t.to_csv(out/'all_stratified_prevalence.tsv',sep='\t',index=False)
assign=b[['person_id','study_name','sample_id','selection_rule']].copy()
for key,(_,groups) in features.items():
 for label,mask in groups.items():assign[key+'::'+label]=mask.astype(int)
assign.to_csv(out/'individual_group_membership.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
book=PdfPages(out/'GUS_characteristics_all_panels.pdf')
for key,(title,groups) in features.items():
 tt=t[(t.feature==key)&t.query_id.isin(queries)];order=[g for g in groups if g in set(tt.group)]
 v=tt.pivot(index='group',columns='query_id',values='prevalence_percent').loc[order,queries];ns=tt[tt.query_id=='any_GUS'].set_index('group').included_people
 fig,ax=plt.subplots(figsize=(19,max(9,len(order)*.39+6)));fig.subplots_adjust(left=.29,right=.91,bottom=.40 if len(order)<10 else .30,top=.85)
 im=ax.imshow(v,cmap='YlGnBu',vmin=0,vmax=100,aspect='auto');ax.set_yticks(range(len(order)),[f'{COUNTRY_NAMES[g] if key == "07_country" else g} | n={ns[g]:,}'+(' †' if ns[g]<30 else '') for g in order])
 ax.set_xticks(range(12),[gus_label(q) for q,c,_ in meta]+['Any GUS'],rotation=45,ha='right',rotation_mode='anchor')
 for i in range(len(order)):
  for j in range(12):ax.text(j,i,f'{v.iloc[i,j]:.1f}',ha='center',va='center',fontsize=8,color='white' if v.iloc[i,j]>=55 else '#172A35')
 for j,(_,c,_) in enumerate(meta):ax.get_xticklabels()[j].set_color(colors[clades.index(c)])
 ax.axvline(10.5,color='white',lw=3);ax.tick_params(length=0);fig.colorbar(im,ax=ax,fraction=.025,pad=.02,label='Prevalence (%)')
 ax.set_title(title+'\nOne selected sample per person; metadata and detection from the same sample',loc='left',pad=20,fontsize=13)
 fig.text(.29,.03,'Descriptive, unadjusted. Cells include zeros; missing metadata are shown separately.\n† n < 30. Disease groups overlap; country/host differences may reflect study composition.\nBMI: adults only; flagged implausible values excluded from BMI strata. No causal interpretation.',fontsize=9)
 for ext in ['png','pdf','svg']:fig.savefig(out/f'{key}.{ext}',dpi=220,bbox_inches='tight')
 book.savefig(fig,bbox_inches='tight');plt.close(fig)
# Eight-panel summary with Wilson intervals, same denominators as heatmaps.
fig,axs=plt.subplots(2,4,figsize=(24,13));fig.subplots_adjust(left=.09,right=.98,bottom=.11,top=.90,wspace=.90,hspace=.40)
for ax,(key,(title,groups)) in zip(axs.flat,features.items()):
 tt=t[(t.feature==key)&(t.query_id=='any_GUS')].set_index('group');order=[g for g in groups if g in tt.index];tt=tt.loc[order];y=np.arange(len(tt));x=tt.prevalence_percent.to_numpy()
 ax.errorbar(x,y,xerr=np.vstack([x-tt.ci95_lower.to_numpy(),tt.ci95_upper.to_numpy()-x]),fmt='o',ms=3,color='#126782',capsize=2,lw=1)
 ax.set_yticks(y,[f'{COUNTRY_NAMES[g] if key == "07_country" else g} (n={int(n):,})' for g,n in zip(tt.index,tt.included_people)],fontsize=7 if len(tt)>10 else 8);ax.set_ylim(len(tt)-.4,-.6);ax.set_xlim(0,100);ax.set_xticks([0,50,100]);ax.grid(axis='x',alpha=.2);ax.set_title(title,fontsize=12,loc='left');ax.set_xlabel('Any-GUS prevalence (%)');ax.spines[['top','right']].set_visible(False)
fig.suptitle('Any-GUS prevalence across individual characteristics',fontsize=21,fontweight='bold')
fig.text(.09,.025,'One selected sample per person; 4,904/6,559 positive overall (74.77%). Points: observed proportions; bars: nominal Wilson 95% intervals.\nUnadjusted descriptive estimates, not population-representative or causal effects. Disease groups overlap. Missing-data groups retained.',fontsize=12)
for ext in ['png','pdf','svg']:fig.savefig(out/f'00_any_GUS_overview.{ext}',dpi=200,bbox_inches='tight')
book.savefig(fig,bbox_inches='tight');plt.close(fig);book.close()
# Within-study disease vs study-defined controls. Descriptive only; common study is required.
contrasts=[]
for condition in ['CRC','adenoma','IBD','T2D','RA','fatty_liver','cirrhosis']:
 for study,g in b.groupby('study_name'):
  case=g[g.study_condition.eq(condition)];ctrl=g[g.study_condition.eq('control')]
  if len(case)<10 or len(ctrl)<10:continue
  for q in allq:
   contrasts.append(dict(condition=condition,study=study,query_id=q,case_n=len(case),case_positive=int(case[q].sum()),control_n=len(ctrl),control_positive=int(ctrl[q].sum()),case_prevalence=case[q].mean()*100,control_prevalence=ctrl[q].mean()*100,difference_percentage_points=(case[q].mean()-ctrl[q].mean())*100,control_definition='study_condition=control; not necessarily free of comorbidity'))
ct=pd.DataFrame(contrasts);ct.to_csv(out/'disease_within_study_descriptive.tsv',sep='\t',index=False)
# Re-read exported table and check all source counts against explicit membership.
rt=pd.read_csv(out/'all_stratified_prevalence.tsv',sep='\t')
for row in rt.itertuples():
 mask=features[row.feature][1][row.group];assert row.included_people==int(mask.sum());assert row.positive_people==int(b.loc[mask,row.query_id].sum());assert row.ci95_lower-1e-10<=row.prevalence_percent<=row.ci95_upper+1e-10
(out/'validation.json').write_text(json.dumps({'status':'passed','individuals':len(b),'single_sample_positive':int(b.any_GUS.sum()),'BMI_flagged':int(bad.sum()),'table_rows':len(rt),'figures':9,'sample_presence_crosscheck':'all 18 outcomes matched','mutually_exclusive_group_coverage':'passed except intentionally overlapping disease groups','inputs':{n:hashlib.sha256((src/n).read_bytes()).hexdigest() for n in ['baseline_sample_selection.tsv','gus_sample_presence.tsv']}},indent=2))
print(t[t.query_id=='any_GUS'][['feature','group','positive_people','included_people','prevalence_percent']].to_string(index=False))
