from pathlib import Path
from gus_plot_labels import gus_label, CLADE_LABELS
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);o=p.parse_args().out;o.mkdir(exist_ok=True,parents=True)
r=Path('/quobyte/angelazgrp/2026-ysjin-GUS_public')
a=r/'results/investigut_gus/id95_cov90_evalue1e-5_24754089_r0/postprocess_24754212_r0/analysis'
d=pd.read_csv(a/'gus_prevalence_stratified.tsv',sep='\t'); d=d[(d.stratum=='country')&(d.rule=='any_eligible_sample_within_stratum')]
clades=['L1','FMN','mL1','hybrid','NL','L2'];colors=['#0072B2','#E69F00','#009E73','#CC79A7','#D55E00','#6B65A8']
m=[]
for l in (r/'reference/active_gus.faa').read_text().splitlines():
 if l.startswith('>'):
  q,c,s=l[1:].split();m.append((q,c.split('=')[1],s.split('=')[1]))
m.sort(key=lambda v:clades.index(v[1]));qs=[x[0] for x in m]+["any_GUS"]
n=d[d.query_id=='any_GUS'].set_index('stratum_value').included_people.sort_values(ascending=False,kind='stable');n=n[n>=30];countries=list(n.index)
v=d.pivot(index='stratum_value',columns='query_id',values='prevalence_percent').loc[countries,qs]
pv=d.pivot(index='stratum_value',columns='query_id',values='positive_people').loc[countries,qs]
assert v.shape==(23,12) and v.notna().all().all() and n.sum()==6524
assert np.allclose(v.to_numpy(),pv.to_numpy()/n.to_numpy()[:,None]*100)
names=dict(AUT='Austria',BGD='Bangladesh',CAN='Canada',CHN='China',DEU='Germany',DNK='Denmark',ESP='Spain',EST='Estonia',FIN='Finland',FJI='Fiji',FRA='France',GBR='United Kingdom',HUN='Hungary',ISL='Iceland',ISR='Israel',ITA='Italy',JPN='Japan',KAZ='Kazakhstan',MDG='Madagascar',MNG='Mongolia',NLD='Netherlands',NOR='Norway',PER='Peru',RUS='Russia',SVK='Slovakia',SWE='Sweden',TZA='Tanzania',USA='United States')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(18,16));fig.subplots_adjust(left=.25,right=.90,bottom=.32,top=.89)
im=ax.imshow(v.to_numpy(),cmap='YlGnBu',vmin=0,vmax=100,aspect='auto')
ax.set_yticks(range(len(countries)),[f'{names[c]}  |  n = {n[c]:,}'+(' †' if n[c]<30 else '') for c in countries])
ax.set_xticks(range(12),[gus_label(q) for q,c,s in m]+["Any GUS"],rotation=50,ha='right',rotation_mode='anchor')
for j,(q,c,s) in enumerate(m): ax.get_xticklabels()[j].set_color(colors[clades.index(c)])
for i in range(len(countries)):
 for j in range(12):
  x=v.iloc[i,j];ax.text(j,i,f'{x:.1f}',ha='center',va='center',fontsize=8,color='white' if x>=55 else '#182A35')
ax.set_xticks(np.arange(-.5,12,1),minor=True);ax.set_yticks(np.arange(-.5,len(countries),1),minor=True);ax.grid(which='minor',color='white',linewidth=.6);ax.tick_params(which='both',length=0)
ax.axvline(10.5,color='white',linewidth=3)
for sp in ax.spines.values():sp.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.025);cb.set_label('Prevalence (%)');cb.set_ticks([0,20,40,60,80,100])
fig.suptitle('Reference GUS homolog prevalence by country',fontsize=19,fontweight='bold',y=.965)
fig.text(.25,.926,'6,524 people | 23 countries with n ≥ 30 | Ordered by included people',fontsize=12)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=c,label=CLADE_LABELS[k]) for k,c in zip(clades,colors)],loc='lower left',bbox_to_anchor=(.245,.072),ncol=6,frameon=False,title='Reference GUS clade')
fig.text(.25,.026,'Cells: positive people / all included people in that country (%); any eligible sample positive.\nCountries with n < 30 excluded from this view. Country differences may reflect study composition.\nReference organism names label query sequences, not hit taxonomy.',fontsize=9,linespacing=1.5)
for ext in ['png','pdf','svg']:fig.savefig(o/f'GUS_country_prevalence_n_ge_30.{ext}',dpi=300,bbox_inches='tight')
v.to_csv(o/'country_prevalence_percent.tsv',sep='\t');pv.to_csv(o/'country_positive_people.tsv',sep='\t');n.to_csv(o/'country_denominators.tsv',sep='\t')
d[d.query_id.isin(qs) & d.stratum_value.isin(countries)].to_csv(o/'country_prevalence_source.tsv',sep='\t',index=False)
(o/'validation.json').write_text(json.dumps({'status':'passed','countries':len(countries),'minimum_country_n':30,'queries':12,'sum_country_denominators':int(n.sum()),'rule':'any_eligible_sample_within_stratum','all_percentages_checked_against_counts':True},indent=2))
