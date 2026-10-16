"""Readable query labels; original reference names retained for consistency."""
NAMES={
'NP_416134.1':('Escherichia coli','Loop 1'),
'WP_003467686.1':('Clostridium perfringens','Loop 1'),
'6MVG_A':('Ruminococcus gnavus','FMN-binding'),
'WP_022501306':('Eubacterium sp.','FMN-binding'),
'WP_000966715':('Streptococcus agalactiae','Loop 1'),
'WP_012740861':('Eubacterium eligens','Loop 1'),
'3CMG_A':('Bacteroides fragilis','mini-Loop 1'),
'6NCY_A':('Fusicatenibacter saccharivorans','Hybrid'),
'WP_270534885':('Bacteroides massiliensis','mini-Loop 1'),
'WP_007841259':('Bacteroides dorei','No Loop'),
'WP_217761416':('Prevotella copri','Loop 2')}
CLADE_LABELS={'L1':'Loop 1','FMN':'FMN-binding','mL1':'mini-Loop 1','hybrid':'Hybrid','NL':'No Loop','L2':'Loop 2'}
FPF_NOTE = "* R. gnavus reference GUS: matched representative annotated Roseburia inulinivorans; not R. gnavus abundance."
def gus_label(q,wrap=True,fpf=False):
 if q in NAMES:
  name,clade=NAMES[q];return name+('\n' if wrap else ' ')+f'GUS ({clade})'+(' *' if fpf and q=='6MVG_A' else '')
 if q=='any_GUS':return 'Any reference GUS'
 if q.startswith('clade_'):return CLADE_LABELS[q[6:]]+' GUS union'
 return q
