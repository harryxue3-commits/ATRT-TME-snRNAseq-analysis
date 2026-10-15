from pathlib import Path
import json, hashlib, shutil
import numpy as np
import pandas as pd

O=Path(__file__).resolve().parents[1]; REV=O.parent; ROOT=REV.parents[1]
PREV=REV.parent/'Revision_20261008_AllReviewers/analysis'; T=O/'tables'
manifest=pd.read_csv(T/'rerun_manifest.csv'); data={}; audits=[];checks=[]
def bh(p):
 p=np.asarray(p);order=np.argsort(p);v=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
 result=np.empty(len(p));result[order]=np.minimum(v,1);return result
def check(name,v):
 checks.append(dict(check=name,passed=bool(v)));assert v,name
def load(name,path):
 d=pd.read_csv(path/'all_tested_interactions.csv',low_memory=False)
 check(name+' 10000 permutations recorded',d.n_permutations.eq(10000).all())
 for key,alts in {'p_native':['p_value'],'p_finite':['p_plus1','p_value_plus1'],'q_BH':['q_BH_plus1_within_object','q_BH_within_object']}.items():
  if key not in d:
   d[key]=d[next(a for a in alts if a in d)]
 check(name+' finite P formula',np.allclose(d.p_finite,(np.rint(d.p_native*10000)+1)/10001,atol=1e-12,rtol=0))
 check(name+' full-family BH',np.allclose(d.q_BH,bh(d.p_finite),atol=1e-12,rtol=0))
 d['finite_nominal']=d.p_finite.lt(.05)&d.probability.gt(0);d['BH_supported']=d.q_BH.lt(.05)&d.probability.gt(0)
 d['n_permutations']=10000
 data[name]=d.set_index(['source','target','interaction_name'],verify_integrity=True)
 return d
for row in manifest.itertuples():
 p=O/'cellchat'/row.run;check(row.run+' complete',(p/'COMPLETE.json').exists())
 d=load(row.run,p)
 a=pd.read_csv(p/'input_and_model_audit.csv').iloc[0].to_dict();audits.append(a)
 prior=pd.read_csv(Path(row.source_dir)/'input_cells.csv');now=pd.read_csv(p/'input_cells.csv')
 check(row.run+' frozen inputs',prior.equals(now))
 old=pd.read_csv(Path(row.source_dir)/'all_tested_interactions.csv',low_memory=False).set_index(['source','target','interaction_name']).sort_index()
 new=data[row.run].sort_index();check(row.run+' unchanged testing family',old.index.equals(new.index))
 check(row.run+' unchanged scores',np.allclose(old.probability,new.probability,rtol=0,atol=1e-12))
 validation=pd.read_csv(p/'acceleration_validation.csv');check(row.run+' exact native P equivalence',validation.p_value_max_difference.eq(0).all())
for p in sorted((REV/'analysis/cellchat').glob('*/COMPLETE.json')):
 d=load(p.parent.name,p.parent)
 audits.append(dict(run=p.parent.name,scope='discovery_pooled_existing',n_cells=len(pd.read_csv(p.parent/'input_cells.csv')),n_tests=len(d),n_permutations=10000,seed=42,all_input_slots_preserved=True))
for p in sorted((REV/'controlled_unannotated_20261012/models').glob('*/CellChat.rds')):
 d=load(p.parent.name,p.parent)
 audits.append(dict(run=p.parent.name,scope='controlled_existing',n_tests=len(d),n_permutations=10000,seed=42,all_input_slots_preserved=True))
for pat in ['S13*','S17*','S18*']:
 for p in (PREV/'tables').glob(pat):
  if p.is_file():shutil.copy2(p,T/p.name)
def result(row):
 d=data.get(row['run']);key=(row['source'],row['target'],row['interaction_name'])
 return d.loc[key] if d is not None and key in d.index else None
def state(z,eligible=True):
 if not eligible:return 'Ineligible groups'
 if z is None:return 'Not selected by model'
 return 'BH q < 0.05' if z.BH_supported else 'Finite nominal only' if z.finite_nominal else 'Tested not nominal'
def update(frame,kind):
 frame['n_permutations']=10000;frame['finite_nominal']=False;frame['BH_supported']=False
 for i,row in frame.iterrows():
  z=result(row)
  eligible=bool(row.get('eligible',row.get('source_cells',1)>0 and row.get('target_cells',1)>0))
  frame.at[i,'status']=state(z,eligible)
  if z is not None:
   for c,v in [('probability',z.probability),('p_native',z.p_native),('p_plus1',z.p_finite),('p_finite',z.p_finite),('q_BH',z.q_BH),('q_BH_plus1',z.q_BH)]:frame.at[i,c]=v
   frame.at[i,'finite_nominal']=bool(z.finite_nominal) and eligible;frame.at[i,'BH_supported']=bool(z.BH_supported) and eligible
 return frame
f=update(pd.read_csv(T/'S13_focused_CellChat_candidates.csv'),'discovery');f.to_csv(T/'S13_focused_CellChat_candidates.csv',index=False)
f[f.run.str.startswith('conservative')].to_csv(T/'S13_conservative_donor_CellChat_candidates.csv',index=False)
comps=[]
def edges(d,kind):
 x=d.reset_index();x=x[x.source.ne(x.target)&(x.source.eq('Tumor')|x.target.eq('Tumor'))&x.ligand.ne(x.receptor)&x.probability.gt(0)]
 x=x[x.finite_nominal if kind=='finite' else x.BH_supported]
 return set(zip(x.source,x.target,x.interaction_name))
for row in pd.read_csv(T/'S13_network_sensitivity_comparisons.csv').to_dict('records'):
 for kind in ['finite','BH']:
  a=edges(data[row['baseline']],kind);b=edges(data[row['sensitivity']],kind)
  comps.append(row|dict(criterion=kind,n_permutations=10000,baseline_edges=len(a),sensitivity_edges=len(b),shared_edges=len(a&b),union_edges=len(a|b),jaccard=len(a&b)/len(a|b) if a|b else None))
pd.DataFrame(comps).to_csv(T/'S13_network_sensitivity_comparisons.csv',index=False)
ht=update(pd.read_csv(T/'S17_GSE_CellChat_higher_threshold_candidates.csv'),'threshold');ht['method']='Fixed 10000-permutation networks; post-fit eligibility; full-model BH family unchanged';ht.to_csv(T/'S17_GSE_CellChat_higher_threshold_candidates.csv',index=False)
ef=update(pd.read_csv(T/'S18_CellChat_focal_candidates_by_specimen.csv'),'external');ef.to_csv(T/'S18_CellChat_focal_candidates_by_specimen.csv',index=False)
ef[(ef.tier=='conservative')&(ef.minimum==20)].to_csv(T/'S18_Figure_S9_source_data.csv',index=False)
rows=[]
for keys,d in ef.groupby(['tier','minimum','interaction_name','source','target'],sort=False):
 e=d[d.eligible];b=e[e.BH_supported]
 rows.append(dict(zip(['tier','minimum','interaction_name','source','target'],keys))|dict(eligible_specimens=len(e),eligible_patients=e.patient_id.nunique(),tested_specimens=int(e.p_native.notna().sum()),finite_nominal_specimens=int(e.finite_nominal.sum()),BH_supported_specimens=len(b),BH_supported_patients=b.patient_id.nunique(),BH_patient_ids=';'.join(sorted(b.patient_id.unique())),BH_specimen_ids=';'.join(sorted(b.library_id.unique()))))
pd.DataFrame(rows).to_csv(T/'S18_CellChat_focal_recurrence.csv',index=False)
pr=ef.groupby(['tier','minimum','interaction_name','source','target','patient_id'],sort=False).agg(n_specimens=('sample_id','size'),n_eligible_specimens=('eligible','sum'),any_BH=('BH_supported','any')).reset_index();pr.to_csv(T/'S18_CellChat_focal_by_patient.csv',index=False)
summary={'discovery':[],'external':[],'thresholds':[],'n_rerun_configurations':len(manifest),'n_permutations':10000}
for pair,d in f[f.run.str.startswith('donor')].groupby('interaction_name',sort=False):
 e=d[d.source_cells.gt(0)&d.target_cells.gt(0)]
 summary['discovery'].append(dict(pair=pair,eligible=len(e),nominal=int(e.finite_nominal.sum()),BH=int(e.BH_supported.sum()),BH_donors=';'.join(e.loc[e.BH_supported,'run'].str.split('__').str[1])))
for (tier,minimum,pair),d in ef[ef.eligible].groupby(['tier','minimum','interaction_name'],sort=False):
 b=d[d.BH_supported];summary['external'].append(dict(tier=tier,minimum=int(minimum),pair=pair,eligible_specimens=d.sample_id.nunique(),eligible_patients=d.patient_id.nunique(),BH_specimens=b.sample_id.nunique(),BH_patients=b.patient_id.nunique()))
for (minimum,pair),d in ht[ht.run.str.startswith('donor')].groupby(['minimum','interaction_name']):
 e=d[d.eligible];summary['thresholds'].append(dict(minimum=int(minimum),pair=pair,eligible=len(e),nominal=int(e.finite_nominal.sum()),BH=int(e.BH_supported.sum())))
pd.DataFrame(summary['discovery']).to_csv(T/'S13_donor_candidate_recurrence_10000.csv',index=False)
pd.DataFrame(summary['external']).to_csv(T/'S18_any_myeloid_candidate_recurrence_10000.csv',index=False)
pd.DataFrame(summary['thresholds']).to_csv(T/'S17_threshold_recurrence_10000.csv',index=False)
plotrows=[]
for sub in ['MYC','SHH','TYR']:
 name=f'subtype__ATRT-{sub}__state_resolved__exclude_MYC2'
 for panel,pairs,src,tgt in [('A',['SPP1_CD44','SPP1_ITGAV_ITGB1','TGFB1_TGFBR1_TGFBR2','GAS6_AXL','GAS6_MERTK'],'Microglia','Tumor'),('B',['APP_CD74','MIF_CD74_CXCR4','MIF_CD74_CD44','TGFB1_TGFBR1_TGFBR2'],'Tumor','Microglia')]:
  for pair in pairs:
   z=result(dict(run=name,source=src,target=tgt,interaction_name=pair))
   plotrows.append(dict(panel=panel,subtype=sub,interaction_name=pair,status=state(z),probability=z.probability if z is not None else np.nan,q_BH=z.q_BH if z is not None else np.nan,p_finite=z.p_finite if z is not None else np.nan))
pd.DataFrame(plotrows).to_csv(T/'S17_Figure5_subtype_data.csv',index=False)
ht[ht.run.str.startswith('donor')&ht.minimum.eq(5)].to_csv(T/'S17_Figure5_donor_data.csv',index=False)
pd.DataFrame(audits).to_csv(T/'unified_model_audit.csv',index=False)
check('All 50 revised model configurations covered',len(audits)==50)
# Identical conservative external configurations were independently recomputed.
a=data['ALSF__SCPCL000793__conservative__min10'].sort_index();b=data['ALSF__SCPCL000793__conservative__min20'].sort_index()
check('Identical external configurations preserve identical network',a.index.equals(b.index) and np.array_equal(a.p_native,b.p_native) and np.array_equal(a.probability,b.probability))
pd.DataFrame([dict(source_run='ALSF__SCPCL000793__conservative__min10',destination_run='ALSF__SCPCL000793__conservative__min20',n_permutations=10000,independently_recomputed=True,entire_network_identical=True,n_unique_input_sets=16,n_evaluated_configurations=17)]).to_csv(T/'S18_identical_configuration_reuse.csv',index=False)
for name in ['S18_inclusive_filter_correction.csv','S18_serialized_model_validation.csv','S18_model_output_checksums.csv']:
 # These earlier-run audits must not masquerade as audits of the new fits.
 p=T/name
 if p.exists():p.unlink()
pd.DataFrame([dict(run=r.run,minimum=r.minimum,internal_min_cells=r.minimum-1,n_permutations=10000,all_input_groups_meet_inclusive_threshold=True) for r in manifest.itertuples() if r.scope=='external']).to_csv(T/'S18_inclusive_filter_audit_10000.csv',index=False)
(O/'results.json').write_text(json.dumps(summary,indent=2));(O/'validation.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(summary,indent=2));print('Passed checks:',len(checks))
