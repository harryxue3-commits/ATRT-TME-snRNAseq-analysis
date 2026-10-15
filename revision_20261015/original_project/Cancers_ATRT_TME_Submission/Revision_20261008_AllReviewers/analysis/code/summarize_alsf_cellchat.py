from pathlib import Path
import pandas as pd,numpy as np,json,sys
O=Path(__file__).resolve().parents[1];ROOT=O.parents[2];A=ROOT/'ALSF-ATRT-snRNAseq/analysis_20261006';T=O/'tables';D=O/'cellchat'
meta=pd.read_csv(A/'tables/cell_annotations.csv.gz',low_memory=False);manifest=pd.read_csv(A/'cellchat/input_manifest.csv');status=[];eligrows=[];focal=[]
myeloid=['Microglia','Macrophage_like','BAM_like','Myeloid_unresolved'];tum='Tumor_CNV_supported'
for lib in manifest.library_id:
 m=meta[meta.library_id.eq(lib)];pat=m.patient_id.iloc[0];sample=m.sample_id.iloc[0];timing=m.patient_group.iloc[0]
 for tier,minimum in [('conservative',20),('broad',20),('conservative',10)]:
  name=f'ALSF__{lib}__{tier}__min{minimum}';d=D/name;counts=m['cellchat_'+tier].value_counts().drop('Exclude',errors='ignore')
  eligible=counts[counts>=minimum].index.tolist();expected=len(eligible)>=2 and tum in eligible
  sta=pd.read_csv(d/'status.csv').status.iloc[0] if (d/'status.csv').exists() else 'Pending'
  status.append(dict(run=name,library_id=lib,sample_id=sample,patient_id=pat,patient_group=timing,tier=tier,minimum=minimum,status=sta,expected_model=expected,n_eligible_groups=len(eligible),n_input_cells=int(counts.reindex(eligible).sum())))
  for gp,n in counts.items():eligrows.append(dict(run=name,library_id=lib,sample_id=sample,patient_id=pat,tier=tier,minimum=minimum,group=gp,n_cells=int(n),eligible=gp in eligible))
  dt=pd.read_csv(d/'all_tested_interactions.csv',low_memory=False) if (d/'all_tested_interactions.csv').exists() else pd.DataFrame()
  jobs=[]
  for pop in myeloid:jobs.extend([('APP_CD74',tum,pop),('SPP1_CD44',pop,tum),('SPP1_ITGAV_ITGB1',pop,tum)])
  jobs.extend([('NCAM1_FGFR1','Astrocyte',tum),('PTN_PTPRZ1','Astrocyte',tum)])
  for pair,src,tgt in jobs:
   e=src in eligible and tgt in eligible
   z=dt[dt.interaction_name.eq(pair)&dt.source.eq(src)&dt.target.eq(tgt)] if len(dt) else pd.DataFrame()
   state='Ineligible groups' if not e else ('Pending' if sta=='Pending' else ('Not selected by model' if not len(z) else ('BH q < 0.05' if bool(z.iloc[0].BH_q_lt_05) else ('Native nominal only' if bool(z.iloc[0].nominal_p_lt_05) else 'Tested not nominal'))))
   focal.append(dict(run=name,library_id=lib,sample_id=sample,patient_id=pat,patient_group=timing,tier=tier,minimum=minimum,interaction_name=pair,source=src,target=tgt,n_source=int(counts.get(src,0)),n_target=int(counts.get(tgt,0)),eligible=e,status=state,probability=float(z.iloc[0].probability) if len(z) else np.nan,p_native=float(z.iloc[0].p_native) if len(z) else np.nan,p_plus1=float(z.iloc[0].p_plus1) if len(z) else np.nan,q_BH=float(z.iloc[0].q_BH_within_object) if len(z) else np.nan))
st=pd.DataFrame(status);st.to_csv(T/'S18_CellChat_run_status.csv',index=False);pd.DataFrame(eligrows).to_csv(T/'S18_CellChat_group_eligibility.csv',index=False);fo=pd.DataFrame(focal);fo.to_csv(T/'S18_CellChat_focal_candidates_by_specimen.csv',index=False)
rows=[]
for keys,d in fo.groupby(['tier','minimum','interaction_name','source','target'],sort=False):
 e=d[d.eligible];nom=e[e.status.isin(['BH q < 0.05','Native nominal only'])];bh=e[e.status.eq('BH q < 0.05')]
 rows.append(dict(zip(['tier','minimum','interaction_name','source','target'],keys))|dict(eligible_specimens=len(e),eligible_patients=e.patient_id.nunique(),tested_specimens=int(e.p_native.notna().sum()),native_nominal_specimens=len(nom),BH_supported_specimens=len(bh),BH_supported_patients=bh.patient_id.nunique(),BH_patient_ids=';'.join(sorted(bh.patient_id.unique())),BH_specimen_ids=';'.join(sorted(bh.library_id.unique()))))
su=pd.DataFrame(rows);su.to_csv(T/'S18_CellChat_focal_recurrence.csv',index=False)
# Patient aggregation is any eligible specimen passing within-model BH, not an independent statistical test.
pr=fo.groupby(['tier','minimum','interaction_name','source','target','patient_id'],sort=False).agg(n_specimens=('sample_id','size'),n_eligible_specimens=('eligible','sum'),any_BH=('status',lambda x:(x=='BH q < 0.05').any())).reset_index();pr.to_csv(T/'S18_CellChat_focal_by_patient.csv',index=False)
complete=(st.status!='Pending').all() and not st.status.str.startswith('ERROR').any()
summary={'complete':bool(complete),'n_specifications':len(st),'n_complete_models':int(st.status.eq('Complete').sum()),'n_skipped':int(st.status.str.startswith('Skipped').sum()),'n_pending':int(st.status.eq('Pending').sum()),'n_errors':int(st.status.str.startswith('ERROR').sum()),'results':su.to_dict('records')}
if complete:
 main=fo[(fo.tier=='conservative')&(fo.minimum==20)&fo.eligible]
 texts=[]
 for pair,dir in [('APP_CD74','tumor-to-myeloid APP–CD74'),('SPP1_ITGAV_ITGB1','myeloid-to-tumor SPP1–(ITGAV + ITGB1)'),('SPP1_CD44','myeloid-to-tumor SPP1–CD44')]:
  d=main[main.interaction_name.eq(pair)];b=d[d.status.eq('BH q < 0.05')]
  texts.append(f'{dir} passed within-model BH correction in at least one eligible myeloid population in {b.sample_id.nunique()}/{d.sample_id.nunique()} evaluable specimens, representing {b.patient_id.nunique()}/{d.patient_id.nunique()} evaluable patients')
 txt=f'Of 24 prespecified external model configurations, {summary["n_complete_models"]} eligible configurations were evaluated (16 distinct input sets) and {summary["n_skipped"]} were skipped because fewer than two groups or no CNV-supported tumor group met eligibility. In the primary conservative models (minimum 20 nuclei per group), '+ '; '.join(texts)+'. These are within-specimen predictions, and a patient is counted once when at least one specimen passes. Broader-annotation and minimum-ten sensitivities are reported with all eligible denominators in Table S18. Astrocyte candidates were not assessable. Computational recurrence does not establish receptor activation or replicate a microglia-specific mechanism.'
 summary['manuscript_results']=txt
 (O/'data/external_cellchat_results.txt').write_text(txt)
(O/'data/alsf_cellchat_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='results'},indent=2))
if '--require-complete' in sys.argv:assert complete,'External models are incomplete'
