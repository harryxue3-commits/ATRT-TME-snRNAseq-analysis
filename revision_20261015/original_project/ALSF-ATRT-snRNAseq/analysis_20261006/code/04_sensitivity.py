from pathlib import Path
import importlib.util,json
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
OUT=Path(__file__).resolve().parents[1]
def module(file,name):
 spec=importlib.util.spec_from_file_location(name,OUT/'code'/file);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
S=module('02_scores.py','scores');A=module('03_annotate.py','annotation')
e=pd.read_csv(OUT/'tables/cell_annotations.csv.gz',index_col=0,low_memory=False)
fit=np.load(OUT/'reference/null_calibration.npz');cal=S.Calibration(fit['genes'],fit['detection'],fit['mean_log1p'])
records=[]
mapping={'Astrocyte':'Astrocyte','Microglia':'Microglia','Macrophage':'Macrophage_like','BAM':'BAM_like','Monocyte':'Monocyte_like'}
for lib,d in e.groupby('library_id',sort=False):
 use=d.qc_pass&(d.cnv_reference_like|d.normal_reference_without_cnv);idx=np.flatnonzero(use);ev=d.iloc[idx].copy()
 if not len(ev):continue
 w=sparse.load_npz(OUT/'data'/f'{lib}_rankweights.npz')[idx];a=ad.read_h5ad(OUT/'data'/f'{lib}_counts.h5ad',backed='r')
 for panel,label in mapping.items():
  for gene in S.PANELS[panel]:
   r=cal.score(w,[g for g in S.PANELS[panel] if g!=gene],panel+'_minus_'+gene);p=ev.copy()
   for key,val in r.items():p[panel+'__'+key]=val
   detected=(a[idx,[gene]].X.toarray().ravel()>0).astype(int);p[panel+'__n_detected']-=detected
   for anchor_col,genes in [('astro_anchor_detected',S.C['rules']['astro_anchor']),('microglia_anchor_detected',S.C['rules']['microglia_anchor']),('macrophage_anchor_detected',['CD163','MRC1','ITGA4','F13A1','STAB1'])]:
    if (panel=='Astrocyte' and anchor_col=='astro_anchor_detected') or (panel=='Microglia' and anchor_col=='microglia_anchor_detected') or (panel=='Macrophage' and anchor_col=='macrophage_anchor_detected'):
     if gene in genes:p[anchor_col]-=detected
   pred=A.classify_tme(p)
   z=pd.DataFrame({'cell_id':ev.index,'panel':panel,'omitted_gene':gene,'baseline':ev.tme_marker_k2.values,'perturbed':pred.values});z['unchanged']=z.baseline.eq(z.perturbed);records.append(z)
 a.file.close();print(lib,'leave-one-marker-out completed',len(ev),'eligible cells',flush=True)
r=pd.concat(records);r.to_csv(OUT/'tables/leave_one_marker_out.csv.gz',index=False)
summary=[]
for panel,label in mapping.items():
 z=r[r.panel.eq(panel) & r.baseline.eq(label)];st=z.groupby('cell_id').unchanged.mean()
 e['LOO_'+panel+'_stability']=e.index.map(st)
 summary.append({'panel':panel,'n_conservative_cells':len(st),'mean_own_panel_LOO_stability':st.mean(),'n_stable_at_least_80pct':int(st.ge(.8).sum()),'limitation':'Marker score/anchor perturbation only; reference, CNV gates and neighbors fixed'})
pd.DataFrame(summary).to_csv(OUT/'tables/marker_sensitivity_summary.csv',index=False)
counts=[]
for lib,d in e.groupby('library_id'):
 for k in [1,2,3]:
  for label in list(mapping.values())+['Myeloid_unresolved']:
   counts.append({'library_id':lib,'minimum_markers':k,'population':label,'n_cells':int(d[f'tme_marker_k{k}'].eq(label).sum())})
pd.DataFrame(counts).to_csv(OUT/'tables/minimum_marker_sensitivity.csv',index=False)
t=pd.read_csv(OUT/'tables/CNV_operational_thresholds.csv').set_index('library_id');tum=[]
for lib,d in e.groupby('library_id'):
 hi=d.infercnv_total_cnv.gt(t.loc[lib,'infercnv_total_cnv_q99'])&d.cnv_mean_arm_fraction.gt(t.loc[lib,'cnv_mean_arm_fraction_q99'])
 e.loc[d.index,'tumor_passes_q99_sensitivity']=d.tumor_cnv_supported&hi
 tum.append({'library_id':lib,'conservative_q95':int(d.tumor_cnv_supported.sum()),'broad_q90':int(d.tumor_cnv_broad.sum()),'stricter_q99':int((d.tumor_cnv_supported&hi).sum())})
pd.DataFrame(tum).to_csv(OUT/'tables/tumor_CNV_threshold_sensitivity.csv',index=False)
e.to_csv(OUT/'tables/cell_annotations.csv.gz');e.to_csv(OUT/'tables/cell_annotations.csv')
print(pd.DataFrame(summary).to_string(index=False),flush=True)
