"""Same capped-rank marker framework as GSE283839, recalibrated in this dataset's candidate TME.
Scores on tumor candidates are diagnostics only; normal labels require a compartment gate.
"""
from pathlib import Path
import json,hashlib,warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import rankdata
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
OUT=Path(__file__).resolve().parents[1]
C=json.loads((OUT/'reference/previous_analysis_config.json').read_text());PANELS=C['panels']
REFS=json.loads((OUT/'reference/reference_gene_sets.json').read_text())
PANELS={**PANELS,**{'reference_'+k:REFS[k] for k in ['Astrocytes','Microglia','CAMs','Monocytes']}}
SEED=20261005;RMAX=1500;NNULL=2000

def rank_weights(x):
 w=x.astype(np.float32).copy();w.eliminate_zeros()
 for i in range(w.shape[0]):
  a,b=w.indptr[i:i+2];w.data[a:b]=np.maximum(0,RMAX-rankdata(-w.data[a:b],method='average'))
 w.eliminate_zeros();return w

class Calibration:
 def __init__(self,genes,det,mean):
  self.genes=genes;self.gi={g:i for i,g in enumerate(genes)};self.det=det;self.mean=mean
  exclude=set(g for gs in list(PANELS.values())+list(REFS.values())+list(C['astro_programs'].values()) for g in gs)
  self.valid=np.array([i for i,g in enumerate(genes) if det[i]>0 and g not in exclude and not g.startswith(('MT-','RPL','RPS')) and g not in ['MALAT1','NEAT1']])
  self.scales=[max(1e-8,np.std(det[self.valid])),max(1e-8,np.std(mean[self.valid]))];self.pools={};self.controls={}
 def design(self,gs,name):
  ix=np.array([self.gi[g] for g in gs if g in self.gi]);n=len(ix);assert n>0
  rng=np.random.default_rng(SEED+int.from_bytes(hashlib.sha256(name.encode()).digest()[:4],'little'))
  draws=np.zeros((NNULL,n),dtype=np.int32)
  for j,i in enumerate(ix):
   if i not in self.pools:
    d=((self.det[self.valid]-self.det[i])/self.scales[0])**2+((self.mean[self.valid]-self.mean[i])/self.scales[1])**2
    self.pools[i]=self.valid[np.argsort(d,kind='stable')[:200]]
   pool=self.pools[i];draws[:,j]=rng.choice(pool,size=NNULL)
   if j:
    bad=(draws[:,j,None]==draws[:,:j]).any(1)
    while bad.any():
     draws[bad,j]=rng.choice(pool,size=bad.sum());bad=(draws[:,j,None]==draws[:,:j]).any(1)
  ctrl=sparse.csr_matrix((np.ones(draws.size,dtype=np.float32),(draws.ravel(),np.repeat(np.arange(NNULL),n))),shape=(len(self.genes),NNULL))
  return ix,ctrl,n*RMAX-n*(n+1)/2
 def score(self,w,gs,name):
  if name not in self.controls:self.controls[name]=self.design(gs,name)
  ix,ctrl,den=self.controls[name];v=np.minimum(w[:,ix].sum(1).A.ravel()/den,1)
  null=np.minimum((w@ctrl).toarray()/den,1)
  p=(1+(null>=v[:,None]-1e-8).sum(1))/(NNULL+1);sd=null.std(1)
  z=np.divide(v-null.mean(1),sd,out=np.zeros(len(v)),where=sd>1e-8)
  return {'ucell':v,'empirical_p':p,'enrichment_z':z}

def main():
 o=pd.read_csv(OUT/'tables/prepared_metadata.csv.gz',index_col=0,low_memory=False)
 thresholds=[];o['cnv_reference_like']=False;o['cnv_elevated_q90']=False;o['cnv_elevated_q95']=False
 for lib,d in o.groupby('library_id',sort=False):
  r=d[d.normal_reference_marker_supported]
  row={'library_id':lib,'n_marker_supported_reference':len(r),'available':bool(d.infercnv_total_cnv.notna().any())}
  if row['available']:
   assert len(r)>=5
   row['limited_reference_n_lt20']=len(r)<20
   for metric in ['infercnv_total_cnv','cnv_mean_arm_fraction']:
    for q in [.9,.95,.99]:row[f'{metric}_q{int(q*100)}']=float(r[metric].max() if len(r)<20 and q>=.95 else r[metric].quantile(q))
   low=(d.infercnv_total_cnv<=row['infercnv_total_cnv_q90'])&(d.cnv_mean_arm_fraction<=row['cnv_mean_arm_fraction_q90'])
   o.loc[d.index,'cnv_reference_like']=low
   for q in [90,95]:
    high=(d.infercnv_total_cnv>row[f'infercnv_total_cnv_q{q}'])&(d.cnv_mean_arm_fraction>row[f'cnv_mean_arm_fraction_q{q}'])
    o.loc[d.index,f'cnv_elevated_q{q}']=high
  thresholds.append(row)
 pd.DataFrame(thresholds).to_csv(OUT/'tables/CNV_operational_thresholds.csv',index=False)
 o['tme_calibration_eligible']=o.qc_pass&o.cnv_reference_like
 # Candidate normal cells determine expression/detection matching; tumor cells do not train these controls.
 smdet=None;smmean=None;n=0
 for path in sorted((OUT/'data').glob('*_counts.h5ad')):
  a=ad.read_h5ad(path);sel=o.loc[a.obs_names,'tme_calibration_eligible'].to_numpy();x=a.X[sel]
  if smdet is None:genes=a.var_names.to_numpy();smdet=np.zeros(a.n_vars);smmean=np.zeros(a.n_vars)
  assert np.array_equal(genes,a.var_names)
  if len(x.indptr)>1:
   smdet+=(x>0).sum(0).A.ravel();y=x.astype(np.float32).multiply(1e4/np.maximum(x.sum(1),1));y.data=np.log1p(y.data);smmean+=y.sum(0).A.ravel();n+=x.shape[0]
  del a,x
 assert n>=100
 np.savez_compressed(OUT/'reference/null_calibration.npz',genes=genes.astype(str),detection=smdet/n,mean_log1p=smmean/n,n_cells=n)
 cal=Calibration(genes,smdet/n,smmean/n);all_ev=[];avail=[]
 for path in sorted((OUT/'data').glob('*_counts.h5ad')):
  a=ad.read_h5ad(path);ev=o.loc[a.obs_names].copy();x=a.X;w=rank_weights(x);lib=ev.library_id.iloc[0]
  sparse.save_npz(OUT/'data'/f'{lib}_rankweights.npz',w)
  for name,gs in PANELS.items():
   out=cal.score(w,gs,name)
   ix=[cal.gi[g] for g in gs if g in cal.gi];ev[name+'__n_detected']=(x[:,ix]>0).sum(1).A.ravel()
   for k,v in out.items():ev[name+'__'+k]=v
   for g in gs:
    i=cal.gi.get(g);avail.append({'library_id':lib,'panel':name,'gene':g,'present':i is not None,'n_detected_all':int((x[:,i]>0).sum()) if i is not None else 0,'n_detected_candidate_TME':int((x[ev.tme_calibration_eligible.values,i]>0).sum()) if i is not None else 0})
  ev.to_csv(OUT/'tables'/f'{lib}_marker_scores.csv.gz');all_ev.append(ev)
  print(lib,'scored','TME candidates',int(ev.tme_calibration_eligible.sum()),'CNV elevated',int(ev.cnv_elevated_q95.sum()),flush=True)
 pd.concat(all_ev).to_csv(OUT/'tables/per_cell_scores.csv.gz')
 pd.DataFrame(avail).to_csv(OUT/'tables/marker_availability.csv',index=False)
 print('Scoring complete. Calibration cells:',n,flush=True)
if __name__=='__main__':main()
