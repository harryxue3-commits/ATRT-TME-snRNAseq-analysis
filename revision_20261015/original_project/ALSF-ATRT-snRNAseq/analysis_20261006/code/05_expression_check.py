"""Exploratory expression check for the CNV-unavailable sample. Never upgrades primary labels."""
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score,balanced_accuracy_score
OUT=Path(__file__).resolve().parents[1];rng=np.random.default_rng(20261005)
e=pd.read_csv(OUT/'tables/cell_annotations.csv.gz',index_col=0,low_memory=False)
blocks=[];meta=[];query=None;queryobs=None
for lib,d in e.groupby('library_id',sort=False):
 a=ad.read_h5ad(OUT/'data'/f'{lib}_counts.h5ad')
 if d.infercnv_total_cnv.isna().all():query=a.X;queryobs=d;continue
 for group,mask in [('Tumor',d.tumor_cnv_supported),('TME',d.compartment_conservative.eq('TME_marker_supported'))]:
  ix=np.flatnonzero(mask);ix=rng.choice(ix,size=min(400,len(ix)),replace=False)
  blocks.append(a.X[ix]);dd=d.iloc[ix][['patient_id','library_id']].copy();dd['training_class']=group;meta.append(dd)
 genes=a.var_names.to_numpy()
x=sparse.vstack(blocks).tocsr();m=pd.concat(meta);y=m.training_class.eq('Tumor').to_numpy()
def norm(z):
 z=z.astype(np.float32).multiply(1e4/np.maximum(z.sum(1),1)).tocsr();z.data=np.log1p(z.data);return z
lx=norm(x);lq=norm(query);mean=lx.mean(0).A.ravel();var=lx.power(2).mean(0).A.ravel()-mean**2
valid=np.array([not g.startswith(('MT-','RPL','RPS')) and g not in ['MALAT1','NEAT1'] for g in genes]);var[~valid]=-1
hvg=np.argsort(var)[-2000:];xx=lx[:,hvg].toarray();qq=lq[:,hvg].toarray()
def fit_predict(train,test):
 scaler=StandardScaler();xt=scaler.fit_transform(train);pca=PCA(n_components=30,random_state=20261005,svd_solver='randomized');return pca.fit_transform(xt),pca.transform(scaler.transform(test))
validation=[]
for pid in sorted(m.patient_id.unique()):
 tr=m.patient_id.ne(pid).to_numpy();te=~tr;z,v=fit_predict(xx[tr],xx[te]);model=LogisticRegression(C=1,max_iter=2000,class_weight='balanced').fit(z,y[tr]);p=model.predict_proba(v)[:,1]
 validation.append({'held_out_patient':pid,'n':int(te.sum()),'balanced_agreement_with_operational_labels':balanced_accuracy_score(y[te],p>=.5),'AUROC_against_operational_labels':roc_auc_score(y[te],p),'note':'Internal agreement, not independently verified tumor accuracy. HVG screening used all training-reference cells.'})
z,v=fit_predict(xx,qq);model=LogisticRegression(C=1,max_iter=2000,class_weight='balanced').fit(z,y)
nn=NearestNeighbors(n_neighbors=20,metric='cosine').fit(z);dist,ix=nn.kneighbors(v);td,ti=nn.kneighbors();support=y[ix].mean(1)
threshold=float(np.quantile(td[y].mean(1),.95));prob=model.predict_proba(v)[:,1]
o=queryobs[['cell_barcode','library_id','sample_id','patient_id','disease_timing','annotation_conservative','qc_pass','tumor_normal_marker_conflict']].copy()
o['expression_classifier_tumor_score']=prob;o['fraction_tumor_reference_neighbors']=support;o['mean_neighbor_cosine_distance']=dist.mean(1);o['training_tumor_distance_q95']=threshold
o['exploratory_tumor_candidate']=o.qc_pass & o.annotation_conservative.eq('Unresolved') & ~o.tumor_normal_marker_conflict & (prob>=.95) & (support>=.9) & (dist.mean(1)<=threshold)
o['evidence_status']='Expression similarity only; CNV unavailable. Excluded from primary and broad CNV-supported tumor sets.'
o.to_csv(OUT/'tables/CNV_unavailable_expression_candidates.csv');pd.DataFrame(validation).to_csv(OUT/'tables/expression_check_internal_validation.csv',index=False)
pd.DataFrame({'gene_symbol':genes[hvg]}).to_csv(OUT/'reference/expression_check_features.csv',index=False)
print('Exploratory expression candidates',int(o.exploratory_tumor_candidate.sum()),'of',len(o),'distance cutoff',threshold,flush=True);print(pd.DataFrame(validation).to_string(index=False),flush=True)
