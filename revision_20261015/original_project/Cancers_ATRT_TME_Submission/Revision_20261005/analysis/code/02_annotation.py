#!/usr/bin/env python3
"""TME-restricted marker calibration, reference-signature comparison and refinement."""
import os
os.environ.setdefault('NUMBA_CACHE_DIR','/private/tmp/atrt_numba')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/atrt_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/private/tmp/atrt_cache')
from pathlib import Path
import json, importlib.metadata
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, mannwhitneyu
from sklearn.metrics import adjusted_rand_score
import anndata as ad
import scanpy as sc
import networkx as nx

OUT=Path(__file__).resolve().parents[1]
C=json.loads((OUT/'code/analysis_config.json').read_text())
PANELS=C['panels']; RNG=np.random.default_rng(C['seed'])

def bh(p):
 p=np.asarray(p);o=np.argsort(p);q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];r=np.empty(len(p));r[o]=np.minimum(q,1);return r

def ucell(ranks,indices):
 n=len(indices); lo=n*(n+1)/2;return np.clip(1-(ranks[:,indices].sum(1)-lo)/(n*C['max_rank']-lo),0,1)

class Calibration:
 def __init__(self,x,genes,exclude):
  self.x=x;self.g=genes;self.gi={g:i for i,g in enumerate(genes)}
  self.r=np.minimum(rankdata(-x,axis=1,method='average'),C['max_rank']).astype(np.float32)
  self.det=(x>0).mean(0); self.mean=np.log1p(x*(1e4/x.sum(1))[:,None]).mean(0)
  self.valid=np.array([i for i,g in enumerate(genes) if g not in exclude and not g.startswith(('MT-','RPL','RPS')) and g not in ['MALAT1','NEAT1'] and self.det[i]>0])
  self.scales=np.array([np.std(self.det[self.valid]),np.std(self.mean[self.valid])]);self.pools={};self.cache={}
 def score(self,gs,name):
  gs=[g for g in gs if g in self.gi];ix=np.array([self.gi[g] for g in gs]);n=len(ix)
  if not n:return None
  observed=ucell(self.r,ix); detected=(self.x[:,ix]>0).sum(1)
  # Reproducible per-panel RNG; no identical controls within a sampled signature.
  import hashlib
  seed=int.from_bytes(hashlib.sha256(name.encode()).digest()[:4],'little')+C['seed'];rng=np.random.default_rng(seed)
  pools=[]
  for i in ix:
   if i not in self.pools:
    d=((self.det[self.valid]-self.det[i])/self.scales[0])**2+((self.mean[self.valid]-self.mean[i])/self.scales[1])**2
    self.pools[i]=self.valid[np.argsort(d,kind='stable')[:C['null_pool_size']]]
   pools.append(self.pools[i])
  draws=np.zeros((C['n_null'],n),dtype=np.int32)
  for j,pool in enumerate(pools):
   draws[:,j]=rng.choice(pool,size=C['n_null'])
   if j:
    bad=np.any(draws[:,j,None]==draws[:,:j],axis=1)
    while bad.any():
     draws[bad,j]=rng.choice(pool,size=bad.sum());bad=np.any(draws[:,j,None]==draws[:,:j],axis=1)
  null=np.empty((self.x.shape[0],C['n_null']),dtype=np.float32)
  lo=n*(n+1)/2
  for start in range(0,C['n_null'],100):
   dr=draws[start:start+100]
   null[:,start:start+len(dr)]=np.clip(1-(self.r[:,dr].sum(2)-lo)/(n*C['max_rank']-lo),0,1)
  p=(1+(null>=observed[:,None]-1e-8).sum(1))/(C['n_null']+1)
  sd=null.std(1);z=np.divide(observed-null.mean(1),sd,out=np.zeros(len(observed)),where=sd>1e-8)
  out={'ucell':observed,'n_detected':detected,'null_p95':np.quantile(null,.95,axis=1,method='higher'),'empirical_p':p,'enrichment_z':z}
  self.cache[name]=out
  return out

def embedding(a,prefix):
 b=a.copy();sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b)
 # Expression is never taken from this scaled object for marker testing.
 sc.pp.highly_variable_genes(b,flavor='seurat',n_top_genes=min(2000,b.n_vars))
 use=b.var.highly_variable & ~b.var_names.str.startswith(('MT-','RPL','RPS'))
 b=b[:,use].copy();sc.pp.scale(b,max_value=10);sc.pp.pca(b,n_comps=min(30,b.n_obs-2),random_state=C['seed'])
 sc.pp.neighbors(b,n_neighbors=min(15,b.n_obs-1),n_pcs=min(30,b.n_obs-2),metric='cosine',random_state=C['seed'])
 sc.tl.umap(b,min_dist=.35,random_state=C['seed'])
 conn=b.obsp['connectivities']; graph=nx.from_scipy_sparse_array(conn)
 labels={}
 for res in [.4,.8,1.2]:
  groups=nx.community.louvain_communities(graph,resolution=res,seed=C['seed'],weight='weight')
  lab=np.zeros(b.n_obs,dtype=int)
  for j,inds in enumerate(sorted(groups,key=lambda g:min(g))):lab[list(inds)]=j
  labels[f'{prefix}_louvain_{res}']=lab
 coord=pd.DataFrame(b.obsm['X_umap'],index=a.obs_names,columns=[prefix+'_umap1',prefix+'_umap2'])
 for k,v in labels.items():coord[k]=v
 a.obsm['X_'+prefix+'_pca']=b.obsm['X_pca'];a.obsm['X_'+prefix+'_umap']=b.obsm['X_umap']
 return coord,conn,labels

def main():
 a=ad.read_h5ad(OUT/'data/TME_raw_counts.h5ad');x=a.X.toarray();genes=np.array(a.var_names);gi={g:i for i,g in enumerate(genes)}
 assert a.n_obs==516 and a.obs.Final_Annotation_Focus.eq('TME').all()
 refs=json.loads((OUT/'reference/reference_gene_sets.json').read_text())
 relevant={k:v for k,v in refs.items() if k in ['CAMs','Microglia','Monocytes','Astrocytes','Oligodendrocytes','Neurons','Endothelial cells','Mural cells'] or k.startswith('developmental_')}
 excluded=set(g for gs in list(PANELS.values())+list(refs.values())+list(C['astro_programs'].values()) for g in gs)
 microsets=pd.read_csv(OUT.parents[2]/'Supplementary/Code/04_microglia_state_analysis/microglia_state_program_gene_sets.tsv',sep='\t')
 excluded.update(microsets.gene)
 cal=Calibration(x,genes,excluded)
 evidence=a.obs.copy();evidence.index.name='Cell'
 scores=[]; availability=[]
 for name,gs in {**PANELS,**{'reference_'+k:v for k,v in relevant.items()}}.items():
  result=cal.score(gs,name)
  for field,values in result.items():evidence[name+'__'+field]=values
  for g in gs:
   i=gi.get(g);availability.append({'panel':name,'gene':g,'in_matrix':i is not None,'detected_TME':int((x[:,i]>0).sum()) if i is not None else 0,'detected_published_astrocytes':int((x[a.obs.Final_Annotation.eq('Astrocytes'),i]>0).sum()) if i is not None else 0,'detected_published_microglia':int((x[a.obs.Final_Annotation.eq('Microglia'),i]>0).sum()) if i is not None else 0,'used_in_identity_score':i is not None,'exclusion_reason':'' if i is not None else 'absent from raw matrix'})
  print('Calibrated',name,'p<.05',int((result['empirical_p']<.05).sum()),flush=True)
 pd.DataFrame(availability).to_csv(OUT/'tables/marker_availability.csv',index=False)
 coord,conn,clusters=embedding(a,'TME');evidence=evidence.join(coord)
 detection=lambda gs:np.sum(x[:,[gi[g] for g in gs if g in gi]]>0,axis=1)
 evidence['astro_anchor_detected']=detection(C['rules']['astro_anchor'])
 evidence['astro_nonmetabolic_anchor_detected']=detection(C['rules']['astro_nonmetabolic_anchor'])
 evidence['microglia_anchor_detected']=detection(C['rules']['microglia_anchor'])
 evidence['macrophage_anchor_detected']=detection(['CD163','MRC1','ITGA4','F13A1','STAB1'])
 evidence['qc_exclusion']= (evidence.raw_n_genes<200)|(evidence.raw_library_size<=0)|(evidence.raw_percent_mt>10)
 evidence['qc_flag']=np.where(evidence.qc_exclusion,'Fails count/mitochondrial safeguard','Pass prior QC; no added hard exclusion')
 # Neighbors are computed without published labels or reference-score labels.
 preliminary_astro=(evidence['Astrocyte__empirical_p']<.05)&(evidence['Astrocyte__n_detected']>=1)&(evidence.astro_anchor_detected>=1)
 preliminary_my=evidence['Myeloid__n_detected']>=2
 weights=conn.copy();weights.setdiag(0);den=np.asarray(weights.sum(1)).ravel()
 evidence['astro_neighborhood_support']=np.asarray(weights@preliminary_astro.astype(float)).ravel()/np.maximum(den,1e-10)
 preliminary_nonmet=(evidence['Astrocyte_no_metabolic__empirical_p']<.05)&(evidence['Astrocyte_no_metabolic__n_detected']>=1)&(evidence.astro_nonmetabolic_anchor_detected>=1)
 evidence['astro_nonmetabolic_neighborhood_support']=np.asarray(weights@preliminary_nonmet.astype(float)).ravel()/np.maximum(den,1e-10)
 evidence['myeloid_neighborhood_support']=np.asarray(weights@preliminary_my.astype(float)).ravel()/np.maximum(den,1e-10)
 alternatives=['OPC','Oligodendrocyte','Neuron','Endothelial','Pericyte','T_cell','NK','B_cell','Dendritic']
 def supported(name,k=2):return (evidence[name+'__empirical_p']<.05)&(evidence[name+'__n_detected']>=k)
 def competing(target):
  z=evidence[target+'__enrichment_z'];bad=np.zeros(len(evidence),dtype=bool)
  for p in alternatives+(['Myeloid'] if target.startswith('Astrocyte') else ['Astrocyte']):
   bad |= supported(p,2)&(evidence[p+'__enrichment_z']>z+1)
  return bad
 def classify(k=2,astro_panel='Astrocyte',overrides=None,strict_backbone=False):
  prior=evidence.Final_Annotation.astype(str)
  labels=prior.copy();confidence=pd.Series('Original TME control label',index=evidence.index)
  reason=pd.Series('Retained published control population; competing markers audited',index=evidence.index)
  support=lambda n,kk=2:supported(n,kk) if not overrides or n not in overrides else ((overrides[n]['empirical_p']<.05)&(overrides[n]['n_detected']>=kk))
  anchor=evidence.astro_nonmetabolic_anchor_detected if astro_panel=='Astrocyte_no_metabolic' else evidence.astro_anchor_detected
  mg_anchor=evidence.microglia_anchor_detected.copy()
  if overrides:
   for pn,ov in overrides.items():
    omitted=ov.get('omitted_gene')
    if pn==astro_panel and omitted in (C['rules']['astro_nonmetabolic_anchor'] if astro_panel=='Astrocyte_no_metabolic' else C['rules']['astro_anchor']):anchor=anchor-(x[:,gi[omitted]]>0)
    if pn=='Microglia' and omitted in C['rules']['microglia_anchor']:mg_anchor=mg_anchor-(x[:,gi[omitted]]>0)
  astro_context=(evidence['astro_nonmetabolic_neighborhood_support' if astro_panel=='Astrocyte_no_metabolic' else 'astro_neighborhood_support']>=.5)|support('reference_Astrocytes',2)
  astro=support(astro_panel,k)&(anchor>=1)&~competing(astro_panel)&astro_context&~evidence.qc_exclusion
  my_context=(evidence.myeloid_neighborhood_support>=.5)|support('reference_Microglia',2)|support('reference_CAMs',2)|support('reference_Monocytes',2)
  # A broad backbone is a detection safeguard, not a second correlated null test.
  # Subtype calls below still require calibrated panel enrichment.
  backbone=support('Myeloid',2) if strict_backbone else evidence['Myeloid__n_detected']>=2
  my=backbone&my_context&~evidence.qc_exclusion&~competing('Myeloid')
  mg=my&support('Microglia',k)&(mg_anchor>=1)
  mac_anchor=evidence.macrophage_anchor_detected.copy()
  if overrides and 'Macrophage' in overrides and overrides['Macrophage'].get('omitted_gene') in ['CD163','MRC1','ITGA4','F13A1','STAB1']:
   mac_anchor=mac_anchor-(x[:,gi[overrides['Macrophage']['omitted_gene']]]>0)
  mac=my&support('Macrophage',max(2,k))&(mac_anchor>=1)
  bam=my&support('BAM',max(2,k))
  mono=my&support('Monocyte',max(2,k))
  candidates=prior.eq('Microglia')|my
  labels[candidates]='Myeloid unresolved';confidence[candidates]='Uncertain'
  reason[candidates]='Myeloid source label or backbone; insufficient or overlapping subtype evidence'
  labels[prior.eq('Astrocytes')]='Astrocyte uncertain';confidence[prior.eq('Astrocytes')]='Uncertain';reason[prior.eq('Astrocytes')]='Insufficient coherent astrocyte evidence for conservative subset'
  labels[astro]='Astrocyte';confidence[astro]='High (operational marker criteria)';reason[astro]='Astrocyte panel, anchor, neighborhood and competing-lineage safeguards pass'
  # Shared activated profiles are not interpreted as proof of infiltrating macrophages.
  micro=mg&~(mac|bam|mono)
  labels[micro]='Microglia';confidence[micro]='High (operational marker criteria)';reason[micro]='Microglial panel and anchor with myeloid backbone; no competing macrophage/monocyte support'
  mon=mono&~mg&~bam
  labels[mon]='Monocyte-like';confidence[mon]='Moderate';reason[mon]='Positive monocyte panel within myeloid lineage'
  b=bam&~mg&~mono
  labels[b]='Border-associated macrophage-like';confidence[b]='Moderate';reason[b]='Positive BAM panel; anatomical location and origin unproven'
  m=mac&~mg&~bam&~mono
  labels[m]='Macrophage-like';confidence[m]='Moderate';reason[m]='Positive macrophage panel; origin unproven'
  refcam=support('reference_CAMs',2);refmono=support('reference_Monocytes',2)
  confidence[(b|m)&refcam]='High (marker and reference signatures)'
  confidence[mon&refmono]='High (marker and reference signatures)'
  both=mg&(mac|bam|mono);reason[both]='Microglial and macrophage/monocyte panels both supported; identity unresolved'
  conflict=astro&my;labels[conflict]='Mixed lineage uncertain';confidence[conflict]='Uncertain';reason[conflict]='Coherent astrocytic and myeloid evidence; review potential doublet/contamination'
  labels[evidence.qc_exclusion]='QC uncertain';confidence[evidence.qc_exclusion]='Uncertain';reason[evidence.qc_exclusion]='Additional count/mitochondrial safeguard fails'
  return labels,confidence,reason
 for k in [1,2,3]:
  lab,conf,reason=classify(k);evidence[f'annotation_k{k}']=lab
  if k==2:evidence['tme_annotation_revised']=lab;evidence['annotation_confidence']=conf;evidence['annotation_reason']=reason
 nonmet,_,_=classify(2,'Astrocyte_no_metabolic');evidence['astrocyte_nonmetabolic_supported']=nonmet.eq('Astrocyte')
 strict,_,_=classify(2,strict_backbone=True);evidence['annotation_strict_backbone_null']=strict
 evidence['included_in_conservative_astrocytes']=evidence.Final_Annotation.eq('Astrocytes')&evidence.tme_annotation_revised.eq('Astrocyte')
 evidence['included_in_conservative_microglia']=evidence.Final_Annotation.eq('Microglia')&evidence.tme_annotation_revised.eq('Microglia')
 evidence['new_TME_astrocyte_candidate']=~evidence.Final_Annotation.eq('Astrocytes')&evidence.tme_annotation_revised.eq('Astrocyte')
 evidence['new_TME_myeloid_candidate']=~evidence.Final_Annotation.eq('Microglia')&evidence.tme_annotation_revised.isin(['Microglia','Macrophage-like','Monocyte-like','Border-associated macrophage-like'])
 predictions=pd.read_csv(OUT.parents[2]/'PathwayAnalysis/PathwayAnalysisV2/GSE283839_ATRT_RNA_V3_non_malignant_TME_27B_cell_type_predictions_opt.csv',usecols=['cell_name','predicted_cell_type']).set_index('cell_name')
 assert predictions.index.is_unique and set(predictions.index)==set(evidence.index)
 evidence=evidence.join(predictions);evidence.rename(columns={'predicted_cell_type':'C2S_original_prediction'},inplace=True)
 old=pd.read_csv(OUT.parents[2]/'ATRT_TME_Manuscript_Synthesis_2026-07-15/artifacts/astrocyte_state_ucell/cell_program_scores.csv')
 none=~old.groupby('cell_id').positive.any();evidence['below_all_original_astro_state_thresholds']=evidence.index.map(none).fillna(False)
 # Full myeloid subclustering, with the author-labeled population retained for audit.
 myidx=np.flatnonzero(evidence.Final_Annotation.eq('Microglia')|preliminary_my)
 b=a[myidx].copy();co,_,myclusters=embedding(b,'myeloid')
 evidence=evidence.join(co)
 clustering=[]
 for subset,labs in [('TME',clusters),('myeloid',myclusters)]:
  for k,v in labs.items():
   for k2,v2 in labs.items():
    if k<k2:clustering.append({'subset':subset,'partition_1':k,'partition_2':k2,'adjusted_rand_index':adjusted_rand_score(v,v2)})
 pd.DataFrame(clustering).to_csv(OUT/'tables/clustering_stability.csv',index=False)
 # Descriptive cluster markers; cell-wise p-values are not donor-level inference.
 log=np.log1p(x*(1e4/x.sum(1))[:,None]);markerrows=[]
 for subset,idx,lab in [('TME',np.arange(len(a)),clusters['TME_louvain_0.8']),('myeloid',myidx,myclusters['myeloid_louvain_0.8'])]:
  for group in np.unique(lab):
   ii=idx[lab==group];jj=idx[lab!=group]
   if len(ii)<3 or len(jj)<3:continue
   use=np.flatnonzero((x[ii]>0).mean(0)>=.1)
   test=mannwhitneyu(log[ii][:,use],log[jj][:,use],axis=0,alternative='two-sided',method='asymptotic')
   q=bh(test.pvalue);delta=log[ii][:,use].mean(0)-log[jj][:,use].mean(0)
   pick=np.argsort(delta)[::-1][:30]
   for z in pick:markerrows.append({'subset':subset,'cluster':int(group),'n_cells':len(ii),'gene':genes[use[z]],'mean_log1p_difference':delta[z],'fraction_in':(x[ii,use[z]]>0).mean(),'fraction_out':(x[jj,use[z]]>0).mean(),'p_descriptive':test.pvalue[z],'q_within_cluster':q[z]})
 pd.DataFrame(markerrows).to_csv(OUT/'tables/descriptive_cluster_markers.csv',index=False)
 # Leave-one-marker-out calibration changes the score AND the matched null.
 loo=[]
 for name in ['Astrocyte','Microglia','Macrophage','BAM','Monocyte']:
  for g in PANELS[name]:
   r=cal.score([h for h in PANELS[name] if h!=g],name+'_minus_'+g);r['omitted_gene']=g
   lab,_,_=classify(2,overrides={name:r})
   for cell,base,pert in zip(evidence.index,evidence.tme_annotation_revised,lab):loo.append({'Cell':cell,'perturbed_panel':name,'omitted_gene':g,'primary_annotation':base,'perturbed_annotation':pert,'unchanged':base==pert})
 pd.DataFrame(loo).to_csv(OUT/'tables/leave_one_marker_out.csv',index=False)
 stability=pd.DataFrame(loo).groupby('Cell').unchanged.mean();evidence['leave_one_marker_out_stability']=stability
 evidence.to_csv(OUT/'tables/per_cell_annotation_evidence.csv')
 smallcols=['ID','subtype','Final_Annotation','Final_Annotation_Focus','tme_annotation_revised','annotation_confidence','annotation_reason','qc_flag','included_in_conservative_astrocytes','included_in_conservative_microglia','astrocyte_nonmetabolic_supported','C2S_original_prediction','annotation_k1','annotation_k3','leave_one_marker_out_stability']
 evidence[smallcols].to_csv(OUT/'tables/TME_revised_annotations.csv')
 a.obs=evidence.copy();a.write_h5ad(OUT/'data/TME_annotated.h5ad',compression='gzip')
 full=pd.read_csv(OUT/'data/full_original_metadata.csv').set_index('Cell',drop=False)
 full['analysis_annotation']=full.Final_Annotation.astype(str);full.loc[evidence.index,'analysis_annotation']=evidence.tme_annotation_revised
 assert full.loc[full.Final_Annotation_Focus.ne('TME'),'analysis_annotation'].equals(full.loc[full.Final_Annotation_Focus.ne('TME'),'Final_Annotation'])
 full.to_csv(OUT/'tables/full_cohort_labels_preserved.csv',index=False)
 pd.crosstab(evidence.Final_Annotation,evidence.tme_annotation_revised).to_csv(OUT/'tables/annotation_transitions.csv')
 evidence.groupby(['ID','subtype','tme_annotation_revised']).size().rename('n').reset_index().to_csv(OUT/'tables/revised_counts_by_donor.csv',index=False)
 evidence.groupby(['subtype','tme_annotation_revised']).size().rename('n').reset_index().to_csv(OUT/'tables/revised_counts_by_subtype.csv',index=False)
 counts=evidence.tme_annotation_revised.value_counts().to_dict()
 versions={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','scikit-learn','scanpy','anndata','networkx','umap-learn']}
 summary={'counts':counts,'conservative_astrocytes':int(evidence.included_in_conservative_astrocytes.sum()),'conservative_microglia':int(evidence.included_in_conservative_microglia.sum()),'additional_qc_exclusions':int(evidence.qc_exclusion.sum()),'reference_comparison':'Published signature enrichment, not whole-transcriptome reference mapping or independent validation','versions':versions,'malignant_labels_changed':0,'TME_n':516,'seed':C['seed']}
 (OUT/'tables/annotation_summary.json').write_text(json.dumps(summary,indent=2))
 print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
