"""Preserve supplied annotations; extract integer total RNA counts and auditable CNV metrics."""
from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

OUT=Path(__file__).resolve().parents[1]
OLD=OUT.parents[1]/'Cancers_ATRT_TME_Submission/Revision_20261005/analysis'
C=json.loads((OLD/'code/analysis_config.json').read_text())
for src,dst in [(OLD/'code/analysis_config.json',OUT/'reference/previous_analysis_config.json'),(OLD/'reference/reference_gene_sets.json',OUT/'reference/reference_gene_sets.json')]:shutil.copy2(src,dst)
panels=C['panels']; panels['Immune_anchor']=['PTPRC','SPI1','CSF1R'];panels['Tumor_support_only']=['SOX2','NES','MKI67','TOP2A','SMARCB1','SMARCA4']
obs=[];rows=[];gene_map=[]
for path in sorted((OUT.parent/'SCPCP000010_single-cell').rglob('*processed_rna.h5ad')):
 a=ad.read_h5ad(path);o=a.obs.copy();lib=str(o.library_id.iloc[0]);o.index=lib+'|'+o.index.astype(str);o.index.name='cell_id'
 o['cell_barcode']=a.obs_names.astype(str);o['patient_id']=o.participant_id.astype(str);o['patient_group']=o.disease_timing.astype(str);o['molecular_subgroup']='Not provided'
 o['source_cluster']=o.cluster.astype(str);o['source_cluster_id']=lib+'|'+o.source_cluster
 x=a.raw.X.tocsr();x.eliminate_zeros();assert np.all(x.data==np.round(x.data)) and x.min()>=0
 assert np.allclose(x.sum(1).A.ravel(),o['sum'])
 assert np.array_equal(a.raw.var_names,a.var_names)
 symbols=a.raw.var.gene_symbol.astype(str).to_numpy();symbols=np.array([s if s not in ['','nan','None','NA'] else g for s,g in zip(symbols,a.raw.var_names)])
 genes=np.unique(symbols);gi={g:i for i,g in enumerate(genes)}
 agg=sparse.csr_matrix((np.ones(len(symbols),dtype=np.int32),(np.arange(len(symbols)),[gi[g] for g in symbols])),shape=(len(symbols),len(genes)))
 x=(x.astype(np.int32)@agg).tocsr();x.eliminate_zeros()
 o['raw_library_size']=x.sum(1).A.ravel();o['raw_n_genes']=(x>0).sum(1).A.ravel()
 o['qc_pass']=o.scDblFinder_class.astype(str).eq('singlet') & (o.raw_n_genes>=200) & (o.subsets_mito_percent<=10) & (o.raw_library_size>0)
 for name,gs in panels.items():
  inds=[gi[g] for g in gs if g in gi];o[name+'__n_detected']=(x[:,inds]>0).sum(1).A.ravel()
 for name,gs in {'astro_anchor':C['rules']['astro_anchor'],'microglia_anchor':C['rules']['microglia_anchor'],'macrophage_anchor':['CD163','MRC1','ITGA4','F13A1','STAB1']}.items():
  o[name+'_detected']=(x[:,[gi[g] for g in gs if g in gi]]>0).sum(1).A.ravel()
 o['normal_reference_marker_supported']=o.is_infercnv_reference.eq(True)&o.qc_pass&(
  ((o['Myeloid__n_detected']>=3)&(o['Immune_anchor__n_detected']>=1))|
  (o['T_cell__n_detected']>=2)|(o['B_cell__n_detected']>=2)|
  (o['Endothelial__n_detected']>=3)|(o['Pericyte__n_detected']>=3))
 o['cnv_mean_arm_fraction']=np.nan;o['cnv_large_arms']=np.nan
 if a.uns.get('infercnv_status')=='success':
  f=a.uns['infercnv_table'].set_index('barcodes').loc[o.cell_barcode];f.index=o.index
  pc=f.filter(regex='^proportion_cnv_chr');o['cnv_mean_arm_fraction']=pc.mean(1);o['cnv_large_arms']=(pc>=.5).sum(1)
  o['cnv_subcluster']=f.subcluster
  assert np.allclose(f.filter(regex='^has_cnv_chr').sum(1),o.infercnv_total_cnv)
  for kind,sel in [('reference_all',o.is_infercnv_reference.eq(True)),('reference_marker_supported',o.normal_reference_marker_supported),('query',o.is_infercnv_reference.eq(False))]:
   for metric in ['infercnv_total_cnv','cnv_mean_arm_fraction','cnv_large_arms']:
    z=o.loc[sel,metric];rows.append({'library_id':lib,'reference_selection':kind,'metric':metric,'n':len(z),'min':z.min(),'q50':z.quantile(.5),'q75':z.quantile(.75),'q90':z.quantile(.9),'q95':z.quantile(.95),'q99':z.quantile(.99),'max':z.max()})
 else:o['cnv_subcluster']='Not available'
 sr=a.uns['singler_results'].reindex(o.cell_barcode);sr.index=o.index
 o['singler_delta_next']=sr['delta.next'];o['singler_pruned']=sr['pruned.labels'].isna()
 # Rank scores are diagnostic; source labels remain verbatim.
 scols=sr.filter(regex='^scores\\.');sel=['scores.'+str(t).replace(':','.') for t in sr.labels]
 o['singler_delta_median']=[float(scols.iloc[j][c]-scols.iloc[j].median()) if c in scols.columns else np.nan for j,c in enumerate(sel)]
 o['source_umap1']=a.obsm['X_umap'][:,0];o['source_umap2']=a.obsm['X_umap'][:,1]
 b=ad.AnnData(x,obs=o,var=pd.DataFrame(index=pd.Index(genes,name='gene_symbol')))
 b.obsm['X_source_umap']=a.obsm['X_umap'];b.obsm['X_source_pca']=a.obsm['X_pca']
 # Neighbors are within sample and restricted to retained singlets; no source labels are used.
 inds=np.flatnonzero(o.qc_pass);nbr=NearestNeighbors(n_neighbors=16,metric='cosine',n_jobs=2).fit(a.obsm['X_pca'][inds,:30]);di,nn=nbr.kneighbors()
 np.savez_compressed(OUT/'data'/f'{lib}_neighbors.npz',cell_index=inds,neighbors=inds[nn[:,:15]],distances=di[:,:15])
 b.write_h5ad(OUT/'data'/f'{lib}_counts.h5ad',compression='gzip')
 if len(gene_map)==0:
  gene_map=pd.DataFrame({'ensembl_id':a.raw.var_names,'gene_symbol_input':a.raw.var.gene_symbol.astype(str).values,'gene_symbol_output':symbols});gene_map.to_csv(OUT/'tables/gene_id_mapping.csv',index=False)
 obs.append(o);print(lib,'prepared',len(o),'QC',int(o.qc_pass.sum()),'clean reference',int(o.normal_reference_marker_supported.sum()),flush=True)
 del a,b,x
pd.concat(obs).to_csv(OUT/'tables/prepared_metadata.csv.gz')
pd.DataFrame(rows).to_csv(OUT/'tables/CNV_distributions.csv',index=False)
print('Preparation complete',flush=True)
