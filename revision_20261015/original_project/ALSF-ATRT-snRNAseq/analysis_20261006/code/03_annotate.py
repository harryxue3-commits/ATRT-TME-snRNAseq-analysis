"""Compartment-aware conservative and broader annotations, with no forced tumor/normal calls."""
from pathlib import Path
import json,warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
OUT=Path(__file__).resolve().parents[1]
C=json.loads((OUT/'reference/previous_analysis_config.json').read_text())
ALTERNATIVES=['OPC','Oligodendrocyte','Neuron','Endothelial','Pericyte','T_cell','NK','B_cell','Dendritic']
NORMAL_TYPES=['Astrocyte','Microglia','Macrophage_like','BAM_like','Monocyte_like','Myeloid_unresolved']+ALTERNATIVES

def supported(e,n,k=2):return (e[n+'__empirical_p']<.05)&(e[n+'__n_detected']>=k)
def neighbors(e,mask):
 out=pd.Series(0.,index=e.index)
 for lib,d in e.groupby('library_id',sort=False):
  f=np.load(OUT/'data'/f'{lib}_neighbors.npz');vals=mask.loc[d.index].to_numpy();out.loc[d.index[f['cell_index']]]=vals[f['neighbors']].mean(1)
 return out

def classify_tme(e,k=2,context=.5):
 def sup(n,kk=2):return supported(e,n,kk)
 def competing(target):
  bad=pd.Series(False,index=e.index)
  for name in ALTERNATIVES+(['Myeloid'] if target=='Astrocyte' else ['Astrocyte']):
   if name!=target:
    positive=sup(name)
    if name=='Astrocyte':positive &= sup('Astrocyte_no_metabolic') & e.astro_nonmetabolic_anchor_detected.ge(1)
    bad|=positive&(e[name+'__enrichment_z']>e[target+'__enrichment_z']+1)
  return bad
 # Without CNV results only marker-supported immune/vascular source-reference cells enter the normal gate.
 gate=e.qc_pass&(e.cnv_reference_like|e.normal_reference_without_cnv)
 # Audit showed SLC1A3-driven apparent astrocytes inside myeloid populations. Require an independent
 # non-transporter panel and ALDH1L1/AQP4, and veto a coherent myeloid lineage regardless of relative z.
 my_lineage_conflict=((e['Myeloid__n_detected']>=2)&((e.myeloid_neighbor_support>=.5)|sup('reference_Microglia')|sup('reference_CAMs')))|e.source_myeloid_context
 astro=sup('Astrocyte',k)&(e.astro_anchor_detected>=1)&sup('Astrocyte_no_metabolic')&e.astro_nonmetabolic_anchor_detected.ge(1)&~my_lineage_conflict&~competing('Astrocyte')&((e.astro_neighbor_support>=context)|sup('reference_Astrocytes'))&e.cnv_reference_like&e.qc_pass
 my=gate&(e['Myeloid__n_detected']>=2)&((e.myeloid_neighbor_support>=context)|sup('reference_Microglia')|sup('reference_CAMs')|sup('reference_Monocytes'))&~competing('Myeloid')
 mg=my&sup('Microglia',k)&(e.microglia_anchor_detected>=1)
 mac=my&sup('Macrophage',max(2,k))&(e.macrophage_anchor_detected>=1)
 bam=my&sup('BAM',max(2,k));mono=my&sup('Monocyte',max(2,k))
 label=pd.Series('Unresolved',index=e.index);label[my]='Myeloid_unresolved'
 label[astro]='Astrocyte';label[mg&~(mac|bam|mono)]='Microglia'
 label[mac&~(mg|bam|mono)]='Macrophage_like';label[bam&~(mg|mono)]='BAM_like';label[mono&~(mg|bam)]='Monocyte_like'
 label[astro&my]='Mixed_lineage_unresolved'
 # Secondary normal classes require markers, normal-compatible CNV and support from at least one source method.
 for p in ALTERNATIVES:
  gs=gate&sup(p)&e['source_support_'+p]&label.eq('Unresolved')
  if p in ['OPC','Oligodendrocyte','Neuron']:gs &= e.cnv_reference_like
  for alt in ALTERNATIVES:
   if alt!=p:gs &= ~(sup(alt)&(e[alt+'__enrichment_z']>e[p+'__enrichment_z']-1))
  gs &= ~(my|astro)
  label[gs]=p
 label[~e.qc_pass]='QC_excluded'
 return label

def source_harmonize(value):
 t=str(value).lower()
 if 'microgli' in t:return 'Microglia'
 if any(x in t for x in ['macrophage','monocyte','myeloid','phagocyte']):return 'Myeloid'
 if 'astrocy' in t:return 'Astrocyte'
 if 'oligodendrocyte progenitor' in t:return 'OPC'
 if 'oligodendro' in t:return 'Oligodendrocyte'
 if any(x in t for x in ['neuron','pyramidal','purkinje','interneuron']) and 'associated' not in t:return 'Neuron'
 if 'endothel' in t:return 'Endothelial'
 if 'pericyte' in t:return 'Pericyte'
 if 'dendritic' in t and 'follicular' not in t:return 'Dendritic'
 if 'natural killer' in t or t=='nk cells':return 'NK'
 if any(x in t for x in ['t cell','t cells','lymphocyte of t']):return 'T_cell'
 if 'b cell' in t or 'plasma cell' in t or 'plasmablast' in t:return 'B_cell'
 if any(x in t for x in ['glial cell','bergmann','radial glia']):return 'Glial_unspecified'
 return 'Other_or_unresolved'

def main():
 e=pd.read_csv(OUT/'tables/per_cell_scores.csv.gz',index_col=0,low_memory=False)
 e['normal_reference_without_cnv']=e.infercnv_total_cnv.isna()&e.normal_reference_marker_supported
 for m in ['singler','cellassign','scimilarity','consensus']:e[m+'_broad_class']=e[m+'_celltype_annotation'].map(source_harmonize)
 # Read exact raw anchor counts, rather than treating a broad astrocyte score as an independent check.
 import anndata as ad
 e['astro_nonmetabolic_anchor_detected']=0
 for lib,d in e.groupby('library_id',sort=False):
  a=ad.read_h5ad(OUT/'data'/f'{lib}_counts.h5ad',backed='r');e.loc[d.index,'astro_nonmetabolic_anchor_detected']=(a[:,['ALDH1L1','AQP4']].X>0).sum(1).A.ravel();a.file.close()
 for p in ALTERNATIVES:e['source_support_'+p]=e[[m+'_broad_class' for m in ['singler','cellassign','scimilarity']]].eq(p).any(axis=1)
 e['myeloid_neighbor_support']=neighbors(e,e['Myeloid__n_detected'].ge(2)&(e.cnv_reference_like|e.normal_reference_without_cnv))
 source_my=e[[m+'_broad_class' for m in ['singler','cellassign','scimilarity']]].isin(['Microglia','Myeloid']).sum(axis=1).ge(2)
 e['source_myeloid_context']=source_my&neighbors(e,source_my).ge(.5)
 e['astro_neighbor_support']=neighbors(e,supported(e,'Astrocyte',2)&supported(e,'Astrocyte_no_metabolic')&e.astro_nonmetabolic_anchor_detected.ge(1)&e.cnv_reference_like&~e.source_myeloid_context)
 for k in [1,2,3]:e[f'tme_marker_k{k}']=classify_tme(e,k=k)
 e['annotation_conservative']=e.tme_marker_k2
 e['annotation_broad']=classify_tme(e,k=1,context=.25)
 # Preserve nesting: changing k can reveal competing programs; do not silently remove a conservative cell from the sensitivity set.
 core=e.annotation_conservative.isin(NORMAL_TYPES)
 e.loc[core,'annotation_broad']=e.loc[core,'annotation_conservative']
 e['cnv_high_neighbor_fraction']=neighbors(e,e.cnv_elevated_q95&e.qc_pass)
 e['cnv_high_broad_neighbor_fraction']=neighbors(e,e.cnv_elevated_q90&e.qc_pass)
 # Strong immune evidence is a veto on CNV-only tumor calls, even when CNV is high.
 immune=((e['Myeloid__n_detected']>=3)&supported(e,'Myeloid')&(e['Immune_anchor__n_detected']>=1))|supported(e,'T_cell')|supported(e,'B_cell')|supported(e,'Dendritic')
 # Endothelial/pericyte expression can occur in ATRT; record conflicts rather than calling these cells normal.
 vascular=(supported(e,'Endothelial',3)&e.source_support_Endothelial)|(supported(e,'Pericyte',3)&e.source_support_Pericyte)
 e['tumor_normal_marker_conflict']=immune|vascular
 clsize=e[e.qc_pass].groupby('source_cluster_id').size();e['source_cluster_n_qc']=e.source_cluster_id.map(clsize).fillna(0).astype(int)
 e['tumor_cnv_supported']=e.qc_pass&e.cnv_elevated_q95&(e.cnv_large_arms>=2)&(e.cnv_high_neighbor_fraction>=.8)&(e.source_cluster_n_qc>=20)&~(immune|vascular)&~core
 e['tumor_cnv_broad']=e.qc_pass&e.cnv_elevated_q90&(e.cnv_high_broad_neighbor_fraction>=.5)&(e.source_cluster_n_qc>=20)&~immune&~e.annotation_broad.isin(NORMAL_TYPES)
 e['tumor_cnv_broad']|=e.tumor_cnv_supported
 e.loc[e.tumor_cnv_supported,'annotation_conservative']='Tumor_CNV_supported'
 e.loc[e.tumor_cnv_broad,'annotation_broad']='Tumor_CNV_supported'
 e['compartment_conservative']=np.select([e.annotation_conservative.isin(NORMAL_TYPES),e.tumor_cnv_supported,~e.qc_pass],['TME_marker_supported','Tumor_CNV_supported','QC_excluded'],default='Unresolved')
 e['annotation_reason']='Insufficient or conflicting evidence; not used as a named CellChat population'
 e.loc[core,'annotation_reason']='Normal-compatible compartment gate plus coherent marker panel and lineage safeguards'
 e.loc[core&e.normal_reference_without_cnv,'annotation_reason']='Marker-supported source normal reference; CNV unavailable (immune/vascular only)'
 e.loc[e.tumor_cnv_supported,'annotation_reason']='CNV count and mean affected-arm fraction above sample normal-reference threshold, >=2 arms >=50% affected, >=80% high-CNV neighbors, no strong immune/vascular conflict'
 e.loc[~e.qc_pass,'annotation_reason']='Predicted doublet, <200 detected genes, >10% mitochondrial reads, or empty library'
 e['patient_group']=e.disease_timing
 e.to_csv(OUT/'tables/cell_annotations.csv.gz');e.to_csv(OUT/'tables/cell_annotations.csv')
 e.groupby(['library_id','sample_id','patient_id','disease_timing','annotation_conservative'],observed=True).size().rename('n_cells').reset_index().to_csv(OUT/'tables/conservative_counts_by_sample.csv',index=False)
 e.groupby(['library_id','sample_id','patient_id','disease_timing','annotation_broad'],observed=True).size().rename('n_cells').reset_index().to_csv(OUT/'tables/broad_counts_by_sample.csv',index=False)
 print('Conservative',e.annotation_conservative.value_counts().to_dict(),flush=True);print('Broad',e.annotation_broad.value_counts().to_dict(),flush=True)
 assert e.index.is_unique and len(e)==42703
 assert not (e.tumor_cnv_supported&e.annotation_conservative.isin(NORMAL_TYPES)).any()
 assert e.loc[e.annotation_conservative.isin(NORMAL_TYPES+['Tumor_CNV_supported']),'qc_pass'].all()
if __name__=='__main__':main()
