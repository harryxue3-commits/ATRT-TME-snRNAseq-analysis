#!/usr/bin/env python3
"""Extract the frozen TME and fixed-label tumor summaries from deposited counts."""
from pathlib import Path
import csv, hashlib, importlib.util, json, re, sys
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

OUT = Path(__file__).resolve().parents[1]
ROOT = OUT.parents[2]
META = ROOT/'NCBIGEO/snRNA/GSE283839_ATRT_RNA_V3_metadata.tsv'
COUNTS = ROOT/'NCBIGEO/snRNA/GSE283839_ATRT_RNA_V3_counts_raw.tsv'

PANELS = {
 'Astrocyte':['ALDH1L1','AQP4','SLC1A2','SLC1A3','GFAP','SOX9'],
 'Astrocyte_no_metabolic':['ALDH1L1','AQP4','GFAP','SOX9','S100B','FGFR3'],
 'Microglia':['P2RY12','TMEM119','GPR34','SALL1','P2RY13','SLC2A5','OLFML3','CX3CR1'],
 'Macrophage':['CD163','MRC1','MS4A7','CD14','LYZ','FCGR3A','ITGA4','F13A1','STAB1'],
 'Monocyte':['FCN1','VCAN','S100A8','S100A9','CCR2'],
 'BAM':['MRC1','LYVE1','F13A1','CD163','STAB1'],
 'Myeloid':['LST1','TYROBP','FCER1G','AIF1','CSF1R','SPI1','CTSS','CD68','C1QA','C1QB','C1QC'],
 'OPC':['PDGFRA','OLIG1','OLIG2','SOX10','CSPG4'],
 'Oligodendrocyte':['MBP','PLP1','MOG','MOBP'],
 'Neuron':['SNAP25','SYT1','RBFOX3','STMN2'],
 'Endothelial':['CLDN5','PECAM1','VWF','KDR'],
 'Pericyte':['RGS5','PDGFRB','ABCC9','KCNJ8','CSPG4'],
 'T_cell':['CD3D','CD3E','TRAC'],
 'NK':['NKG7','GNLY','KLRD1'],
 'B_cell':['MS4A1','CD79A','CD79B'],
 'Dendritic':['CD1C','FCER1A','CLEC10A','CLEC9A','XCR1'],
}
TUMOR_STATES=['CP-like','Cilia-like','Hypoxic','IPC-like','Mesenchymal-like','NPC-like','OPC-like','RG-like']
LR = [
 ('PTN','PTPRZ1','single receptor'),('MDK','PTPRZ1','single receptor'),
 ('NCAM1','FGFR1','single receptor'),('SPP1','CD44','single receptor'),
 ('SPP1','ITGAV+ITGB1','integrin heterodimer'),('SPP1','ITGAV+ITGB3','integrin heterodimer'),
 ('SPP1','ITGAV+ITGB5','integrin heterodimer'),('SPP1','ITGA4+ITGB1','integrin heterodimer'),
 ('SPP1','ITGA5+ITGB1','integrin heterodimer'),('SPP1','ITGA9+ITGB1','integrin heterodimer'),
 ('APP','CD74','candidate CD74 association; no activation inferred'),
 ('MIF','CD74+CD44','CD74 co-receptor complex'),('MIF','CD74+CXCR4','CD74 co-receptor complex'),
 ('TGFB1','TGFBR1+TGFBR2','heteromeric receptor'),('TGFB2','TGFBR1+TGFBR2','heteromeric receptor'),
 ('TGFB3','TGFBR1+TGFBR2','heteromeric receptor'),('GAS6','AXL','single receptor'),
 ('GAS6','MERTK','single receptor'),('IGF1','IGF1R','single receptor'),
 ('PDGFB','PDGFRB','single receptor gene; dimerization not measured'),
 ('PDGFC','PDGFRA','single receptor gene; dimerization not measured'),
 ('LAMA2','CD44','candidate association'),('MDK','LRP1','single receptor'),
 ('JAG1','NOTCH1','single receptor'),('SEMA3C','PLXND1','candidate pair; accessory receptors not resolved'),
 ('THBS1','CD47','single receptor'),('TNC','ITGAV+ITGB1','integrin heterodimer'),
]

def load_module(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
 for d in ['data','tables','figures','reference','logs','cache']:(OUT/d).mkdir(exist_ok=True)
 meta=pd.read_csv(META,sep='\t').set_index('Cell',drop=False)
 assert len(meta)==17564 and meta.index.is_unique
 tme=meta.loc[meta.Final_Annotation_Focus.eq('TME')].copy(); assert len(tme)==516
 support=load_module(ROOT/'Supplementary/Code/03_astrocyte_support_analysis/analyze_astrocyte_support_modules.py','support')
 astro=load_module(ROOT/'Supplementary/Code/02_astrocyte_state_analysis/score_astrocyte_state_programs.py','astro')
 config={'seed':20261005,'n_null':2000,'null_pool_size':200,'max_rank':1500,'panels':PANELS,'tumor_states':TUMOR_STATES,'support_modules':support.MODULES,'astro_programs':astro.SIGNATURES,'lr_pairs':[dict(ligand=l,receptor=r,note=n) for l,r,n in LR],
 'rules':{'identity_empirical_p':0.05,'minimum_marker_counts':[1,2,3],'primary_minimum_markers':2,'competing_enrichment_z_margin':1.0,'minimum_myeloid_backbone_markers':2,'astro_anchor':['ALDH1L1','AQP4','SLC1A2','SLC1A3'],'astro_nonmetabolic_anchor':['ALDH1L1','AQP4'],'microglia_anchor':['P2RY12','TMEM119','GPR34','SALL1'],'minimum_neighborhood_support':0.5,'qc_exclusion_min_genes':200,'qc_exclusion_max_mito_percent':10,'min_astro_pseudobulk':5,'min_tumor_pseudobulk':20,'communication_group_sizes':[5,10,20]}}
 (OUT/'code/analysis_config.json').write_text(json.dumps(config,indent=2))
 ref=pd.ExcelFile(OUT/'reference/Sankowski_2024_supplementary_tables.xlsx')
 signatures=pd.read_excel(ref,sheet_name='Suppl. Table 21',header=1)
 refsets={c:[str(x).strip() for x in signatures[c].dropna() if str(x).strip()] for c in signatures}
 dev=pd.read_excel(ref,sheet_name='Suppl. Table 13',header=1)
 refrows=[]
 for c,gs in refsets.items():
  for g in gs:refrows.append(dict(reference_panel=c,gene=g,source_table=21))
 for group,df in dev.groupby('cluster'):
  filt=df[(df.p_val_adj<.05)&(df.avg_log2FC>0)&(df['pct.1']>df['pct.2'])].sort_values('avg_log2FC',ascending=False)
  gs=[g for g in filt.gene.astype(str) if g not in ['XIST','MALAT1','NEAT1'] and not g.startswith(('MT-','RPL','RPS'))][:30]
  refsets['developmental_'+str(group)]=gs
  for g in gs:refrows.append(dict(reference_panel='developmental_'+str(group),gene=g,source_table=13))
 (OUT/'reference/reference_gene_sets.json').write_text(json.dumps(refsets,indent=2))
 pd.DataFrame(refrows).to_csv(OUT/'tables/reference_gene_sets.csv',index=False)
 wanted=set(g for gs in support.MODULES.values() for g in gs)|set(g for l,r,n in LR for g in [l]+r.split('+'))
 wanted|={'CD74','APP','SPP1','APOE','APOC1','GPNMB','CD163','ITGAV','ITGB1','CXCR4','CD44'}
 genes=[]; tmecount=[]; bulk=[]; selected=[]; lrgenes=[]
 groupnames=sorted((meta.ID+'|'+meta.Final_Annotation).unique()); groupmap={x:i for i,x in enumerate(groupnames)}
 sha=hashlib.sha256()
 with COUNTS.open('rb') as f:
  line=f.readline();sha.update(line);barcodes=[re.sub(r'^([ACGT]+)\.([0-9]+)',r'\1-\2',x) for x in line.decode().rstrip('\n\r').split('\t')]
  assert len(barcodes)==len(meta) and len(set(barcodes))==len(meta) and set(barcodes)==set(meta.index)
  meta=meta.loc[barcodes];positions={c:i for i,c in enumerate(barcodes)}
  tidx=np.array([positions[c] for c in tme.index]); groupidx=np.array([groupmap[x] for x in meta.ID+'|'+meta.Final_Annotation])
  lib=np.zeros(len(meta),dtype=np.int64);gene_n=np.zeros(len(meta),dtype=np.int32);mt=np.zeros(len(meta),dtype=np.int64)
  for n,line in enumerate(f):
   sha.update(line);g,values=line.decode().rstrip('\n\r').split('\t',1);v=np.fromstring(values,sep='\t',dtype=np.int64)
   assert len(v)==len(meta) and (v>=0).all(),(n,g)
   genes.append(g);tmecount.append(v[tidx].astype(np.int32));lib+=v;gene_n+=(v>0)
   if g.startswith('MT-'):mt+=v
   bulk.append(np.bincount(groupidx,weights=v,minlength=len(groupnames)).astype(np.int64))
   if g in wanted:lrgenes.append(g);selected.append(v.astype(np.int32))
   if n%10000==0:print('Read genes',n,flush=True)
 assert len(set(genes))==len(genes),'Duplicate genes require aggregation'
 x=np.stack(tmecount).T
 assert np.array_equal(x.sum(1),lib[tidx])
 tme['raw_library_size']=lib[tidx];tme['raw_n_genes']=gene_n[tidx];tme['raw_percent_mt']=mt[tidx]/lib[tidx]*100
 obj=ad.AnnData(sparse.csr_matrix(x),obs=tme.drop(columns='Cell'),var=pd.DataFrame(index=pd.Index(genes,name='gene')))
 obj.layers['counts']=obj.X.copy();obj.uns['scope']='Final_Annotation_Focus == TME; original malignant assignments fixed'
 obj.write_h5ad(OUT/'data/TME_raw_counts.h5ad',compression='gzip')
 np.savez_compressed(OUT/'data/fixed_tumor_summaries.npz',counts=np.stack(bulk).T,groupnames=np.array(groupnames),genes=np.array(genes),selected_counts=np.stack(selected).T,selected_genes=np.array(lrgenes),barcodes=np.array(barcodes),library_sizes=lib)
 tme.to_csv(OUT/'tables/TME_original_metadata.csv',index=False)
 meta.to_csv(OUT/'data/full_original_metadata.csv',index=False)
 manifest={'metadata_sha256':hashlib.sha256(META.read_bytes()).hexdigest(),'counts_sha256':sha.hexdigest(),'reference_sha256':hashlib.sha256((OUT/'reference/Sankowski_2024_supplementary_tables.xlsx').read_bytes()).hexdigest(),'reference_url':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41591-023-02673-1/MediaObjects/41591_2023_2673_MOESM2_ESM.xlsx','n_all':len(meta),'n_tme':len(tme),'n_non_tme_unchanged':len(meta)-len(tme),'n_genes':len(genes),'barcode_match':'one-to-one; first barcode numeric suffix only','raw_counts_integer':True}
 (OUT/'data/input_manifest.json').write_text(json.dumps(manifest,indent=2))
 print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__': main()
