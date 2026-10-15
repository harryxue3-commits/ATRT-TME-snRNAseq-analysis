from pathlib import Path
import json,gzip,hashlib,importlib.util,os,warnings
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/atrt_mpl')
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
from threadpoolctl import threadpool_limits
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
OUT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('annot',OUT/'code/03_annotate.py');A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
e=pd.read_csv(OUT/'tables/cell_annotations.csv.gz',index_col=0,low_memory=False)
extra=pd.read_csv(OUT/'tables/CNV_unavailable_expression_candidates.csv',index_col=0)
e['tumor_expression_only_candidate']=e.index.map(extra.exploratory_tumor_candidate).fillna(False).astype(bool)
e['expression_tumor_score']=e.index.map(extra.expression_classifier_tumor_score)
valid=A.NORMAL_TYPES+['Tumor_CNV_supported']
for tier in ['conservative','broad']:
 e['cellchat_'+tier]=e['annotation_'+tier].where(e['annotation_'+tier].isin(valid),'Exclude')
e['cellchat_broad_macrophage_BAM_merged']=e.cellchat_broad.replace({'Macrophage_like':'Macrophage_BAM_like','BAM_like':'Macrophage_BAM_like'})
e['cellchat_exploratory_expression_tumor']=e.cellchat_broad.copy();e.loc[e.tumor_expression_only_candidate,'cellchat_exploratory_expression_tumor']='Tumor_expression_only'
e['conservative_astrocyte']=e.annotation_conservative.eq('Astrocyte');e['conservative_microglia']=e.annotation_conservative.eq('Microglia');e['conservative_macrophage_like']=e.annotation_conservative.eq('Macrophage_like');e['conservative_BAM_like']=e.annotation_conservative.eq('BAM_like')
e['broader_microglia_myeloid_uncertain']=e.annotation_broad.isin(['Microglia','Myeloid_unresolved'])
e['broader_macrophage_BAM_like']=e.annotation_broad.isin(['Macrophage_like','BAM_like'])
cols=['cell_barcode','library_id','sample_id','patient_id','submitter_id','disease_timing','patient_group','molecular_subgroup','annotation_conservative','annotation_broad','compartment_conservative','cellchat_conservative','cellchat_broad','cellchat_broad_macrophage_BAM_merged','cellchat_exploratory_expression_tumor','qc_pass','scDblFinder_class','infercnv_total_cnv','cnv_mean_arm_fraction','is_infercnv_reference','normal_reference_without_cnv','tumor_cnv_supported','tumor_cnv_broad','tumor_expression_only_candidate','tumor_passes_q99_sensitivity','conservative_astrocyte','conservative_microglia','conservative_macrophage_like','conservative_BAM_like','broader_microglia_myeloid_uncertain','broader_macrophage_BAM_like','LOO_Microglia_stability','LOO_Macrophage_stability','LOO_BAM_stability','singler_celltype_annotation','cellassign_celltype_annotation','scimilarity_celltype_annotation','consensus_celltype_annotation','annotation_reason']
registry=e[cols].copy();registry.index.name='cell_id';registry.to_csv(OUT/'ALSF_ATRT_cell_registry.csv')
e.to_csv(OUT/'tables/cell_annotations.csv.gz');e.to_csv(OUT/'tables/cell_annotations.csv')
defs=[]
for tier in ['conservative','broad']:
 for typ in valid:
  d=registry[e['annotation_'+tier].eq(typ)];cid=f'{tier}_{typ}';d.to_csv(OUT/'cell_sets'/f'{cid}.csv')
  defs.append({'cohort_id':cid,'tier':tier,'population':typ,'n_cells':len(d),'n_patients':d.patient_id.nunique(),'n_samples':d.sample_id.nunique(),'definition':'Disjoint annotation within tier; normal compartment gate plus marker criteria, or CNV-supported tumor criteria.'})
for col in ['broader_microglia_myeloid_uncertain','broader_macrophage_BAM_like','tumor_expression_only_candidate']:
 d=registry[e[col]];d.to_csv(OUT/'cell_sets'/f'{col}.csv');defs.append({'cohort_id':col,'tier':'additional_sensitivity','population':col,'n_cells':len(d),'n_patients':d.patient_id.nunique(),'n_samples':d.sample_id.nunique(),'definition':{'broader_microglia_myeloid_uncertain':'Broad Microglia plus lineage-supported Myeloid_unresolved; not a pure microglia population.','broader_macrophage_BAM_like':'Broad Macrophage_like plus BAM_like; does not establish infiltration.','tumor_expression_only_candidate':'Exploratory expression-only candidates from the CNV-unavailable sample; excluded from both CNV-supported tumor tiers.'}[col]})
pd.DataFrame(defs).to_csv(OUT/'tables/cohort_definitions.csv',index=False)
for panel,label in [('Microglia','Microglia'),('Macrophage','Macrophage_like'),('BAM','BAM_like')]:
 registry[e.annotation_conservative.eq(label)&e['LOO_'+panel+'_stability'].ge(.8)].to_csv(OUT/'cell_sets'/f'optional_marker_stable_{label}.csv')

# All sample/type combinations include explicit zero counts.
counts=[]
for lib,d in e.groupby('library_id',sort=False):
 for tier in ['conservative','broad']:
  for typ in valid:
   n=int(d['annotation_'+tier].eq(typ).sum());counts.append({'library_id':lib,'sample_id':d.sample_id.iloc[0],'patient_id':d.patient_id.iloc[0],'disease_timing':d.disease_timing.iloc[0],'tier':tier,'population':typ,'n_cells':n,'eligible_min10':n>=10,'eligible_min20':n>=20})
pd.DataFrame(counts).to_csv(OUT/'tables/population_counts_and_CellChat_eligibility.csv',index=False)

# Tool comparisons are agreement/conflict diagnostics, not estimates of accuracy against independent truth.
methods=[];cross=[]
for m in ['singler','cellassign','scimilarity','consensus']:
 for subset,mask in [('all',np.ones(len(e),bool)),('conservative_TME',e.compartment_conservative.eq('TME_marker_supported')),('CNV_supported_tumor',e.tumor_cnv_supported)]:
  d=e.loc[mask];lab=d[m+'_celltype_annotation'].fillna('Missing');quality=np.nan
  if m=='scimilarity':quality=d.scimilarity_min_distance.gt(.05).mean()
  if m=='cellassign':quality=d.cellassign_max_prediction.ge(.99).mean()
  methods.append({'method':m,'subset':subset,'n_cells':len(d),'top_label':lab.value_counts().index[0],'top_label_count':int(lab.value_counts().iloc[0]),'unknown_or_other_fraction':lab.str.lower().isin(['unknown','other','missing','unclassified cell']).mean(),'scimilarity_distance_gt005_fraction':quality if m=='scimilarity' else np.nan,'cellassign_p_ge099_fraction':quality if m=='cellassign' else np.nan})
  pd.crosstab(d.annotation_conservative,lab).to_csv(OUT/'tables'/f'{m}_{subset}_crosswalk.csv')
  if subset=='conservative_TME':
   expected=d.annotation_conservative.replace({'Microglia':'Myeloid','Macrophage_like':'Myeloid','BAM_like':'Myeloid','Monocyte_like':'Myeloid','Myeloid_unresolved':'Myeloid'})
   observed=d[m+'_broad_class'].replace({'Microglia':'Myeloid'})
   for t in expected.unique():
    s=expected.eq(t);cross.append({'method':m,'marker_defined_lineage':t,'n_cells':int(s.sum()),'n_agree':int((observed[s]==t).sum()),'agreement_fraction':float((observed[s]==t).mean()),'interpretation':'Marker-defined operational labels are not independent ground truth.'})
pd.DataFrame(methods).to_csv(OUT/'tables/method_diagnostics.csv',index=False);pd.DataFrame(cross).to_csv(OUT/'tables/method_lineage_agreement.csv',index=False)
astroaudit=[]
for m in ['singler','cellassign','scimilarity']:
 d=e[e[m+'_broad_class'].eq('Astrocyte')]
 astroaudit.append({'method':m,'n_astrocyte_calls':len(d),'n_qc_pass':int(d.qc_pass.sum()),'n_CNV_elevated_q95':int(d.cnv_elevated_q95.sum()),'n_normal_CNV_compatible':int(d.cnv_reference_like.sum()),'n_primary_tumor':int(d.tumor_cnv_supported.sum()),'n_CNV_unavailable':int(d.infercnv_total_cnv.isna().sum())})
pd.DataFrame(astroaudit).to_csv(OUT/'tables/source_astrocyte_call_audit.csv',index=False)

C=json.loads((OUT/'reference/previous_analysis_config.json').read_text())
citations={'Astrocyte':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4707064/','Astrocyte_no_metabolic':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4707064/','Microglia':'https://pmc.ncbi.nlm.nih.gov/articles/PMC8113099/','Macrophage':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5738907/','BAM':'https://www.nature.com/articles/s41591-023-02673-1','Monocyte':'https://www.nature.com/articles/s41591-023-02673-1'}
rationale={'Astrocyte':'Astrocyte identity panel; SLC1A3 alone is insufficient, particularly in myeloid nuclei.','Astrocyte_no_metabolic':'Independent check avoiding SLC1A2/SLC1A3 transporter-driven calls; requires ALDH1L1 or AQP4.','Microglia':'Resident microglial identity program; requires P2RY12, TMEM119, GPR34 or SALL1. Loss of expression does not prove macrophage identity.','Macrophage':'Macrophage-associated program; requires a positive macrophage anchor, myeloid lineage and no competing microglial/monocyte/BAM call. Infiltration is unproven.','BAM':'Border-associated/CNS-associated macrophage-like program, kept separate. Expression does not prove anatomical location.','Monocyte':'Positive monocyte program; not inferred merely from absence of microglial markers.','Myeloid':'Shared myeloid backbone; distinguishes lineage from activation programs.'}
marker=[]
for p,gs in C['panels'].items():marker.append({'population_or_role':p,'markers':', '.join(gs),'minimum_detected_primary':2,'rationale':rationale.get(p,'Competing normal lineage control; normal-compatible compartment evidence is required.'),'source':citations.get(p,'Prior GSE283839 marker-panel configuration; see reference/previous_analysis_config.json')})
marker.append({'population_or_role':'Tumor_CNV_supported','markers':'No single tumor marker used','minimum_detected_primary':'Not applicable','rationale':'Sample-relative RNA-inferred CNV burden plus neighborhood coherence and normal-lineage conflict checks. SOX2, NES, MKI67 or low SMARCB1 expression alone do not define tumor.','source':'https://scpca.readthedocs.io/en/stable/processing_information.html#cnv-inference'})
pd.DataFrame(marker).to_csv(OUT/'ATRT_marker_table.csv',index=False)

pb=[];pbobs=[];dot=[];tops=[];manifest=[]
(OUT/'data/annotated').mkdir(exist_ok=True)
selected_genes=list(dict.fromkeys(g for p in ['Astrocyte','Microglia','Macrophage','BAM','Myeloid','T_cell','Endothelial'] for g in C['panels'][p]))
all_identity_genes=set(g for gs in C['panels'].values() for g in gs)
all_identity_genes.update(g for gs in json.loads((OUT/'reference/reference_gene_sets.json').read_text()).values() for g in gs)
all_identity_genes.update(['PTPRC'])
for lib,d in e.groupby('library_id',sort=False):
 a=ad.read_h5ad(OUT/'data'/f'{lib}_counts.h5ad');assert a.obs_names.equals(d.index)
 # AnnData requires uniform object/string types; numeric missing values remain numeric NaN.
 o=d.copy()
 for col in o.select_dtypes(include='object').columns:o[col]=o[col].fillna('Not available').astype(str)
 a.obs=o;a.uns['analysis_scope']='All supplied nuclei audited; normal lineage calls restricted to operational TME gate. Raw integer total RNA counts in X.'
 a.uns['normalization_for_downstream']='Use total RNA counts; library size 10000 and natural log1p. Never use scaled/integrated residuals for CellChat.'
 a.uns['source_path']=str(next((OUT.parent/'SCPCP000010_single-cell').rglob(f'{lib}_processed_rna.h5ad')).resolve())
 a.write_h5ad(OUT/'data/annotated'/f'{lib}_annotated.h5ad',compression='gzip')
 dest=OUT/'cellchat'/lib;dest.mkdir(exist_ok=True)
 use=d.cellchat_conservative.ne('Exclude')|d.cellchat_broad.ne('Exclude')|d.tumor_expression_only_candidate
 x=a.X[use.to_numpy()];ids=d.index[use]
 with gzip.open(dest/'counts.mtx.gz','wb',compresslevel=3) as h:
  with threadpool_limits(limits=2):io.mmwrite(h,x.T,field='integer')
 pd.Series(a.var_names).to_csv(dest/'genes.tsv',sep='\t',header=False,index=False);pd.Series(ids).to_csv(dest/'barcodes.tsv',sep='\t',header=False,index=False)
 registry.loc[ids].to_csv(dest/'metadata.csv')
 # Exact export dimension and count invariants; file read-back is performed separately.
 manifest.append({'library_id':lib,'n_export_cells':int(use.sum()),'n_features':a.n_vars,'total_UMI_export':int(x.sum()),'has_CNV':bool(d.infercnv_total_cnv.notna().any())})
 log=a.X.astype(np.float32).multiply(1e4/np.maximum(a.X.sum(1),1)).tocsr();log.data=np.log1p(log.data)
 for tier in ['conservative','broad']:
  for pop in valid:
   mask=d['annotation_'+tier].eq(pop).to_numpy();n=int(mask.sum())
   if not n:continue
   z=a.X[mask];pb.append(sparse.csr_matrix(z.sum(0)));pbobs.append({'pseudobulk_id':f'{lib}|{tier}|{pop}','library_id':lib,'patient_id':d.patient_id.iloc[0],'disease_timing':d.disease_timing.iloc[0],'tier':tier,'population':pop,'n_cells':n})
   if tier=='conservative':
    for gene in selected_genes:
     i=a.var_names.get_loc(gene);dot.append({'library_id':lib,'population':pop,'n_cells':n,'gene':gene,'n_detected':int((z[:,i]>0).sum()),'mean_log1p_CP10k':float(log[mask,i].mean())})
    other=d.qc_pass.to_numpy()&~mask
    if n>=10 and other.sum()>=20:
     mu=log[mask].mean(0).A.ravel();muout=log[other].mean(0).A.ravel();det=(z>0).mean(0).A.ravel();diff=mu-muout;validg=np.flatnonzero(det>=.1);pick=validg[np.argsort(diff[validg])[-30:][::-1]]
     for i in pick:tops.append({'library_id':lib,'patient_id':d.patient_id.iloc[0],'population':pop,'n_cells':n,'gene':a.var_names[i],'mean_log1p_difference_vs_other_QC_nuclei':float(diff[i]),'fraction_detected_in':float(det[i]),'used_in_annotation_panel':a.var_names[i] in all_identity_genes,'inference':'Descriptive only; cells are not independent patients and identity genes are selected by annotation.'})
 print(lib,'exported',int(use.sum()),'CellChat nuclei',flush=True)
pd.DataFrame(manifest).to_csv(OUT/'cellchat/input_manifest.csv',index=False)
pd.DataFrame(dot).to_csv(OUT/'tables/marker_expression_by_sample_population.csv',index=False)
pd.DataFrame(tops).to_csv(OUT/'tables/descriptive_population_markers.csv',index=False)
bm=pd.DataFrame(pbobs).set_index('pseudobulk_id');b=ad.AnnData(sparse.vstack(pb).tocsr().astype(np.int64),obs=bm,var=a.var.copy());b.uns['note']='Separate conservative and broad pseudobulks overlap. Do not treat them as independent replicates. Use patient as biological replicate.';b.write_h5ad(OUT/'data/population_pseudobulk_counts.h5ad',compression='gzip')
bm.to_csv(OUT/'tables/pseudobulk_metadata.csv')
stats={'n_nuclei':len(e),'n_patients':e.patient_id.nunique(),'n_samples':e.sample_id.nunique(),'qc_excluded':int((~e.qc_pass).sum()),'predicted_doublets':int(e.scDblFinder_class.eq('doublet').sum()),'conservative':e.annotation_conservative.value_counts().to_dict(),'broad':e.annotation_broad.value_counts().to_dict(),'source_consensus_unknown':int(e.consensus_celltype_annotation.eq('Unknown').sum()),'SCimilarity_distance_above_005':int(e.scimilarity_min_distance.gt(.05).sum()),'expression_only_tumor_candidates':int(e.tumor_expression_only_candidate.sum())}
(OUT/'tables/final_summary.json').write_text(json.dumps(stats,indent=2))
print(json.dumps(stats,indent=2),flush=True)
