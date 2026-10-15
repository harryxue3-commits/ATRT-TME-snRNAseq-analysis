"""External-cohort expression checks; paired libraries are not independent patients."""
from pathlib import Path
import json,itertools,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parents[1];ROOT=O.parents[2];A=ROOT/'ALSF-ATRT-snRNAseq/analysis_20261006';R=ROOT/'Cancers_ATRT_TME_Submission/Revision_20261007/analysis';T=O/'tables';F=O/'figures'
meta=pd.read_csv(A/'tables/cell_annotations.csv.gz').set_index('cell_id');assert len(meta)==42703 and meta.index.is_unique
use=['library_id','sample_id','patient_id','patient_group','molecular_subgroup','annotation_conservative','annotation_broad','qc_pass','tumor_cnv_supported','tumor_cnv_broad','annotation_reason','cellchat_conservative','cellchat_broad']
meta[use].to_csv(T/'S18_ALSF_cell_ID_groups.csv.gz')
clinical=pd.read_csv(ROOT/'ALSF-ATRT-snRNAseq/SCPCP000010_single-cell/single-cell_metadata.tsv',sep='\t')
clinical.to_csv(T/'S14_ALSF_clinical_metadata_as_provided.csv',index=False)
summary=meta.groupby(['library_id','sample_id','patient_id','patient_group','annotation_conservative'],dropna=False).size().reset_index(name='n_cells');summary.to_csv(T/'S18_ALSF_counts_by_sample.csv',index=False)
for f in ['method_lineage_agreement.csv','method_diagnostics.csv','marker_sensitivity_summary.csv','tumor_CNV_threshold_sensitivity.csv','CNV_operational_thresholds.csv','minimum_marker_sensitivity.csv','source_astrocyte_call_audit.csv','conservative_counts_by_sample.csv','broad_counts_by_sample.csv']:
 shutil.copy2(A/'tables'/f,T/('S18_'+f))
defs=pd.read_csv(R/'tables/S10_microglia_gene_filter_audit.csv');defs=defs[defs.retained];modules={p:list(d.gene[d.gene.ne('CD74')]) for p,d in defs.groupby('module',sort=False)}
genes=list(dict.fromkeys(['APP','CD74','SPP1','CD44','ITGAV','ITGB1','NCAM1','FGFR1','PTN','PTPRZ1','ALDH1L1','AQP4','SLC1A2','SLC1A3','P2RY12','TMEM119','GPR34','SALL1','CD163','MRC1','MS4A7','LYVE1','F13A1']+[g for gs in modules.values() for g in gs]))
rows=[];percell=[];availability=[]
for lib,m in meta.groupby('library_id',sort=True):
 p=A/'data'/f'{lib}_counts.h5ad';a=ad.read_h5ad(p,backed='r');ix=a.var_names.get_indexer(genes);assert (ix>=0).all()
 x=a[:,ix].X.toarray();assert len(x)==len(m)
 obs=a.obs.copy();obs.index=lib+'|'+obs.index.astype(str) if not obs.index.str.startswith(lib+'|').all() else obs.index
 assert set(obs.index)==set(m.index)
 m=m.loc[obs.index];libs=m.raw_library_size.to_numpy();ln=np.log1p(x*1e4/libs[:,None])
 for g,j in zip(genes,ix):availability.append(dict(library_id=lib,gene=g,present=j>=0))
 for tier in ['conservative','broad']:
  for gp in sorted(m['annotation_'+tier].unique()):
   mask=m['annotation_'+tier].eq(gp).to_numpy()
   for j,g in enumerate(genes):rows.append(dict(library_id=lib,sample_id=m.sample_id.iloc[0],patient_id=m.patient_id.iloc[0],patient_group=m.patient_group.iloc[0],tier=tier,population=gp,n_cells=int(mask.sum()),gene=g,detected_cells=int((x[mask,j]>0).sum()),fraction_detected=float((x[mask,j]>0).mean()),mean_logCP10k=float(ln[mask,j].mean())))
 my=m.annotation_conservative.isin(['Microglia','Macrophage_like','BAM_like','Myeloid_unresolved','Monocyte_like']).to_numpy()
 mm=m.loc[my,['library_id','sample_id','patient_id','patient_group','annotation_conservative']].copy()
 for j,g in enumerate(genes):mm[g]=ln[my,j]
 percell.append(mm);a.file.close();print('EXPRESSION COMPLETE',lib,flush=True)
expr=pd.DataFrame(rows);expr.to_csv(T/'S18_ALSF_expression_by_sample_population.csv',index=False)
pd.DataFrame(availability).to_csv(T/'S18_ALSF_gene_availability.csv',index=False)
my=pd.concat(percell);vals=my[genes];z=(vals-vals.mean())/vals.std(ddof=1);z=z.fillna(0)
for p,gs in modules.items():my['module_'+p]=z[gs].mean(1)
my['CD74_expression']=my['CD74'];my.index.name='cell_id';my.to_csv(T/'S18_ALSF_myeloid_expression_and_scores.csv.gz')
mc=[x for x in my if x.startswith('module_')];cols=mc+['CD74_expression'];srows=[]
for pop in ['All operational myeloid','Microglia','Macrophage_like','BAM_like','Myeloid_unresolved']:
 d=my if pop.startswith('All') else my[my.annotation_conservative.eq(pop)]
 sm=d.groupby(['library_id','sample_id','patient_id','patient_group'])[cols].mean().reset_index();sm['n_cells']=d.groupby(['library_id','sample_id','patient_id','patient_group']).size().to_numpy();sm['population']=pop;sm['eligible_min20']=sm.n_cells>=20;srows.append(sm)
samples=pd.concat(srows);samples.to_csv(T/'S18_ALSF_module_scores_by_sample.csv',index=False)
patients=samples[samples.eligible_min20].groupby(['population','patient_id'])[cols].mean().reset_index();patients.to_csv(T/'S18_ALSF_module_scores_by_patient.csv',index=False)
cr=[]
for pop,d in patients.groupby('population'):
 for p in mc:
  n=len(d);rho=float(spearmanr(d.CD74_expression,d[p]).statistic) if n>=3 else np.nan
  pv=np.nan
  if n>=3 and np.isfinite(rho):
   perm=[spearmanr(d.CD74_expression.to_numpy()[list(j)],d[p]).statistic for j in itertools.permutations(range(n))];pv=float(np.mean(np.abs(perm)>=abs(rho)-1e-12))
  cr.append(dict(population=pop,module=p,n_patients=n,rho=rho,p_exact_two_sided=pv,interpretation='Descriptive external cohort; sample means equally weighted within patient; fewer than 3 patients not tested'))
cr=pd.DataFrame(cr)
def bh(v):
 a=np.array(v);ok=np.isfinite(a);out=np.full(len(a),np.nan);ix=np.flatnonzero(ok)[np.argsort(a[ok])];out[ix]=np.minimum(1,np.minimum.accumulate((a[ix]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1]);return out
cr['q_BH_9_modules']=cr.groupby('population').p_exact_two_sided.transform(bh);cr.to_csv(T/'S18_ALSF_CD74_patient_associations.csv',index=False)
# Focal ligand / complete-receptor co-detection in the same library; raw detection fractions retained.
lr=[]
for tier in ['conservative','broad']:
 e=expr[expr.tier.eq(tier)]
 for lib,dd in e.groupby('library_id'):
  d=dd.set_index(['population','gene']);tum='Tumor_CNV_supported'
  if tum not in dd.population.unique():continue
  for pop in ['Microglia','Macrophage_like','BAM_like','Myeloid_unresolved']:
   if pop not in dd.population.unique():continue
   for lig,rec,src,tgt in [('APP','CD74',tum,pop),('SPP1','CD44',pop,tum),('PTN','PTPRZ1',pop,tum)]:
    l=d.loc[(src,lig)];r=d.loc[(tgt,rec)];lr.append(dict(library_id=lib,patient_id=l.patient_id,patient_group=l.patient_group,tier=tier,pair=lig+'_'+rec,source=src,target=tgt,n_source=int(l.n_cells),n_target=int(r.n_cells),source_ligand_fraction=l.fraction_detected,target_receptor_fraction=r.fraction_detected,source_ligand_mean=l.mean_logCP10k,target_receptor_mean=r.mean_logCP10k,eligible_min20=l.n_cells>=20 and r.n_cells>=20))
pd.DataFrame(lr).to_csv(T/'S18_ALSF_focal_LR_expression.csv',index=False)
# Composition and independent patient expression summary.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.3,1]})
groups=['Microglia','Macrophage_like','BAM_like','Myeloid_unresolved'];co=summary.pivot_table(index=['library_id','patient_id','patient_group'],columns='annotation_conservative',values='n_cells',fill_value=0).reindex(columns=groups,fill_value=0)
co.plot.bar(stacked=True,ax=ax[0],color=['#377f99','#ad6f39','#709358','#b7b7b7']);ax[0].set(title='A  Conservative myeloid counts by specimen',ylabel='Nuclei',xlabel='');ax[0].set_xticklabels([f'{l[-3:]} | {p}\n{g}' for l,p,g in co.index],rotation=65,ha='right',fontsize=8);ax[0].legend(fontsize=7,frameon=False)
pp=patients[patients.population.eq('All operational myeloid')];ap='module_antigen_presentation'
ax[1].scatter(pp.CD74_expression,pp[ap],s=55,c='#377f99')
for _,r in pp.iterrows():ax[1].annotate(r.patient_id,(r.CD74_expression,r[ap]),xytext=(5,5),textcoords='offset points',fontsize=8)
ax[1].set(title='B  Four independent patient summaries',xlabel='Mean myeloid CD74, ln(1 + CP10K)',ylabel='Antigen-presentation score excluding CD74')
fig.subplots_adjust(bottom=.35,wspace=.4)
for ext in ['png','pdf','svg']:fig.savefig(F/f'Figure_S8_external_cohort.{ext}',dpi=220,bbox_inches='tight')
summaryout={'n_nuclei':len(meta),'n_samples':meta.sample_id.nunique(),'n_patients':meta.patient_id.nunique(),'counts_conservative':meta.annotation_conservative.value_counts().to_dict(),'counts_broad':meta.annotation_broad.value_counts().to_dict(),'cd74_patient':cr.replace({np.nan:None}).to_dict('records'),'microglia_by_library':meta[meta.annotation_conservative.eq('Microglia')].library_id.value_counts().to_dict()}
(O/'data/alsf_expression_summary.json').write_text(json.dumps(summaryout,indent=2));print(json.dumps(summaryout,indent=2),flush=True)
