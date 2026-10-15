#!/usr/bin/env python3
"""Donor-aware sensitivity analyses with immutable published tumor labels."""
from pathlib import Path
import json, importlib.util
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, kruskal
import anndata as ad

OUT=Path(__file__).resolve().parents[1];ROOT=OUT.parents[2]
C=json.loads((OUT/'code/analysis_config.json').read_text())

def bh(p):
 p=np.asarray(p);o=np.argsort(p);q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];r=np.empty(len(p));r[o]=np.minimum(q,1);return r

def ucell(r,ix):
 n=len(ix);lo=n*(n+1)/2;return np.clip(1-(r[:,ix].sum(1)-lo)/(1500*n-lo),0,1)

def main():
 a=ad.read_h5ad(OUT/'data/TME_annotated.h5ad');e=a.obs.copy();x=a.layers['counts'].toarray();genes=list(a.var_names);gi={g:i for i,g in enumerate(genes)}
 lib=x.sum(1);log=np.log1p(x*(1e4/lib)[:,None]);r=np.minimum(rankdata(-x,axis=1,method='average'),1500)
 astro=e.Final_Annotation.eq('Astrocytes').to_numpy();mg=e.Final_Annotation.eq('Microglia').to_numpy()
 sets={
  'Published astrocytes':astro,
  'Conservative astrocytes':e.included_in_conservative_astrocytes.to_numpy(),
  'Expanded TME astrocytes':e.tme_annotation_revised.eq('Astrocyte').to_numpy(),
  'Nonmetabolic identity astrocytes':(astro&e.astrocyte_nonmetabolic_supported.to_numpy()),
  'Astrocytes identity k1':(astro&e.annotation_k1.eq('Astrocyte').to_numpy()),
  'Astrocytes identity k3':(astro&e.annotation_k3.eq('Astrocyte').to_numpy()),
  'Published microglia':mg,
  'Conservative microglia':e.included_in_conservative_microglia.to_numpy(),
  'Macrophage-like including BAM':e.tme_annotation_revised.isin(['Macrophage-like','Border-associated macrophage-like']).to_numpy(),
  'Myeloid unresolved':e.tme_annotation_revised.eq('Myeloid unresolved').to_numpy(),
 }
 (OUT/'data/analysis_sets.json').write_text(json.dumps({k:e.index[v].tolist() for k,v in sets.items()},indent=2))
 groups=[]
 for k,m in sets.items():
  for donor in sorted(e.ID.unique()):
   ii=m&e.ID.eq(donor).to_numpy();groups.append(dict(analysis_set=k,ID=donor,subtype=e.loc[e.ID.eq(donor),'subtype'].iloc[0],n_cells=int(ii.sum())))
 pd.DataFrame(groups).to_csv(OUT/'tables/analysis_set_donor_counts.csv',index=False)
 # Recompute UCell from the deposited raw counts and verify the saved baseline.
 old=pd.read_csv(ROOT/'ATRT_TME_Manuscript_Synthesis_2026-07-15/artifacts/astrocyte_state_ucell/cell_program_scores.csv')
 astro_scores={p:ucell(r,np.array([gi[g] for g in gs])) for p,gs in C['astro_programs'].items()}
 errors=[]
 for row in old.itertuples():errors.append(abs(astro_scores[row.program][e.index.get_loc(row.cell_id)]-row.ucell))
 assert max(errors)<1e-6,('UCell baseline mismatch',max(errors))
 # The original 10,000-null calibration remains fixed when restricting cells.
 ap=[]
 for name,mask in sets.items():
  if 'astrocyte' not in name.lower() or name.startswith('Expanded'):continue
  for p,gs in C['astro_programs'].items():
   orig=old[old.program.eq(p)].set_index('cell_id')
   for idx in np.flatnonzero(mask):
    cell=e.index[idx];row=orig.loc[cell];det=int((x[idx,[gi[g] for g in gs]]>0).sum())
    for k in [1,2,3]:ap.append({'analysis_set':name,'Cell':cell,'ID':e.iloc[idx].ID,'subtype':e.iloc[idx].subtype,'program':p,'state_minimum_detected_genes':k,'ucell':astro_scores[p][idx],'null_p95':row.null_p95,'detected_genes':det,'positive':bool(astro_scores[p][idx]-row.null_p95>1e-12 and det>=k)})
 # The 1e-12 comparison tolerance treats serialized equality at the cutoff as
 # equality, avoiding false positives caused solely by CSV roundoff.
 aps=pd.DataFrame(ap);aps.to_csv(OUT/'tables/astrocyte_program_sensitivity_per_cell.csv',index=False)
 check=aps[(aps.analysis_set=='Published astrocytes')&(aps.state_minimum_detected_genes==2)].set_index(['Cell','program']).positive
 original=old.set_index(['cell_id','program']).positive
 assert (check.loc[original.index].to_numpy()==original.to_numpy()).all(),'Original astrocyte positivity not reproduced'
 aps.groupby(['analysis_set','program','state_minimum_detected_genes']).agg(n_cells=('Cell','size'),positive_cells=('positive','sum'),positive_fraction=('positive','mean'),mean_ucell=('ucell','mean')).reset_index().to_csv(OUT/'tables/astrocyte_program_sensitivity_summary.csv',index=False)
 aps.groupby(['analysis_set','ID','subtype','program','state_minimum_detected_genes']).agg(n_cells=('Cell','size'),positive_cells=('positive','sum'),positive_fraction=('positive','mean'),mean_ucell=('ucell','mean')).reset_index().to_csv(OUT/'tables/astrocyte_program_donor_summary.csv',index=False)
 # Matched pseudobulk: sum raw counts, use all genes for library denominators.
 bulk=np.load(OUT/'data/fixed_tumor_summaries.npz');bgenes=list(bulk['genes']);assert bgenes==genes
 bnames=list(bulk['groupnames']);bcounts=bulk['counts'];source_meta=pd.read_csv(OUT/'data/full_original_metadata.csv')
 support=[];eligible=[];support_genes=[]
 for name,mask in sets.items():
  if 'astrocyte' not in name.lower():continue
  for donor in sorted(e.ID.unique()):
   ai=np.flatnonzero(mask&e.ID.eq(donor).to_numpy())
   bi=[j for j,n in enumerate(bnames) if n.split('|')[0]==donor and n.split('|')[1] in C['tumor_states']]
   nt=int((source_meta.ID.eq(donor)&source_meta.Final_Annotation.isin(C['tumor_states'])).sum())
   ok=len(ai)>=5 and nt>=20
   eligible.append(dict(analysis_set=name,ID=donor,n_astrocytes=len(ai),n_resolved_tumor=nt,eligible=ok,reason='Eligible' if ok else 'Fewer than 5 astrocytes or 20 resolved tumor nuclei'))
   if not len(ai) or not nt:continue
   ac=x[ai].sum(0);tc=bcounts[bi].sum(0);av=np.log2(ac/max(ac.sum(),1)*1e6+1);tv=np.log2(tc/max(tc.sum(),1)*1e6+1)
   for module,gs in C['support_modules'].items():
    ix=[gi[g] for g in gs if g in gi];missing=[g for g in gs if g not in gi]
    # The supplied baseline script zero-fills absent requested genes. Preserve
    # that denominator for reproduction, and expose a present-gene sensitivity.
    am=av[ix].sum()/len(gs);tmv=tv[ix].sum()/len(gs);ds=am-tmv
    support.append(dict(analysis_set=name,ID=donor,subtype=e.loc[e.ID.eq(donor),'subtype'].iloc[0],n_astrocytes=len(ai),n_tumor=nt,eligible=ok,module=module,n_requested_genes=len(gs),n_present_genes=len(ix),missing_genes=';'.join(missing),astro_mean_log2CPM=float(am),tumor_mean_log2CPM=float(tmv),difference=float(ds),present_genes_only_difference=float((av[ix]-tv[ix]).mean())))
    for g in gs:support_genes.append(dict(analysis_set=name,ID=donor,module=module,gene=g,eligible=ok,gene_in_matrix=g in gi,astro_counts=int(ac[gi[g]]) if g in gi else np.nan,tumor_counts=int(tc[gi[g]]) if g in gi else np.nan,astro_library=int(ac.sum()),tumor_library=int(tc.sum()),astro_log2CPM=float(av[gi[g]]) if g in gi else np.nan,tumor_log2CPM=float(tv[gi[g]]) if g in gi else np.nan))
 elig=pd.DataFrame(eligible);elig.to_csv(OUT/'tables/pseudobulk_eligibility.csv',index=False)
 sup=pd.DataFrame(support);sup.to_csv(OUT/'tables/support_module_paired_differences.csv',index=False)
 pd.DataFrame(support_genes).to_csv(OUT/'tables/support_gene_pseudobulk.csv',index=False)
 sr=[]
 for (name,module),d in sup[sup.eligible].groupby(['analysis_set','module']):
  sr.append(dict(analysis_set=name,module=module,n_donors=len(d),median_difference=d.difference.median(),minimum_difference=d.difference.min(),maximum_difference=d.difference.max(),interpretation='Descriptive; insufficient donor replication for generalization'))
 pd.DataFrame(sr).to_csv(OUT/'tables/support_module_sensitivity_summary.csv',index=False)
 # Microglial scores use the published 288-cell scale for comparability.
 defs=pd.read_csv(ROOT/'Supplementary/Code/04_microglia_state_analysis/microglia_state_program_gene_sets.tsv',sep='\t')
 names={'m1_inflammatory':'inflammatory_associated','m2_immunoregulatory':'immunoregulatory_associated'}
 ms=[];usage=[];zmat={};all_z=np.zeros_like(log)
 means=log[mg].mean(0);std=log[mg].std(0,ddof=1);all_z=np.divide(log-means,std,out=np.zeros_like(log),where=std>0)
 for module,d in defs.groupby('module',sort=False):
  gs=list(d.gene);ret=[g for g in gs if g in gi and (x[mg,gi[g]]>0).sum()>=9];ix=[gi[g] for g in ret]
  z=all_z[:,ix].mean(1);rankscore=ucell(r,np.array(ix));uncut=[gi[g] for g in gs if g in gi];z_all=all_z[:,uncut].mean(1)
  zmat[names.get(module,module)]=z
  for g in gs:usage.append(dict(program=names.get(module,module),gene=g,detected_published_microglia=int((x[mg,gi[g]]>0).sum()) if g in gi else 0,retained_original_3percent=g in ret))
  for idx in np.flatnonzero(mg):ms.append(dict(Cell=e.index[idx],ID=e.iloc[idx].ID,subtype=e.iloc[idx].subtype,revised_annotation=e.iloc[idx].tme_annotation_revised,program=names.get(module,module),n_genes=len(ret),fixed_baseline_z=z[idx],ucell=rankscore[idx],z_without_detection_filter=z_all[idx]))
 mscore=pd.DataFrame(ms);mscore.to_csv(OUT/'tables/myeloid_program_scores.csv',index=False);pd.DataFrame(usage).to_csv(OUT/'tables/myeloid_program_gene_usage.csv',index=False)
 donor=[];test=[]
 for name,mask in sets.items():
  if 'astrocyte' in name.lower():continue
  selected=mscore[mscore.Cell.isin(e.index[mask])]
  for (sample,sub,program),d in selected.groupby(['ID','subtype','program']):
   donor.append(dict(analysis_set=name,ID=sample,subtype=sub,program=program,n_cells=len(d),eligible_5_cells=len(d)>=5,mean_fixed_z=d.fixed_baseline_z.mean(),mean_ucell=d.ucell.mean(),mean_z_all_genes=d.z_without_detection_filter.mean()))
 ds=pd.DataFrame(donor);ds.to_csv(OUT/'tables/myeloid_program_donor_scores.csv',index=False)
 for (name,program),d in ds[ds.eligible_5_cells].groupby(['analysis_set','program']):
  subcounts=d.groupby('subtype').size();ok=len(subcounts)==3 and min(subcounts)>=2
  for method,col in [('fixed_z','mean_fixed_z'),('ucell','mean_ucell'),('no_gene_filter','mean_z_all_genes')]:
   p=kruskal(*[v[col].to_numpy() for _,v in d.groupby('subtype')]).pvalue if ok else np.nan
   test.append(dict(analysis_set=name,program=program,method=method,n_donors=len(d),MYC_donors=int(subcounts.get('ATRT-MYC',0)),SHH_donors=int(subcounts.get('ATRT-SHH',0)),TYR_donors=int(subcounts.get('ATRT-TYR',0)),p_value=p,status='Exploratory donor-level test' if ok else 'Not assessable: fewer than 2 eligible donors in a subtype'))
 tests=pd.DataFrame(test);tests['q_BH_within_set_method']=np.nan
 for _,ix in tests.groupby(['analysis_set','method']).groups.items():
  keep=[i for i in ix if pd.notna(tests.loc[i,'p_value'])]
  if keep:tests.loc[keep,'q_BH_within_set_method']=bh(tests.loc[keep,'p_value'].to_numpy())
 tests.to_csv(OUT/'tables/myeloid_subtype_sensitivity_tests.csv',index=False)
 influence=[]
 for (name,program),d in ds[ds.eligible_5_cells].groupby(['analysis_set','program']):
  baseline=d.mean_fixed_z.mean()
  for sample in d.ID:
   other=d[d.ID.ne(sample)]
   influence.append(dict(analysis_set=name,program=program,omitted_donor=sample,n_remaining=len(other),all_donor_mean=baseline,leave_one_out_mean=other.mean_fixed_z.mean(),change=other.mean_fixed_z.mean()-baseline))
 pd.DataFrame(influence).to_csv(OUT/'tables/myeloid_leave_one_donor_out.csv',index=False)
 # Sensitivity of the old largest-z summary to single-gene omission.
 zdf=pd.DataFrame(zmat,index=e.index);peak=zdf.idxmax(1);peakrows=[]
 for module,d in defs.groupby('module',sort=False):
  key=names.get(module,module);ret=[g for g in d.gene if (x[mg,gi[g]]>0).sum()>=9]
  for g in ret:
   zz=zdf.copy();zz[key]=all_z[:,[gi[h] for h in ret if h!=g]].mean(1)
   changed=zz.idxmax(1)!=peak
   peakrows.append(dict(program=key,omitted_gene=g,n_published_microglia=288,n_largest_z_label_changed=int(changed[mg].sum()),fraction_changed=float(changed[mg].mean()),interpretation='Largest relative score is not dominant biological function'))
 pd.DataFrame(peakrows).to_csv(OUT/'tables/largest_z_label_sensitivity.csv',index=False)
 # Cell-level and donor-level CD74 association, omitting CD74 from the antigen module.
 corr=[];cd=log[:,gi['CD74']]
 for module in ['antigen_presentation','m2_immunoregulatory']:
  gs=[g for g in defs.loc[defs.module.eq(module),'gene'] if g!='CD74' and (x[mg,gi[g]]>0).sum()>=9]
  z=all_z[:,[gi[g] for g in gs]].mean(1)
  for level in ['cell_descriptive','donor_means']:
   dat=pd.DataFrame({'ID':e.loc[mg,'ID'],'cd74':cd[mg],'program':z[mg]})
   if level=='donor_means':dat=dat.groupby('ID')[['cd74','program']].mean()
   rho,p=spearmanr(dat.cd74,dat.program);corr.append(dict(program=names.get(module,module),level=level,n=len(dat),rho=rho,p_value=p if level=='donor_means' else np.nan,note='CD74 excluded from program to avoid self-correlation; cell correlation descriptive only'))
 corr=pd.DataFrame(corr);corr['q_BH_donor_tests']=np.nan;di=corr.level.eq('donor_means');corr.loc[di,'q_BH_donor_tests']=bh(corr.loc[di,'p_value']);corr.to_csv(OUT/'tables/CD74_program_correlations.csv',index=False)
 # Highlighted marker expression by refined identity and donor.
 expr=[]
 highlight=['P2RY12','TMEM119','GPR34','SALL1','CD163','MRC1','MS4A7','F13A1','LYVE1','SPP1','APOE','APOC1','GPNMB','CD74','APP']
 for name,mask in sets.items():
  if 'astrocyte' in name.lower():continue
  for donorid in ['All']+sorted(e.ID.unique()):
   ii=mask if donorid=='All' else mask&e.ID.eq(donorid).to_numpy()
   for g in highlight:expr.append(dict(analysis_set=name,ID=donorid,gene=g,n_cells=int(ii.sum()),n_detected=int((x[ii,gi[g]]>0).sum()),fraction_detected=float((x[ii,gi[g]]>0).mean()) if ii.any() else np.nan,mean_log1pCP10k=float(log[ii,gi[g]].mean()) if ii.any() else np.nan))
 pd.DataFrame(expr).to_csv(OUT/'tables/highlighted_myeloid_marker_expression.csv',index=False)
 # Complex-aware within-donor transcript co-detection; not CellChat probability.
 full=source_meta.set_index('Cell').loc[bulk['barcodes']];sx=bulk['selected_counts'];sg={g:i for i,g in enumerate(bulk['selected_genes'])};sl=np.log1p(sx*(1e4/bulk['library_sizes'])[:,None]);pos={g:i for i,g in enumerate(bulk['barcodes'])}
 lrsets={k:v for k,v in sets.items() if not k.startswith('Astrocytes identity')}
 lrrows=[]
 for name,mask in lrsets.items():
  for sample in sorted(e.ID.unique()):
   aa=np.array([pos[cell] for cell in e.index[mask&e.ID.eq(sample).to_numpy()]],dtype=int)
   for tumor_def in ['Resolved published states','All author non-TME']:
    tm=full.ID.eq(sample)&(full.Final_Annotation.isin(C['tumor_states']) if tumor_def.startswith('Resolved') else full.Final_Annotation_Focus.ne('TME'))
    tt=np.flatnonzero(tm)
    for direction,si,ti in [('TME_to_tumor',aa,tt),('Tumor_to_TME',tt,aa)]:
     for pair in C['lr_pairs']:
      ligand=pair['ligand'];receptor=pair['receptor'];ris=[sg[g] for g in receptor.split('+')]
      ld=(sx[si,sg[ligand]]>0);rd=np.all(sx[ti][:,ris]>0,axis=1)
      lv=sl[si,sg[ligand]];rv=np.prod(sl[ti][:,ris],axis=1)**(1/len(ris))
      row=dict(analysis_set=name,ID=sample,subtype=e.loc[e.ID.eq(sample),'subtype'].iloc[0],tumor_definition=tumor_def,direction=direction,ligand=ligand,receptor=receptor,receptor_note=pair['note'],n_source=len(si),n_target=len(ti),ligand_positive_cells=int(ld.sum()),receptor_complete_positive_cells=int(rd.sum()),ligand_positive_fraction=float(ld.mean()) if len(si) else np.nan,receptor_complete_positive_fraction=float(rd.mean()) if len(ti) else np.nan,ligand_mean_log1pCP10k=float(lv.mean()) if len(si) else np.nan,receptor_geometric_mean_log1pCP10k=float(rv.mean()) if len(ti) else np.nan,expression_product=float(lv.mean()*rv.mean()) if len(si) and len(ti) else np.nan)
      for threshold in [5,10,20]:row[f'eligible_n{threshold}']=len(si)>=threshold and len(ti)>=threshold
      row['codetected_same_donor']=bool(ld.any() and rd.any());row['both_fractions_ge_10percent']=bool(len(si) and len(ti) and ld.mean()>=.1 and rd.mean()>=.1)
      lrrows.append(row)
 lr=pd.DataFrame(lrrows);lr.to_csv(OUT/'tables/ligand_receptor_donor_audit.csv',index=False)
 lrs=[]
 for keys,d in lr.groupby(['analysis_set','tumor_definition','direction','ligand','receptor']):
  for threshold in [5,10,20]:
   v=d[d[f'eligible_n{threshold}']];lrs.append(dict(zip(['analysis_set','tumor_definition','direction','ligand','receptor'],keys),minimum_cells=threshold,n_eligible_donors=len(v),n_codetected_donors=int(v.codetected_same_donor.sum()),n_donors_both_fractions_ge_10percent=int(v.both_fractions_ge_10percent.sum()),median_expression_product=v.expression_product.median(),eligible_donors=';'.join(v.ID),codetected_donors=';'.join(v.loc[v.codetected_same_donor,'ID'])))
 pd.DataFrame(lrs).to_csv(OUT/'tables/ligand_receptor_sensitivity_summary.csv',index=False)
 src=source_meta[source_meta.Final_Annotation_Focus.ne('TME')].groupby(['ID','subtype','Final_Annotation']).size().rename('n_cells').reset_index();src.to_csv(OUT/'tables/fixed_tumor_state_donor_contributions.csv',index=False)
 # Label-threshold stability is a sensitivity analysis, not validation accuracy.
 sensitiv=[]
 for col in ['annotation_k1','annotation_k2','annotation_k3','annotation_strict_backbone_null']:
  for label,n in e[col].value_counts().items():sensitiv.append(dict(rule=col,label=label,n_cells=int(n)))
 pd.DataFrame(sensitiv).to_csv(OUT/'tables/annotation_threshold_sensitivity.csv',index=False)
 summary={'ucell_baseline_max_absolute_error':max(errors),'analysis_set_sizes':{k:int(v.sum()) for k,v in sets.items()},'pseudobulk_eligible_samples':{k:v.loc[v.eligible,'ID'].tolist() for k,v in elig.groupby('analysis_set')},'LR_method':'Within-donor ligand and all receptor-subunit co-detection plus descriptive expression product. Not CellChat; no probability or permutation p-value generated.','CellChat_rerun_status':'Original model objects, analysis script, database version and full settings absent; plotted values were not used as source data.'}
 (OUT/'tables/sensitivity_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
