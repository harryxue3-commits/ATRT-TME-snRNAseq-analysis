from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
OUT=Path(__file__).resolve().parents[1];CC=OUT/'cellchat';T=OUT/'tables'
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def edges(d):
 v=d[(d.p_value<.05)&(d.probability>0)&(d.source.ne(d.target))&((d.source=='Tumor')|(d.target=='Tumor'))&d.ligand.ne(d.receptor)]
 return set(zip(v.source,v.target,v.interaction_name))
def main():
 status=pd.read_csv(CC/'run_status.csv');finished=status[status.status.str.startswith('Complete')]
 data={};focus=[];thresholds=[]
 targets=[('Tumor','Microglia','APP_CD74'),('Microglia','Tumor','SPP1_ITGAV_ITGB1'),('Microglia','Tumor','SPP1_CD44'),('Astrocytes','Tumor','NCAM1_FGFR1'),('Astrocytes','Tumor','PTN_PTPRZ1')]
 for name in finished.run:
  f=CC/name/'all_tested_interactions.csv';d=pd.read_csv(f)
  # Native CellChat reports b/100; zeros denote finite permutation resolution, not P=0.
  d['p_value_plus1']=(np.rint(d.p_value*100)+1)/101
  d['q_BH_plus1_within_object']=bh(d.p_value_plus1)
  d=d.drop(columns=['q_BH_within_object'],errors='ignore')
  d.to_csv(f,index=False);d[d.nominal_p_lt_05].to_csv(CC/name/'nominal_interactions.csv',index=False)
  data[name]=d
  cells=pd.read_csv(CC/name/'input_cells.csv');sizes=cells.group.value_counts()
  for minimum in [5,10,20]:
   for group,n in sizes.items():thresholds.append(dict(run=name,group=group,n_cells=int(n),minimum=minimum,eligible=n>=minimum,meaning='Group-size eligibility only; not a refit of a reduced-group object'))
  for src,tgt,lr in targets:
   z=d[(d.source==src)&(d.target==tgt)&(d.interaction_name==lr)]
   if len(z):
    row=z.iloc[0];focus.append(dict(run=name,source=src,target=tgt,interaction_name=lr,status='Retained at native P<0.05' if row.p_value<.05 and row.probability>0 else 'Tested, not retained',probability=row.probability,p_native=row.p_value,p_plus1=row.p_value_plus1,q_BH_plus1=row.q_BH_plus1_within_object,source_cells=int(sizes.get(src,0)),target_cells=int(sizes.get(tgt,0))))
   else:focus.append(dict(run=name,source=src,target=tgt,interaction_name=lr,status='Group unavailable' if src not in sizes or tgt not in sizes else 'Not selected for inference',probability=np.nan,p_native=np.nan,p_plus1=np.nan,q_BH_plus1=np.nan,source_cells=int(sizes.get(src,0)),target_cells=int(sizes.get(tgt,0))))
 focused=pd.DataFrame(focus);focused.to_csv(T/'S13_focused_CellChat_candidates.csv',index=False)
 pd.DataFrame(thresholds).to_csv(T/'S13_group_size_eligibility.csv',index=False)
 compare=[]
 for sub in ['ATRT-MYC','ATRT-SHH','ATRT-TYR']:
  a=f'subtype__{sub}__state_resolved__exclude_MYC2';b=f'subtype__{sub}__all_author_tumor__exclude_MYC2'
  ea,eb=edges(data[a]),edges(data[b]);compare.append(dict(comparison='Tumor compartment expansion',subtype=sub,baseline=a,sensitivity=b,baseline_edges=len(ea),sensitivity_edges=len(eb),shared_edges=len(ea&eb),union_edges=len(ea|eb),jaccard=len(ea&eb)/len(ea|eb)))
 for scope in ['state_resolved','all_author_tumor']:
  a=f'subtype__ATRT-MYC__{scope}__exclude_MYC2';b=f'subtype__ATRT-MYC__{scope}__include_MYC2'
  ea,eb=edges(data[a]),edges(data[b]);compare.append(dict(comparison='Restore MYC2',subtype='ATRT-MYC',baseline=a,sensitivity=b,baseline_edges=len(ea),sensitivity_edges=len(eb),shared_edges=len(ea&eb),union_edges=len(ea|eb),jaccard=len(ea&eb)/len(ea|eb)))
 comp=pd.DataFrame(compare);comp.to_csv(T/'S13_network_sensitivity_comparisons.csv',index=False)
 contributions=pd.read_csv(CC/'donor_contributions.csv');contributions.to_csv(T/'S13_donor_contributions.csv',index=False)
 status.to_csv(T/'S13_run_status.csv',index=False)
 def count(lr,kind,scope='state_resolved'):
  q=focused[(focused.interaction_name==lr)&focused.run.str.startswith(kind+'__')]
  if kind=='subtype':q=q[q.run.str.contains('__'+scope+'__exclude_MYC2')]
  eligible=q[~q.status.eq('Group unavailable')]
  return int(eligible.status.eq('Retained at native P<0.05').sum()),len(eligible)
 app=count('APP_CD74','subtype');appwide=count('APP_CD74','subtype','all_author_tumor')
 nc=count('NCAM1_FGFR1','subtype');ptn=count('PTN_PTPRZ1','subtype')
 j=comp[comp.comparison.eq('Tumor compartment expansion')]
 myc=focused[(focused.interaction_name=='APP_CD74')&focused.run.isin(['subtype__ATRT-MYC__state_resolved__exclude_MYC2','subtype__ATRT-MYC__state_resolved__include_MYC2'])].set_index('run')
 p0=float(myc.loc['subtype__ATRT-MYC__state_resolved__exclude_MYC2','probability']);p1=float(myc.loc['subtype__ATRT-MYC__state_resolved__include_MYC2','probability'])
 donorapp=count('APP_CD74','donor');donornc=count('NCAM1_FGFR1','donor');donorptn=count('PTN_PTPRZ1','donor')
 cons=focused[focused.run.str.startswith('conservative__')]
 consastro=cons[(cons.source=='Astrocytes')&cons.interaction_name.eq('NCAM1_FGFR1')&cons.source_cells.ge(5)]
 consmicro=cons[(cons.target=='Microglia')&cons.interaction_name.eq('APP_CD74')&cons.target_cells.ge(5)]
 conservative_text=f"The frozen conservative annotations retain 23 source-labeled astrocytes and 15 source-labeled microglial nuclei. At a five-nucleus glial threshold, only {len(consastro)} tumor(s) supported a conservative astrocyte network and {len(consmicro)} supported a conservative microglial network; these cannot assess molecular-subtype reproducibility. Conservative microglia were not eligible in any donor at a ten-nucleus threshold. Macrophage/BAM-like populations were analyzed as a separate expression-defined group. The source malignant assignments were retained throughout."
 methods=("Fresh sensitivity fits used CellChat 2.1.2 (R 4.4.3) and its bundled full human database, with input expression normalized as ln(1 + CP10K), truncatedMean averaging (trim = 0.1), nboot = 100, seed = 42, and population.size = FALSE. These are separately versioned sensitivity fits, not exact reproductions of unavailable original objects. We compared the state-resolved tumor definition with all author-classified tumor nuclei, retaining the same group-specific eligibility thresholds and downsampling procedure, and repeated the MYC models with ATRT-MYC-2 restored to the eligible sample pool. Its tumor and astrocyte groups still failed their original thresholds, whereas its 18 microglial nuclei were eligible; restoring it also changed the per-donor myeloid downsampling target. Additional donor-specific fits used at least 20 state-resolved tumor nuclei, at least 10 source-labeled microglial nuclei, and at least 5 astrocytes. Exact input cell IDs, donor contributions, and database receptor-complex identities were saved. Native CellChat P < 0.05 is retained for comparison; 100 permutations give limited resolution, and native zero values are not exact zero probabilities. We additionally report (b + 1)/101 permutation P values and BH adjustment within each object. Neither measure tests replicated subtype specificity (Supplementary Table S13).")
 response_unannotated=(f"In these fresh fits, tumor-to-microglia APP–CD74 remained nominally retained in {app[0]}/{app[1]} state-resolved and {appwide[0]}/{appwide[1]} expanded subtype objects, but the full selected cross-compartment repertoires changed (Jaccard overlap {j.jaccard.min():.2f}–{j.jaccard.max():.2f}). This argues against assuming that exclusion has no effect. The expansion includes both Unannotated and Immune-like tumor cells and can change donor eligibility; it is a compartment-definition sensitivity, not an isolated causal effect of the Unannotated label.")
 response_myc=(f"Fresh MYC CellChat fits also restored MYC-2 to the sample pool. APP–CD74 remained nominally retained, with modeled probability {p0:.4f} without and {p1:.4f} with MYC-2. Other candidate repertoires changed, as reported in Table S13. These results support reporting sensitivity rather than claiming that this sample is inconsequential.")
 results=(f"Fresh, separately versioned CellChat fits retained tumor-to-microglia APP–CD74 in {app[0]}/{app[1]} state-resolved subtype objects and {appwide[0]}/{appwide[1]} expanded-tumor objects at the native nominal threshold. However, expanding the tumor compartment changed the selected cross-compartment repertoire (Jaccard overlap {j.jaccard.min():.2f}–{j.jaccard.max():.2f}), showing that the initial interactome is scope-dependent. Restoring ATRT-MYC-2 retained APP–CD74 in MYC, while its modeled probability changed from {p0:.4f} to {p1:.4f} and other candidates changed. Astrocyte-to-tumor NCAM1–FGFR1 was retained in {nc[0]}/{nc[1]} fresh state-resolved subtype objects; PTN–PTPRZ1 was retained in {ptn[0]}/{ptn[1]}, illustrating candidate-specific sensitivity. Donor-specific fits retained APP–CD74 in {donorapp[0]}/{donorapp[1]} eligible tumors, NCAM1–FGFR1 in {donornc[0]}/{donornc[1]} astrocyte-eligible tumors, and PTN–PTPRZ1 in {donorptn[0]}/{donorptn[1]}. These within-dataset results remain expression-based and do not establish functional signaling, independent validation, or subtype specificity (Supplementary Table S13).")
 primary_focus=focused[focused.run.str.contains('__state_resolved__exclude_MYC2')]
 assert not primary_focus.q_BH_plus1.lt(.05).any()
 fdr_note=' None of the five highlighted ligand–receptor candidates passed within-object BH correction in the fresh state-resolved subtype fits when finite-permutation P values were used. At 100 permutations, the minimum corrected permutation P is 1/101; nominal retention should therefore not be interpreted as multiplicity-controlled evidence.'
 results+=fdr_note
 response_unannotated+=' In the expanded TYR object, APP–CD74 was not selected for inference; this is not evidence of biological absence.'+fdr_note
 conservative_text=conservative_text.replace('1 tumor(s)','one tumor')
 summary_rows=[dict(comparison='Expanded tumor compartment',finding=f"APP–CD74 retained in {appwide[0]}/{appwide[1]} subtype objects; repertoire Jaccard overlap {j.jaccard.min():.2f}–{j.jaccard.max():.2f}."),dict(comparison='Restore ATRT-MYC-2',finding=f"MYC APP–CD74 probability {p0:.4f} to {p1:.4f}; nominal retention unchanged, other candidates changed."),dict(comparison='Same-donor models',finding=f"APP–CD74 {donorapp[0]}/{donorapp[1]}; NCAM1–FGFR1 {donornc[0]}/{donornc[1]}; PTN–PTPRZ1 {donorptn[0]}/{donorptn[1]} eligible tumors."),dict(comparison='Conservative populations',finding=f"Only {len(consastro)} eligible conservative astrocyte donor and {len(consmicro)} conservative microglial donor at the five-nucleus threshold; insufficient for subtype assessment."),dict(comparison='Multiple-testing sensitivity',finding='None of the five highlighted pairs passes within-object BH correction in the fresh state-resolved subtype fits using finite-permutation P values. Results are reported at the native nominal threshold for comparability.')]
 summary=dict(methods_text=methods,results_text=results,response_unannotated=response_unannotated,response_myc2=response_myc,conservative_text=conservative_text,summary_rows=summary_rows,n_completed=len(finished),status_rows=len(status),network_comparisons=compare)
 (T/'cellchat_revision_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False));print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
