"""Auditable all-reviewer analysis. Original malignant labels never changed."""
from pathlib import Path
import json, shutil, hashlib, itertools
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
O=Path(__file__).resolve().parents[1]; ROOT=O.parents[2]
G=ROOT/'Cancers_ATRT_TME_Submission/Revision_20261005/analysis'
R=ROOT/'Cancers_ATRT_TME_Submission/Revision_20261007/analysis'
A=ROOT/'ALSF-ATRT-snRNAseq/analysis_20261006'
T=O/'tables';F=O/'figures';T.mkdir(exist_ok=True);F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
C=json.loads((G/'code/analysis_config.json').read_text())
def save(d,n):d.to_csv(T/n,index=False)
def figsave(fig,n):
 for ext in ['png','pdf','svg']:fig.savefig(F/f'{n}.{ext}',dpi=220,bbox_inches='tight')
 plt.close(fig)
def dotplot(ax,df,groups,genes,title):
 for i,gp in enumerate(groups):
  d=df[df.group.eq(gp)].set_index('gene').reindex(genes)
  sc=ax.scatter(np.arange(len(genes)),np.full(len(genes),i),s=100*d.fraction, c=d.mean_logCP10k,cmap='viridis',vmin=0,vmax=3,edgecolors='none')
 ax.set(xticks=np.arange(len(genes)),xticklabels=genes,yticks=np.arange(len(groups)),yticklabels=[g.replace("Border-associated macrophage-like", "BAM-like") for g in groups],title=title)
 ax.tick_params(axis='x',labelrotation=65);ax.invert_yaxis();ax.set_xlim(-.7,len(genes)-.3)
 return sc

def main():
 # Preserve existing output rows exactly, not recomputed differently.
 copies={'astrocyte_program_sensitivity_summary.csv':'S16_astrocyte_annotation_programs.csv','astrocyte_program_donor_summary.csv':'S16_astrocyte_programs_by_donor.csv','support_module_sensitivity_summary.csv':'S16_astrocyte_support_summary.csv','support_module_paired_differences.csv':'S16_astrocyte_support_by_donor.csv','pseudobulk_eligibility.csv':'S16_pseudobulk_eligibility.csv','myeloid_subtype_sensitivity_tests.csv':'S16_myeloid_scoring_tests.csv','myeloid_program_donor_scores.csv':'S16_myeloid_scores_by_donor.csv','largest_z_label_sensitivity.csv':'S16_largest_relative_score_LOO.csv','ligand_receptor_donor_audit.csv':'S17_ligand_receptor_expression_by_donor.csv','ligand_receptor_sensitivity_summary.csv':'S17_ligand_receptor_expression_summary.csv'}
 for src,dst in copies.items():shutil.copy2(G/'tables'/src,T/dst)
 for p in (R/'tables').glob('*.csv'):shutil.copy2(p,T/p.name)
 a=ad.read_h5ad(G/'data/TME_annotated.h5ad');sets=json.loads((G/'data/analysis_sets.json').read_text())
 meta=pd.read_csv(G/'data/full_original_metadata.csv').set_index('Cell')
 info=a.obs.copy();info.index.name='cell_id';info.reset_index().to_csv(T/'S15_GSE_TME_cell_annotations.csv.gz',index=False)
 rows=[]
 for name,ids in sets.items():
  for donor in sorted(meta.ID.unique()):
   d=info.loc[info.index.isin(ids)&info.ID.eq(donor)]
   rows.append(dict(analysis_set=name,ID=donor,subtype=meta[meta.ID.eq(donor)].subtype.iloc[0],n_cells=len(d)))
 save(pd.DataFrame(rows),'S15_GSE_annotation_counts_by_donor.csv')
 markers=[]
 rationale={'Astrocyte':'Glutamate transport and differentiated astroglial identity; require coherent panel plus anchor and lineage exclusions.','Astrocyte_no_metabolic':'Sensitivity panel omits SLC1A2/SLC1A3 to assess circular selection of glutamate-associated expression.','Microglia':'Homeostatic microglial identity; activation and nuclear dropout reduce sensitivity.','Macrophage':'Macrophage-associated expression; overlapping activation markers do not prove blood-derived origin.','BAM':'Border-associated macrophage-like expression; anatomical identity is not established without spatial data.','Monocyte':'Inflammatory monocyte-associated expression; not interchangeable with macrophage identity.','Myeloid':'Broad myeloid backbone; these genes do not distinguish microglia from macrophages.'}
 for p,gs in C['panels'].items():
  for g in gs:markers.append(dict(panel=p,gene=g,rationale=rationale.get(p,'Competing lineage or ancillary TME panel; evaluated jointly with focal lineage markers.'),source_references='[11–13,24,27–28,30,39]; operational thresholds in Methods 2.11'))
 save(pd.DataFrame(markers),'S15_marker_panel_and_rationale.csv')
 # Show marker sharing with malignant glial-like states while never using it to relabel them.
 full=ad.read_h5ad(R/'data/GSE283839_all_counts.h5ad')
 assert full.obs_names.is_unique and full.n_obs==17564 and full.n_vars==36601
 markers1=['ALDH1L1','AQP4','SLC1A2','SLC1A3','GFAP','SOX9','PDGFRA','OLIG2','SOX10']
 markers2=['P2RY12','TMEM119','GPR34','SALL1','LST1','TYROBP','CD163','MRC1','MS4A7','LYVE1','F13A1','SPP1','APOE','APOC1','GPNMB']
 gs=list(dict.fromkeys(markers1+markers2+['APP','CD74','PTN','PTPRZ1','NCAM1','FGFR1','CD44','ITGAV','ITGB1']))
 ix=full.var_names.get_indexer(gs);assert (ix>=0).all()
 x=full.X[:,ix].toarray();libs=np.asarray(full.X.sum(1)).ravel();ln=np.log1p(x*1e4/libs[:,None]);obs=full.obs
 group_masks={'Source astrocytes (47)':obs.Final_Annotation.eq('Astrocytes').to_numpy(),'Conservative astrocytes (23)':obs.index.isin(sets['Conservative astrocytes']), 'TME OPC (79)':obs.Final_Annotation.eq('OPC').to_numpy(), 'Malignant OPC-like':obs.Final_Annotation.eq('OPC-like').to_numpy(), 'Malignant RG-like':obs.Final_Annotation.eq('RG-like').to_numpy(), 'Other author tumor':obs.Final_Annotation_Focus.ne('TME').to_numpy()&~obs.Final_Annotation.isin(['OPC-like','RG-like']).to_numpy(), 'Source microglia (288)':obs.Final_Annotation.eq('Microglia').to_numpy()}
 for lab in ['Microglia','Macrophage-like','Border-associated macrophage-like','Myeloid unresolved']:
  exact=[z for z in info.tme_annotation_revised.unique() if str(z).lower()==lab.lower()]
  if not exact and lab=='BAM-like':exact=[z for z in info.tme_annotation_revised.unique() if 'BAM' in str(z)]
  ids=info.index[info.tme_annotation_revised.isin(exact)];group_masks[f'{lab} ({len(ids)})']=obs.index.isin(ids)
 mr=[]
 for gp,m in group_masks.items():
  for j,g in enumerate(gs):mr.append(dict(group=gp,gene=g,n_cells=int(m.sum()),detected_cells=int((x[m,j]>0).sum()),fraction=float((x[m,j]>0).mean()) if m.any() else np.nan,mean_logCP10k=float(ln[m,j].mean()) if m.any() else np.nan))
 mr=pd.DataFrame(mr);save(mr,'S15_GSE_marker_expression.csv')
 fig,axs=plt.subplots(2,1,figsize=(11,9),gridspec_kw={'height_ratios':[1,1.1]})
 sc=dotplot(axs[0],mr,list(group_masks)[:6],markers1,'A  Astroglial markers and unchanged malignant labels')
 sc=dotplot(axs[1],mr,list(group_masks)[6:],markers2,'B  Source myeloid label and conservative TME populations')
 fig.subplots_adjust(hspace=.85,right=.88);cb=fig.colorbar(sc,ax=axs,fraction=.025,pad=.03,extend="max");cb.set_label('Mean ln(1 + CP10K)')
 handles=[axs[0].scatter([],[],s=100*f,color='#777777',label=f'{int(f*100)}%') for f in [.1,.5,1]]
 fig.legend(handles=handles,title='Detected nuclei',loc='center',bbox_to_anchor=(.5,.48),ncol=3,frameon=False)
 fig.text(.02,.01,'Dot area: transcript detection fraction. Marker panels informed annotation; these plots are an audit, not independent accuracy validation.',fontsize=8)
 figsave(fig,'Figure_S5_marker_audit')
 # Explicitly increasing eligibility thresholds on fixed previously fitted networks.
 # This filters unsupported groups; it does not refit group-dependent overexpression selection.
 rr=[];ef=[]
 pairs=[('APP_CD74','Tumor','Microglia'),('SPP1_CD44','Microglia','Tumor'),('SPP1_ITGAV_ITGB1','Microglia','Tumor'),('NCAM1_FGFR1','Astrocytes','Tumor'),('PTN_PTPRZ1','Astrocytes','Tumor')]
 for d in sorted((R/'cellchat').glob('donor__*'))+sorted((R/'cellchat').glob('conservative__*')):
  cells=pd.read_csv(d/'input_cells.csv');gc=cells.group.value_counts() if 'group' in cells else cells.groupby('labels').size()
  tested=pd.read_csv(d/'all_tested_interactions.csv') if (d/'all_tested_interactions.csv').exists() else pd.DataFrame()
  for minimum in [5,10,20]:
   for gp,n in gc.items():ef.append(dict(run=d.name,minimum=minimum,group=gp,n_cells=int(n),eligible=n>=max(20,minimum) if gp=='Tumor' else n>=minimum))
   for pair,src,tgt in pairs:
    elig=gc.get(src,0)>=max(20,minimum) if src=='Tumor' else gc.get(src,0)>=minimum
    elig &= gc.get(tgt,0)>=max(20,minimum) if tgt=='Tumor' else gc.get(tgt,0)>=minimum
    q=tested[tested.interaction_name.eq(pair)&tested.source.eq(src)&tested.target.eq(tgt)] if len(tested) else pd.DataFrame()
    status='Ineligible group' if not elig else ('No completed model' if not len(tested) else ('Not selected by model' if not len(q) else ('Native nominal' if bool(q.iloc[0].nominal_p_lt_05) else 'Tested not nominal')))
    rr.append(dict(run=d.name,minimum=minimum,interaction_name=pair,source=src,target=tgt,n_source=int(gc.get(src,0)),n_target=int(gc.get(tgt,0)),eligible=bool(elig),status=status,probability=float(q.iloc[0].probability) if len(q) else np.nan,q_BH=float(q.iloc[0].q_BH_plus1_within_object) if len(q) else np.nan,method='Higher-threshold post-fit eligibility filter; no refitting or BH family reduction'))
 save(pd.DataFrame(rr),'S17_GSE_CellChat_higher_threshold_candidates.csv');save(pd.DataFrame(ef),'S17_GSE_CellChat_higher_threshold_eligibility.csv')
 # Source group eligibility for all 12 donors, including those with no model.
 co=pd.read_csv(T/'S9_sample_counts_and_eligibility.csv');out=[]
 for _,r in co.iterrows():
  for k in [5,10,20]:out.append(dict(ID=r.ID,subtype=r.subtype,minimum=k,n_astro=r.Astrocytes,astro_eligible=r.Astrocytes>=k,n_micro=r.Microglia,micro_eligible=r.Microglia>=max(10,k),tumor_eligible=r.state_resolved_tumor>=20))
 save(pd.DataFrame(out),'S17_GSE_all_donor_threshold_eligibility.csv')
 # Replacement main Figure 4: continuous scores, no zero-cutoff polarization categories.
 scores=pd.read_csv(R/'tables/S10_microglia_scores.csv');programs=list(scores.columns[3:])
 print('scorecolumns',scores.columns.tolist(),flush=True)
 programs=[z for z in scores if z not in ['ID','subtype','Cell']]
 display={p:p.replace('_',' ').replace('M1 ','').replace('M2 ','') for p in programs}
 display.update({'inflammatory_M1':'Inflammatory associated','immunoregulatory_M2':'Immunoregulatory associated','M1_inflammatory':'Inflammatory associated','M2_immunoregulatory':'Immunoregulatory associated'})
 inf=next(p for p in programs if 'inflammatory' in p.lower());imm=next(p for p in programs if 'immunoregulatory' in p.lower())
 display[inf]='Inflammatory associated';display[imm]='Immunoregulatory associated'
 dom=scores[programs].idxmax(axis=1).value_counts();don=scores.groupby('ID')[programs].mean();dz=(don-don.mean())/don.std(ddof=1)
 display.update({inf:'Inflammatory\nassociated',imm:'Immunoregulatory\nassociated'})
 fig=plt.figure(figsize=(8.8,7.8));grid=fig.add_gridspec(2,2,wspace=1.0,hspace=.52);aa=fig.add_subplot(grid[0,0]);bb=fig.add_subplot(grid[0,1]);cc=fig.add_subplot(grid[1,0]);dd=fig.add_subplot(grid[1,1])
 colors={'ATRT-MYC':'#307e9e','ATRT-SHH':'#b0603c','ATRT-TYR':'#6f649b'}
 for sub,d in scores.groupby('subtype'):aa.scatter(d[inf],d[imm],s=13,c=colors[sub],alpha=.65,label=sub)
 aa.set(xlabel='Inflammatory-associated score',ylabel='Immunoregulatory-associated score',title='A  Continuous program scores');aa.legend(fontsize=8,frameon=False)
 aa.text(.02,.97,f'Spearman rho = {spearmanr(scores[inf],scores[imm]).statistic:.3f}\n288 nuclei; descriptive',transform=aa.transAxes,va='top',fontsize=8)
 bb.barh(np.arange(len(dom)),dom.values,color='#4d8494');bb.set(yticks=np.arange(len(dom)),yticklabels=[display[p] for p in dom.index],xlabel='Nuclei',title='B  Highest relative score');bb.invert_yaxis()
 for j,n in enumerate(dom):bb.text(n+.3,j,str(n),va='center',fontsize=8)
 im=cc.imshow(dz.T,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);cc.set(xticks=np.arange(len(don)),xticklabels=[z.replace('ATRT-','') for z in don.index],yticks=np.arange(len(programs)),yticklabels=[display[p] for p in programs],title='C  Standardized donor means');cc.tick_params(axis='x',labelrotation=65);fig.colorbar(im,ax=cc,shrink=.6,label='Relative donor score')
 mgs=['CSF1R','SLC2A5','CD86','STAT1','MSR1','CD163','CTSB','ABCA1','B2M','CD74','HIF1A','C1QB','TOP2A','CENPF']
 xx=full.X[obs.Final_Annotation.eq('Microglia').to_numpy()][:,full.var_names.get_indexer(mgs)].toarray();ll=libs[obs.Final_Annotation.eq('Microglia').to_numpy()]
 frac=(xx>0).mean(0)*100;mean=np.log1p(xx*1e4/ll[:,None]).mean(0)
 sc=dd.scatter(frac,np.arange(len(mgs)),c=mean,cmap='viridis',s=45);dd.set(yticks=np.arange(len(mgs)),yticklabels=mgs,xlabel='Nuclei with detected transcript (%)',title='D  Marker expression',xlim=(0,105));dd.invert_yaxis();fig.colorbar(sc,ax=dd,shrink=.6,label='Mean ln(1 + CP10K)')
 for j,n in enumerate((xx>0).sum(0)):dd.text(frac[j]+2,j,str(n),va='center',fontsize=8)
 figsave(fig,'Figure_4_revised')
 # Donor-specific ligand and complete-receptor detection, showing attrition not significance stars.
 lr=pd.read_csv(T/'S17_ligand_receptor_expression_by_donor.csv');parts=[]
 for aset,lig,rec,direction in [('Published microglia','APP','CD74','Tumor_to_TME'),('Published astrocytes','NCAM1','FGFR1','TME_to_tumor'),('Published astrocytes','PTN','PTPRZ1','TME_to_tumor')]:
  d=lr[(lr.analysis_set==aset)&(lr.ligand==lig)&(lr.receptor==rec)&(lr.direction==direction)&(lr.tumor_definition=='Resolved published states')].copy();d['pair']=lig+'–'+rec;parts.append(d)
 l=pd.concat(parts);save(l,'S17_highlighted_expression_plot_data.csv')
 fig,axs=plt.subplots(1,3,figsize=(8.5,5.8),sharey=True)
 for ax,(pair,d) in zip(axs,l.groupby('pair',sort=False)):
  d=d.sort_values('ID');yy=np.arange(len(d));ax.scatter(d.ligand_positive_fraction*100,yy,label='Source ligand',c='#227f9d',s=35);ax.scatter(d.receptor_complete_positive_fraction*100,yy,label='Target receptor',marker='s',c='#b26240',s=30)
  ax.set(yticks=yy,yticklabels=[r.ID.replace('ATRT-','') for _,r in d.iterrows()],xlim=(-3,105),xlabel='Detected nuclei (%)',title=pair);ax.invert_yaxis()
 axs[0].set_ylabel('Donor');fig.subplots_adjust(wspace=.20,bottom=.20);fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False); fig.suptitle('Same-donor expression; small groups retained descriptively',y=1.04);figsave(fig,'Figure_S7_LR_expression')
 # Annotation sensitivity visualization.
 counts=pd.read_csv(T/'S15_GSE_annotation_counts_by_donor.csv');ss=pd.read_csv(T/'S16_astrocyte_support_by_donor.csv');print('setnames',list(sets),flush=True)
 fig,axs=plt.subplots(1,2,figsize=(11,4.8),gridspec_kw={'width_ratios':[1,1.4]})
 sn=[s for s in sets if s in ['Published astrocytes','Conservative astrocytes','Expanded TME astrocytes','Nonmetabolic-identity astrocytes']]
 if len(sn)<4:sn=[s for s in sets if 'astro' in s.lower()][:4]
 vals=[len(sets[s]) for s in sn];axs[0].barh(range(len(sn)),vals,color='#568c9b');axs[0].set(yticks=range(len(sn)),yticklabels=[x.replace('Nonmetabolic identity','Non-transporter identity') for x in sn],xlabel='Nuclei',title='A  Annotation-dependent astrocyte recovery');axs[0].invert_yaxis()
 for j,n in enumerate(vals):axs[0].text(n+.5,j,str(n),va='center')
 d=ss[ss.module.str.contains('Glutamate',case=False)&ss.eligible];print('supportcols',d.columns.tolist(),flush=True)
 for aset,g in d.groupby('analysis_set'):
  j=list(d.analysis_set.unique()).index(aset);axs[1].scatter(g.difference,np.linspace(j-.09,j+.09,len(g)),s=45,label=aset)
  for _,r in g.iterrows():
   if aset=='Published astrocytes':axs[1].annotate(r.ID.replace('ATRT-',''),(r.difference,j),xytext=(2,10 if 'MYC' in r.ID else -15),textcoords='offset points',fontsize=7)
 axs[1].axvline(0,c='gray',lw=.8);axs[1].set(yticks=range(d.analysis_set.nunique()),yticklabels=d.analysis_set.unique(),xlabel='Astrocyte minus tumor, mean log2(CPM + 1)',title='B  Eligible donor glutamate-module differences');axs[1].invert_yaxis();fig.subplots_adjust(wspace=1.2);figsave(fig,'Figure_S6_annotation_sensitivity')
 print('GSE analysis finished',flush=True)

if __name__=='__main__':main()
