"""Reviewer 1 sensitivity analyses; source malignant labels remain immutable."""
from pathlib import Path
import json, hashlib, shutil, platform
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import rankdata, spearmanr, kruskal
from scipy import sparse, io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parents[1]
ROOT=OUT.parents[2]
OLD=ROOT/'Cancers_ATRT_TME_Submission/Revision_20261005/analysis'
def bh(v):
    p=np.asarray(v,float); result=np.full(len(p),np.nan); ok=np.flatnonzero(np.isfinite(p))
    order=ok[np.argsort(p[ok])]
    result[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return result
def save(df,name): df.to_csv(OUT/'tables'/name,index=False)
def main():
    a=ad.read_h5ad(OLD/'data/TME_annotated.h5ad'); obs=a.obs.copy()
    genes=list(a.var_names); gi={g:i for i,g in enumerate(genes)}
    x=a.layers['counts'].toarray().astype(float)
    assert len(obs)==516 and (obs.Final_Annotation_Focus=='TME').all()
    meta=pd.read_csv(OLD/'data/full_original_metadata.csv')
    assert len(meta)==17564 and meta.Cell.is_unique
    lib=x.sum(1); log=np.log1p(x/lib[:,None]*1e4)
    am=obs.Final_Annotation.eq('Astrocytes').to_numpy(); mm=obs.Final_Annotation.eq('Microglia').to_numpy()
    config=json.loads((OLD/'code/analysis_config.json').read_text())
    original=pd.read_csv(ROOT/'ATRT_TME_Manuscript_Synthesis_2026-07-15/artifacts/astrocyte_state_ucell/cell_program_scores.csv')
    ranks=np.minimum(rankdata(-x[am],axis=1,method='average'),1500)
    cells=list(obs.index[am]); rows=[]; errors=[]
    for p,gs in config['astro_programs'].items():
        ix=[gi[g] for g in gs]; n=len(ix); low=n*(n+1)/2
        scores=np.clip(1-(ranks[:,ix].sum(1)-low)/(1500*n-low),0,1)
        old=original[original.program.eq(p)].set_index('cell_id').loc[cells]
        errors.extend(np.abs(scores-old.ucell.to_numpy()).tolist())
        detected=(x[am][:,ix]>0).sum(1)
        for j,c in enumerate(cells):
            for k in [1,2,3]:
                positive=bool(scores[j]-old.iloc[j].null_p95>1e-12 and detected[j]>=k)
                if k==2: assert positive==old.iloc[j].positive
                rows.append(dict(Cell=c,ID=obs.loc[c,'ID'],subtype=obs.loc[c,'subtype'],program=p,k=k,ucell=scores[j],null_p95=old.iloc[j].null_p95,detected_genes=int(detected[j]),positive=positive))
    assert max(errors)<1e-12
    ap=pd.DataFrame(rows); save(ap,'S8_astrocyte_threshold_per_cell.csv')
    aps=ap.groupby(['program','k']).agg(n_cells=('Cell','size'),n_positive=('positive','sum'),fraction=('positive','mean')).reset_index()
    save(aps,'S8_astrocyte_threshold_summary.csv')
    combos=ap.pivot(index=['Cell','ID','k'],columns='program',values='positive').reset_index()
    combos['positive_program_count']=combos[list(config['astro_programs'])].sum(1)
    save(combos,'S8_astrocyte_threshold_overlap.csv')
    # Cohort composition and published CellChat input eligibility, both sample scopes.
    rows=[]
    for (sample,sub),d in meta.groupby(['ID','subtype']):
        states=d.Final_Annotation.isin(config['tumor_states'])
        counts=d.Final_Annotation.value_counts(); tme=int(d.Final_Annotation_Focus.eq('TME').sum())
        rows.append(dict(ID=sample,subtype=sub,total=len(d),TME=tme,TME_percent=100*tme/len(d),Astrocytes=int(counts.get('Astrocytes',0)),Microglia=int(counts.get('Microglia',0)),state_resolved_tumor=int(states.sum()),all_author_tumor=int(d.Final_Annotation_Focus.ne('TME').sum()),Unannotated=int(counts.get('Unannotated',0))))
    cohort=pd.DataFrame(rows); save(cohort,'S9_sample_counts_and_eligibility.csv')
    crows=[]
    for scope,sel in [('All 12 samples',meta),('Without ATRT-MYC-2',meta[meta.ID.ne('ATRT-MYC-2')])]:
        for sub,d in sel.groupby('subtype'):
            t=d.Final_Annotation_Focus.eq('TME'); counts=d.Final_Annotation.value_counts()
            crows.append(dict(scope=scope,subtype=sub,n_samples=d.ID.nunique(),total=len(d),TME=int(t.sum()),TME_percent=100*t.mean(),astrocytes=int(counts.get('Astrocytes',0)),microglia=int(counts.get('Microglia',0))))
    save(pd.DataFrame(crows),'S9_composition_with_without_MYC2.csv')
    # Exact gene-level z-score reconstruction, with the original 288-cell scale fixed.
    defs=pd.read_csv(ROOT/'Supplementary/Code/04_microglia_state_analysis/microglia_state_program_gene_sets.tsv',sep='\t')
    m=obs.loc[mm,['ID','subtype']].copy(); m['Cell']=m.index
    mx=x[mm]; ml=log[mm]; sd=ml.std(0,ddof=1)
    z=np.divide(ml-ml.mean(0),sd,out=np.zeros_like(ml),where=sd>0)
    usage=[]; modules={}; modules_loo={}
    for name,d in defs.groupby('module',sort=False):
        retained=[]
        for g in d.gene:
            count=int((mx[:,gi[g]]>0).sum()) if g in gi else 0
            keep=count>=9
            usage.append(dict(module=name,gene=g,detected_cells=count,percent=100*count/288,retained=keep))
            if keep: retained.append(g)
        modules[name]=retained
        modules_loo[name]=[g for g in retained if g!='CD74']
        m[name]=z[:,[gi[g] for g in retained]].mean(1)
    save(pd.DataFrame(usage),'S10_microglia_gene_filter_audit.csv')
    save(m.reset_index(drop=True),'S10_microglia_scores.csv')
    donor=m.groupby(['ID','subtype'],observed=True).agg(n_cells=('Cell','size'),**{name:(name,'mean') for name in modules}).reset_index()
    assert len(donor)==12 and donor.n_cells.sum()==288
    save(donor,'S9_donor_program_scores.csv')
    tests=[]; baselines=[]
    for exclude in [False,True]:
        d=donor[donor.ID.ne('ATRT-MYC-2')] if exclude else donor
        for name in modules:
            test=kruskal(*[g[name] for _,g in d.groupby('subtype')])
            tests.append(dict(exclude_MYC2=exclude,module=name,n_donors=len(d),p_value=float(test.pvalue)))
            for sub,g in d.groupby('subtype'):
                baselines.append(dict(exclude_MYC2=exclude,module=name,subtype=sub,n_donors=len(g),mean_donor_score=g[name].mean(),median_donor_score=g[name].median()))
    tests=pd.DataFrame(tests); tests['q_BH']=tests.groupby('exclude_MYC2').p_value.transform(bh)
    save(tests,'S9_subtype_tests_with_without_MYC2.csv'); save(pd.DataFrame(baselines),'S9_subtype_baselines_with_without_MYC2.csv')
    # CD74 is excluded from every correlated module, preventing part-whole circularity.
    cd=pd.DataFrame({'Cell':m.Cell,'ID':m.ID,'subtype':m.subtype,'CD74_logCP10k':ml[:,gi['CD74']],'CD74_detected':mx[:,gi['CD74']]>0,'log_library_size':np.log1p(mx.sum(1))})
    for name,gs in modules_loo.items(): cd[name]=z[:,[gi[g] for g in gs]].mean(1)
    save(cd.reset_index(drop=True),'S11_CD74_cell_data.csv')
    dc=cd.groupby(['ID','subtype'],observed=True).agg(n_cells=('Cell','size'),CD74_logCP10k=('CD74_logCP10k','mean'),CD74_detected_fraction=('CD74_detected','mean'),log_library_size=('log_library_size','mean'),**{n:(n,'mean') for n in modules}).reset_index()
    assert len(dc)==12 and dc.n_cells.sum()==288 and not dc.isna().any().any()
    save(dc,'S11_CD74_donor_data.csv')
    cr=[]
    for scope,d,do_test in [('Donor means (all)',dc,True),('Donor means (>=5 nuclei)',dc[dc.n_cells>=5],True),('Donor means (without MYC2)',dc[dc.ID.ne('ATRT-MYC-2')],True),('Pooled nuclei (descriptive)',cd,False)]:
        if do_test:
            rng=np.random.default_rng(20261007)
            rx=rankdata(d.CD74_logCP10k);rx=(rx-rx.mean())/np.linalg.norm(rx-rx.mean())
            permuted=rx[rng.random((99999,len(d))).argsort(axis=1)]
        for name in modules:
            r,p=spearmanr(d.CD74_logCP10k,d[name])
            if do_test:
                ry=rankdata(d[name]);ry=(ry-ry.mean())/np.linalg.norm(ry-ry.mean())
                p=(1+np.count_nonzero(np.abs(permuted@ry)>=abs(r)-1e-12))/100000
            cr.append(dict(scope=scope,module=name,n=len(d),rho=float(r),p_value=float(p) if do_test else np.nan,p_method='Two-sided donor permutation; 99,999 permutations; seed 20261007' if do_test else 'Descriptive only'))
    cor=pd.DataFrame(cr);cor['q_BH_9_modules']=cor.groupby('scope').p_value.transform(bh)
    save(cor,'S11_CD74_correlations.csv')
    # Descriptive partial rank correlation removing donor intercepts and library depth.
    nuisance=np.column_stack([np.ones(len(cd)),pd.get_dummies(cd.ID,drop_first=True,dtype=float).to_numpy(),rankdata(cd.log_library_size)])
    rx=rankdata(cd.CD74_logCP10k); rx=rx-nuisance@np.linalg.lstsq(nuisance,rx,rcond=None)[0]
    residual=[]
    for name in modules:
        ry=rankdata(cd[name]);ry=ry-nuisance@np.linalg.lstsq(nuisance,ry,rcond=None)[0]
        residual.append(dict(module=name,n_cells=len(cd),rho_partial=float(np.corrcoef(rx,ry)[0,1]),adjustment='Donor indicators and rank(log library size); descriptive, no cell-level inferential P value'))
    save(pd.DataFrame(residual),'S11_CD74_depth_donor_adjusted.csv')
    within=[]
    for sample,d in cd.groupby('ID'):
        for name in modules:
            r=spearmanr(d.CD74_logCP10k,d[name]).statistic if len(d)>=5 and d.CD74_logCP10k.nunique()>1 else np.nan
            within.append(dict(ID=sample,module=name,n_cells=len(d),rho=r))
    save(pd.DataFrame(within),'S11_CD74_within_donor_correlations.csv')
    supp=pd.read_csv(OLD/'tables/support_module_paired_differences.csv')
    save(supp[(supp.analysis_set=='Published astrocytes')&supp.eligible],'S12_paired_astrocyte_support_modules.csv')
    # Exact original source-data provenance; new scores are not independent validation.
    manifest=json.loads((OLD/'data/input_manifest.json').read_text())
    manifest.update(analysis_date='2026-10-07',original_ucell_max_error=max(errors),python=platform.python_version(),scope='Reviewer 1; malignant annotations unchanged; original 10,000-null astrocyte calibration frozen')
    (OUT/'data/reviewer1_manifest.json').write_text(json.dumps(manifest,indent=2))
    # Export all source counts for a fresh, explicitly versioned CellChat sensitivity run.
    fullpath=OUT/'data/GSE283839_all_counts.h5ad'
    if not fullpath.exists():
        source=ROOT/'NCBIGEO/snRNA/GSE283839_ATRT_RNA_V3_counts_raw.tsv'
        import re
        vv=[]; ii=[]; jj=[]; gn=[]
        with source.open() as f:
            barcodes=[re.sub(r'^([ACGT]+)\.([0-9]+)',r'\1-\2',b) for b in f.readline().rstrip('\r\n').split('\t')]
            for j,line in enumerate(f):
                g,values=line.rstrip('\r\n').split('\t',1); v=np.fromstring(values,sep='\t',dtype=np.int32)
                assert len(v)==len(barcodes)
                nz=np.flatnonzero(v); vv.append(v[nz]); ii.append(nz.astype(np.int32)); jj.append(np.full(len(nz),j,dtype=np.int32));gn.append(g)
        full=sparse.csr_matrix((np.concatenate(vv),(np.concatenate(ii),np.concatenate(jj))),shape=(len(barcodes),len(gn)))
        mo=meta.set_index('Cell').loc[barcodes].copy()
        ad.AnnData(full,obs=mo,var=pd.DataFrame(index=gn)).write_h5ad(fullpath,compression='gzip')
        io.mmwrite(OUT/'data/GSE283839_counts.mtx',full.T)
        pd.Series(gn).to_csv(OUT/'data/GSE283839_genes.tsv',index=False,header=False)
        mo.reset_index().to_csv(OUT/'data/GSE283839_metadata.csv',index=False)
    fig=plt.figure(figsize=(7,6.5),layout='constrained')
    grid=fig.add_gridspec(2,2,height_ratios=[.85,1])
    axes=[fig.add_subplot(grid[0,:]),fig.add_subplot(grid[1,0]),fig.add_subplot(grid[1,1])]
    names=list(config['astro_programs']); short=['Homeostatic','Pan-reactive','Inflammatory','ECM-remodeling']
    for k,offset,color in [(1,-.24,'#a5bacc'),(2,0,'#346886'),(3,.24,'#152f46')]:
        vals=aps[aps.k==k].set_index('program').loc[names]
        axes[0].bar(np.arange(4)+offset,vals.fraction*100,width=.23,label=f'≥{k} genes',color=color)
    axes[0].set_xticks(np.arange(4),short,rotation=25,ha='right');axes[0].set_ylabel('Program-positive nuclei (%)');axes[0].legend(fontsize=8);axes[0].set_title('A  Astrocyte threshold sensitivity (n=47)',fontsize=10)
    colors={'ATRT-MYC':'#9567ae','ATRT-SHH':'#de8b48','ATRT-TYR':'#408b85'}
    for ax,name,title in zip(axes[1:],['antigen_presentation','m2_immunoregulatory'],['B  Antigen presentation\n(CD74 excluded)','C  Immunoregulatory-associated\nprogram']):
        for sub,d in dc.groupby('subtype'): ax.scatter(d.CD74_logCP10k,d[name],s=35,c=colors[sub],label=sub.replace('ATRT-',''))
        row=cor[(cor.scope=='Donor means (all)')&(cor.module==name)].iloc[0]
        ax.set_title(title,fontsize=10);ax.set_xlabel('Mean CD74 [ln(1 + CP10k)]',fontsize=9);ax.set_ylabel('Mean relative module score',fontsize=9);ax.text(.03,.97,f'n={row.n}; Spearman rho={row.rho:.3f}\nBH q={row.q_BH_9_modules:.3g}',transform=ax.transAxes,va='top',fontsize=8)
    axes[2].legend(fontsize=8,loc='lower right')
    for ax in axes: ax.spines[['top','right']].set_visible(False)
    fig.savefig(OUT/'figures/Supplementary_Figure_S3.png',dpi=250);fig.savefig(OUT/'figures/Supplementary_Figure_S3.pdf');plt.close(fig)
    summary={'astro_threshold_counts':aps.to_dict('records'),'excluded_microglia_markers':pd.DataFrame(usage).query('not retained').to_dict('records'),'CD74_correlations':cor[cor.scope=='Donor means (all)'].to_dict('records'),'MYC2_counts':cohort[cohort.ID.eq('ATRT-MYC-2')].to_dict('records'),'minimum_BH_with_MYC2':tests[~tests.exclude_MYC2].q_BH.min(),'minimum_BH_without_MYC2':tests[tests.exclude_MYC2].q_BH.min()}
    (OUT/'tables/reviewer1_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
