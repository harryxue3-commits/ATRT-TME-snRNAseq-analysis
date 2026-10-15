from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

def false_discovery_control(values):
    # Independent NumPy implementation of BH, checked against R p.adjust below.
    values=np.asarray(values,dtype=float)
    order=np.argsort(values,kind='stable')
    scaled=values[order]*len(values)/np.arange(1,len(values)+1)
    adjusted=np.minimum(1,np.minimum.accumulate(scaled[::-1])[::-1])
    result=np.empty_like(adjusted);result[order]=adjusted
    return result

W=Path(__file__).resolve().parents[1]
REV=W.parent
PAIRS=[('APP_CD74','APP–CD74','Tumor','Microglia'),('NCAM1_FGFR1','NCAM1–FGFR1','Astrocytes','Tumor'),('SPP1_ITGAV_ITGB1','SPP1–(ITGAV + ITGB1)','Microglia','Tumor'),('SPP1_CD44','SPP1–CD44','Microglia','Tumor'),('PTN_PTPRZ1','PTN–PTPRZ1','Astrocytes','Tumor')]
checks=[]; focal=[]; runs=[]; comparisons=[]; sampling=[]; models={}; hashes=[]
def check(name,value):
    checks.append({'check':name,'passed':bool(value)})
    assert value,name
def edges(d,kind):
    mask=d.source.ne(d.target)&(d.source.eq('Tumor')|d.target.eq('Tumor'))&d.ligand.ne(d.receptor)&d.probability.gt(0)
    mask &= d.p_finite.lt(.05) if kind=='finite_nominal' else d.q_BH.lt(.05)
    return set(zip(d.loc[mask,'source'],d.loc[mask,'target'],d.loc[mask,'interaction_name']))
family=None
for sub in ['MYC','SHH','TYR']:
    inputs={}
    for scope in ['restricted','plus_unannotated']:
        name=f'{sub}__{scope}'; p=W/'models'/name
        check(name+' complete',(p/'COMPLETE.json').exists())
        d=pd.read_csv(p/'all_tested_interactions.csv'); m=pd.read_csv(p/'input_cells.csv')
        models[name]=d;inputs[scope]=m
        check(name+' unique cells',m.Cell.is_unique)
        check(name+' no Immune-like',(m.Final_Annotation!='Immune-like').all())
        check(name+' 10000 permutations',d.n_permutations.eq(10000).all())
        check(name+' integer exceedances',d.permutation_exceedances.between(0,10000).all())
        check(name+' finite P reproduced',np.allclose(d.p_finite,(d.permutation_exceedances+1)/10001,atol=1e-12,rtol=0))
        check(name+' zero observed P=1',d.loc[d.probability.eq(0),'p_finite'].eq(1).all())
        check(name+' BH independently reproduced',np.allclose(d.q_BH,false_discovery_control(d.p_finite.to_numpy()),atol=1e-12,rtol=0))
        check(name+' BH support reproduced',np.array_equal(d.BH_supported,d.q_BH.lt(.05)&d.probability.gt(0)))
        keys=set(zip(d.source,d.target,d.interaction_name))
        if family is None:family=keys
        check(name+' same tested family across all six',keys==family and len(keys)==len(d))
        sizes=m.group.value_counts()
        for pid,label,src,tgt in PAIRS:
            z=d[(d.interaction_name==pid)&(d.source==src)&(d.target==tgt)]
            check(name+' '+pid+' tested',len(z)==1)
            r=z.iloc[0]
            focal.append({'subtype':sub,'scope':scope,'interaction_name':pid,'candidate':label,'source':src,'target':tgt,'n_source':int(sizes[src]),'n_target':int(sizes[tgt]),'probability':r.probability,'exceedances':int(r.permutation_exceedances),'p_native':r.p_native,'p_finite':r.p_finite,'q_BH':r.q_BH,'BH_supported':bool(r.BH_supported),'n_tests':len(d),'n_permutations':10000})
        runs.append({'subtype':sub,'scope':scope,'n_cells':len(m),'n_tumor':int(sizes.Tumor),'n_unannotated':int(m.Final_Annotation.eq('Unannotated').sum()),'n_LR':d.interaction_name.nunique(),'n_tests':len(d),'n_permutations':10000,'nominal_finite_entries':int(d.nominal_finite.sum()),'BH_entries':int(d.BH_supported.sum())})
        tumor=m[m.group.eq('Tumor')]
        for donor,v in tumor.groupby('ID',sort=True):
            sampling.append({'subtype':sub,'scope':scope,'patient':donor,'n_tumor':len(v),'n_unannotated':int(v.Final_Annotation.eq('Unannotated').sum()),'n_immune_like':int(v.Final_Annotation.eq('Immune-like').sum())})
        if scope=='restricted':
            oldpath=REV/'analysis/cellchat'/f'subtype__ATRT-{sub}__state_resolved__exclude_MYC2'
            oldm=pd.read_csv(oldpath/'input_cells.csv')
            check(name+' exact original baseline cells',oldm.equals(m))
            prior=pd.read_csv(oldpath/'all_tested_interactions.csv').set_index(['source','target','interaction_name']).sort_index()
            here=d.set_index(['source','target','interaction_name']).loc[prior.index]
            check(name+' original observed scores unchanged',np.array_equal(here.probability.to_numpy(),prior.probability.to_numpy()))
            check(name+' original permutation P unchanged',np.array_equal(here.p_native.to_numpy(),prior.p_value.to_numpy()))
    a,b=inputs['restricted'],inputs['plus_unannotated']
    check(sub+' same donor/group order',a[['ID','group']].equals(b[['ID','group']]))
    check(sub+' exact glial metadata and cells',a[a.group.ne('Tumor')].equals(b[b.group.ne('Tumor')]))
    check(sub+' same per-donor group counts',a.groupby(['ID','group']).size().equals(b.groupby(['ID','group']).size()))
    check(sub+' only permitted tumor annotations',set(b.loc[b.group.eq('Tumor'),'Final_Annotation'])<=set(a.loc[a.group.eq('Tumor'),'Final_Annotation'])|{'Unannotated'})
    for kind in ['finite_nominal','BH_supported']:
        ea=edges(models[f'{sub}__restricted'],kind);eb=edges(models[f'{sub}__plus_unannotated'],kind)
        comparisons.append({'subtype':sub,'selection':kind,'restricted_edges':len(ea),'plus_unannotated_edges':len(eb),'shared':len(ea&eb),'union':len(ea|eb),'jaccard':len(ea&eb)/len(ea|eb)})
validation=pd.read_csv(W/'tables/S13_controlled_acceleration_validation.csv')
check('native acceleration equivalence',validation.all_probabilities_identical.all() and validation.all_pvalues_identical.all())
audit=pd.read_csv(W/'tables/S13_controlled_input_audit.csv')
check('preparation audits pass',audit[['patient_group_order_fixed','glial_cells_fixed','source_annotations_preserved']].to_numpy().all())
screen=pd.read_csv(W/'tables/S13_APP_CD74_gene_screen_audit.csv')
for scope,expected in [('historical_state_resolved',0.0100974596),('historical_all_author_tumor',0.0833429714)]:
    p=screen[(screen.subtype=='TYR')&(screen.scope==scope)&(screen.gene=='APP')&(screen.group=='Tumor')].screening_P.iloc[0]
    check('historical TYR APP screening '+scope,abs(p-expected)<1e-9)
pd.DataFrame(focal).to_csv(W/'tables/S13_controlled_focal_candidates.csv',index=False)
pd.DataFrame(runs).to_csv(W/'tables/S13_controlled_run_summary.csv',index=False)
pd.DataFrame(comparisons).to_csv(W/'tables/S13_controlled_interaction_overlap.csv',index=False)
pd.DataFrame(sampling).to_csv(W/'tables/S13_controlled_tumor_sampling.csv',index=False)
for p in sorted((W/'models').rglob('*')):
    if p.is_file():
        with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        hashes.append({'path':str(p.relative_to(W)),'bytes':p.stat().st_size,'sha256':digest})
pd.DataFrame(hashes).to_csv(W/'tables/S13_controlled_model_checksums.csv',index=False)
summary={'complete':True,'n_checks':len(checks),'all_checks_passed':all(x['passed'] for x in checks),'runs':runs,'focal_candidates':focal,'overlap':comparisons,'tumor_sampling':sampling,'checks':checks}
(W/'data/results.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
print(pd.DataFrame(focal)[['subtype','scope','candidate','probability','p_finite','q_BH','BH_supported']].to_string(index=False))
print(pd.DataFrame(comparisons).to_string(index=False))
print(f'{len(checks)} independent checks passed')
