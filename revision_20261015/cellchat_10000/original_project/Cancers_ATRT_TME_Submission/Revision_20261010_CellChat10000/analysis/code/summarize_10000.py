from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

REV=Path(__file__).resolve().parents[2];A=REV/'analysis';T=A/'tables'
OLD=REV.parent/'Revision_20261007/analysis/cellchat'
PLAN=json.loads((A/'data/analysis_plan.json').read_text())
PAIRS=[('APP_CD74','APP–CD74','Tumor','Microglia'),
       ('NCAM1_FGFR1','NCAM1–FGFR1','Astrocytes','Tumor'),
       ('SPP1_ITGAV_ITGB1','SPP1–(ITGAV + ITGB1)','Microglia','Tumor'),
       ('SPP1_CD44','SPP1–CD44','Microglia','Tumor'),
       ('PTN_PTPRZ1','PTN–PTPRZ1','Astrocytes','Tumor')]
checks=[];focus=[];comparisons=[];data={};runs=[];hashes=[]
def check(name,result):
    checks.append({'check':name,'passed':bool(result)});assert result,name
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def edges(d,kind):
    v=d[d.source.ne(d.target)&((d.source=='Tumor')|(d.target=='Tumor'))&d.ligand.ne(d.receptor)&d.probability.gt(0)]
    v=v[v.p_value.lt(.05) if kind=='native' else v.q_BH_plus1_within_object.lt(.05)]
    return set(zip(v.source,v.target,v.interaction_name))
for model in PLAN['models']:
    name=model['run'];path=A/'cellchat'/name;sub=name.split('__')[1].replace('ATRT-','');scope=name.split('__')[2]
    check(name+' complete', (path/'COMPLETE.json').exists())
    d=pd.read_csv(path/'all_tested_interactions.csv');old=pd.read_csv(OLD/name/'all_tested_interactions.csv')
    cells=pd.read_csv(path/'input_cells.csv');prior=pd.read_csv(OLD/name/'input_cells.csv')
    check(name+' source RDS unchanged',sha(OLD/name/'CellChat.rds')==model['source_rds_sha256'])
    check(name+' exact input cells, labels, patient IDs and order',cells.equals(prior))
    check(name+' unique cell IDs',cells.Cell.is_unique)
    key=['source','target','interaction_name'];x=d.set_index(key).sort_index();y=old.set_index(key).sort_index()
    check(name+' tested family unchanged',x.index.equals(y.index))
    check(name+' observed scores unchanged',np.array_equal(x.probability.to_numpy(),y.probability.to_numpy()))
    check(name+' exactly 10000 permutations',d.n_permutations.eq(10000).all() and d.seed.eq(42).all())
    check(name+' finite permutation correction',np.allclose(d.p_value_plus1,(np.rint(d.p_value*10000)+1)/10001,rtol=0,atol=1e-12))
    check(name+' independent full-family BH reproduction',np.allclose(d.q_BH_plus1_within_object,false_discovery_control(d.p_value_plus1.to_numpy()),rtol=0,atol=1e-12))
    check(name+' nominal classification',np.array_equal(d.nominal_p_lt_05,(d.p_value<.05)&(d.probability>0)))
    check(name+' BH classification',np.array_equal(d.BH_q_lt_05,(d.q_BH_plus1_within_object<.05)&(d.probability>0)))
    data[name]=d;sizes=cells.group.value_counts()
    runs.append({'run':name,'subtype':sub,'scope':scope,'n_cells':len(cells),'n_tests':len(d),'n_permutations':10000,'seed':42,'nominal_entries':int(d.nominal_p_lt_05.sum()),'BH_entries':int(d.BH_q_lt_05.sum())})
    for pid,label,src,tgt in PAIRS:
        z=x.loc[(src,tgt,pid)] if (src,tgt,pid) in x.index else None
        v=y.loc[(src,tgt,pid)] if (src,tgt,pid) in y.index else None
        row={'run':name,'subtype':sub,'scope':scope,'interaction_name':pid,'label':label,'source':src,'target':tgt,'n_source':int(sizes.get(src,0)),'n_target':int(sizes.get(tgt,0)),'n_tests':len(d),'n_permutations':10000,'selected':z is not None}
        if z is not None:
            row.update(probability=z.probability,exceedances=int(z.permutation_exceedances),p_native=z.p_value,p_plus1=z.p_value_plus1,q_BH=z.q_BH_plus1_within_object,nominal=bool(z.nominal_p_lt_05),BH_supported=bool(z.BH_q_lt_05),p_plus1_100=v.p_value_plus1,q_BH_100=v.q_BH_plus1_within_object)
            row['status']='BH supported' if row['BH_supported'] else 'Native nominal only' if row['nominal'] else 'Tested not nominal'
        else:row.update(probability=np.nan,exceedances=np.nan,p_native=np.nan,p_plus1=np.nan,q_BH=np.nan,nominal=False,BH_supported=False,p_plus1_100=np.nan,q_BH_100=np.nan,status='Not selected by model')
        focus.append(row)
    for p in path.iterdir():
        if p.is_file():hashes.append({'path':str(p.relative_to(REV)),'bytes':p.stat().st_size,'sha256':sha(p)})
for sub in ['MYC','SHH','TYR']:
    a=f'subtype__ATRT-{sub}__state_resolved__exclude_MYC2';b=f'subtype__ATRT-{sub}__all_author_tumor__exclude_MYC2'
    for kind in ['native','BH']:
        ea,eb=edges(data[a],kind),edges(data[b],kind)
        comparisons.append({'subtype':sub,'selection':kind,'n_permutations':10000,'baseline_edges':len(ea),'expanded_edges':len(eb),'shared':len(ea&eb),'union':len(ea|eb),'jaccard':len(ea&eb)/len(ea|eb) if ea|eb else np.nan})
f=pd.DataFrame(focus);f.to_csv(T/'S13_10000_focal_candidates.csv',index=False)
pd.DataFrame(runs).to_csv(T/'S13_10000_run_summary.csv',index=False)
pd.DataFrame(comparisons).to_csv(T/'S13_10000_tumor_scope_comparisons.csv',index=False)
pd.DataFrame(hashes).to_csv(T/'S13_10000_model_checksums.csv',index=False)
summary={'complete':True,'n_models':len(runs),'n_permutations':10000,'n_checks':len(checks),'all_passed':all(x['passed'] for x in checks),'runs':runs,'focal_candidates':f.replace({np.nan:None}).to_dict('records'),'scope_comparisons':comparisons,'checks':checks}
(A/'data/results_10000.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False))
print(f[['scope','subtype','label','p_plus1','q_BH','status']].to_string(index=False))
print(pd.DataFrame(comparisons).to_string(index=False))
print(f'{len(checks)} validation checks passed.')
