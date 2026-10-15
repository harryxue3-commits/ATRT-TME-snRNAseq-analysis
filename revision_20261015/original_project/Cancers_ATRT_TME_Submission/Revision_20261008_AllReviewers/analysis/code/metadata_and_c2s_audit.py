from pathlib import Path
import pandas as pd,json,hashlib,shutil
O=Path(__file__).resolve().parents[1];ROOT=O.parents[2];T=O/'tables';D=O/'data'
f=next((D/'atlas_supplement').rglob('*.xlsx'))
s1=pd.read_excel(f,'Supp. Table S1');s2=pd.read_excel(f,'Supp. Table S2')
m=pd.read_csv(ROOT/'Cancers_ATRT_TME_Submission/Revision_20261005/analysis/data/full_original_metadata.csv');cells=m.Cell
q=s2[s2.Cell.isin(cells)].set_index('Cell').loc[cells]
assert len(q)==17564 and q.index.is_unique
assert (q.Final_Annotation.to_numpy()==m.Final_Annotation.to_numpy()).all()
assert (q.Final_Annotation_Focus.to_numpy()==m.Final_Annotation_Focus.to_numpy()).all()
assert (q.ID.to_numpy()==m.ID.to_numpy()).all()
rows=[]
for id,d in q.groupby('ID'):
 s=s1[s1.id.eq(id)].iloc[0];sex=list(d.gender.dropna().unique());orig=d['orig.ident'].unique();assert len(orig)==1
 rows.append(dict(ID=id,source_id=orig[0],molecular_subtype=m[m.ID.eq(id)].subtype.iloc[0],extended_subtype=s.ATRT_subtype,tissue_preparation=s.tissue_type,sequencing_as_reported=s['Sequencing technology'],age_years='Not reported',sex=';'.join(sex) if sex else 'Not reported',tumor_grade='Not reported',site='Not reported',treatment='Not reported',outcome='Not reported',n_nuclei=len(d),source='Blanco-Carmona et al. 2025 Supplementary Tables S1 and S2; exact cell-ID match to GSE283839 subset'))
pd.DataFrame(rows).to_csv(T/'S14_GSE_clinical_metadata.csv',index=False)
text='The analyzed cohort is the 12-tumor 10X v3 subset of the larger 19-tumor atlas. The published supplement reports snap-frozen tissue for all 12 and extended labels SHH-1A for SHH-1 through SHH-4 and SHH-2 for SHH-5; broad labels are MYC (3), SHH (5), and TYR (4). Age, tumor grade, site, treatment, and outcomes are not reported for these 12 tumors in the deposited GEO metadata or the relevant published patient table. The published per-nucleus gender field is missing for all 12; sex values for the seven other atlas tumors were not transferred to this subset. Table S14 lists source identifiers, reported tissue and technology, and explicit missing fields. The 17,564 cell identifiers and malignant/TME assignments match the published supplement exactly.'
(D/'clinical_results.txt').write_text(text)
(D/'clinical_verification.json').write_text(json.dumps({'article_doi':'10.1093/neuonc/noaf179','download_endpoint':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12916745/supplementaryFiles','inner_zip_published_md5':'29b4372a2baebfda9b30de1b70fb8106','inner_zip_observed_md5':hashlib.md5((D/'original_atlas_supplement.zip').read_bytes()).hexdigest(),'workbook_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'exact_cell_matches':len(q),'malignant_labels_preserved':True,'missing_sex_for_12':True},indent=2))
source=ROOT/'Supplementary/Code/05_cell2sentence_analysis';dest=O/'code/original_C2S_implementation';dest.mkdir(exist_ok=True)
for p in source.glob('*.py'):shutil.copy2(p,dest/p.name)
pred=ROOT/'PathwayAnalysis/PathwayAnalysisV2/GSE283839_ATRT_RNA_V3_non_malignant_TME_27B_cell_type_predictions_opt.csv';d=pd.read_csv(pred)
print('C2S prediction shape',d.shape,'columns',d.columns.tolist());assert len(d)==516
out=T/'S19_C2S_saved_predictions.csv.gz';d.to_csv(out,index=False)
settings=[('Model','vandijklab/C2S-Scale-Gemma-2-27B','Recorded in supplied run template'),('Organism','Homo sapiens','Recorded in template'),('Gene ordering','Up to 1,000 genes by decreasing expression; raw-count ranks','Saved cell_sentence and generation script'),('Input token limit','4,096','Supplied template; truncates lowest-ranked genes to fit'),('Output token limit','24','Supplied template'),('Generation','Greedy; do_sample=False','Supplied template'),('Seed','42','Supplied template and ranked-sentence construction'),('Batch size','8','Supplied template; not an execution log'),('Quantization','4-bit NF4; double quantization; bfloat16 compute','Supplied template'),('Padding','Left','Supplied template'),('Cell inclusion','All 516 source TME nuclei; no confidence cutoff','Verified saved output row count and original comparison'),('Label harmonization','Eight broad categories; mapping function archived','Original analyze_cell2sentence_annotation_concordance.py'),('Prompt','The following is a list of 1000 gene names ordered by descending expression level in a Homo sapiens cell. Your task is to give the cell type which this cell belongs to based on its gene expression.\nCell sentence: [ranked genes].\nThe cell type corresponding to these genes is:','Template text; full script supplied verbatim'),('Exact TME execution log','Not retained','Cannot verify every template setting was used in the original run'),('Immutable model revision','Not recorded','Cannot guarantee exact historical model weights'),('Installed package versions','Not pinned in supplied template','Cannot reconstruct exact historical execution environment'),('Input paths in template','Example tumor-cell paths remain in original template','Template provenance only; not asserted as exact TME invocation'),('Saved prediction SHA-256',hashlib.sha256(pred.read_bytes()).hexdigest(),'Original uncompressed CSV'),('Interpretation','Same-data method comparison; not independent validation','Methods 2.5')]
pd.DataFrame(settings,columns=['setting','value','evidence_or_limitation']).to_csv(T/'S19_C2S_implementation_audit.csv',index=False)
print(text)
