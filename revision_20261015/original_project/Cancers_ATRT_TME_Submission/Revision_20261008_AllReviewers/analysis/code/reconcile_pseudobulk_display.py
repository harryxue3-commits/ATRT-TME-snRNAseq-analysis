"""Reconcile displayed module differences to the original available-gene definition."""
from pathlib import Path
import pandas as pd
import numpy as np
import json
O=Path(__file__).resolve().parents[1];ROOT=O.parents[2];T=O/'tables'
v=pd.read_csv(T/'S12_paired_astrocyte_support_modules.csv')
original=pd.read_csv(ROOT/'PathwayAnalysis/results/support_mechanisms/geo_astrocyte_sensitivity/mechanism_module_scores_by_pseudobulk.csv',index_col=0)
v['original_available_gene_difference']=[original.loc[r.ID+'|hybrid_selected',r.module]-original.loc[r.ID+'|geo_tumor_like',r.module] for r in v.itertuples()]
v['absolute_reconstruction_error']=abs(v.present_genes_only_difference-v.original_available_gene_difference)
assert v.absolute_reconstruction_error.max()<1e-12
v['displayed_difference']=v.present_genes_only_difference
v['display_definition']='Mean over genes present in the matrix; unavailable genes omitted'
v['requested_panel_difference']=v.difference
v.to_csv(T/'S12_available_gene_display_and_reconciliation.csv',index=False)
(O/'data/S12_display_reconciliation.json').write_text(json.dumps({'all_24_original_donor_module_differences_reproduced':True,'max_absolute_error':float(v.absolute_reconstruction_error.max()),'unavailable_gene':'GPX1','affected_module':'Antioxidant_glutathione_support','available_genes':14,'requested_genes':15,'displayed_definition':'Available-gene mean, matching original Figure 3A','alternative_definition':'Requested-panel denominator with unavailable genes assigned zero, preserved in difference field of S12 source and S16 audit'},indent=2))
print('Verified all 24 original donor-module differences; reconciled the GPX1 denominator.')
