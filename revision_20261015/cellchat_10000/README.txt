ATRT CellChat: original 10,000-permutation analysis code
Prepared for collaborator sharing, 14 October 2026

WHAT IS INCLUDED
Unmodified archived R/C++ analysis scripts, result-summary Python scripts,
validation tables, and representative sessionInfo records. No expression data,
patient/cell metadata, or saved CellChat RDS inputs are included. The scripts
therefore require the original inputs before they can be rerun. This package
was assembled and integrity-checked; the expensive analyses were not rerun
while preparing it.

START HERE
Inside original_project/Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000/:
1. unified_permutations_20261012/code/rerun.R
   Donor-specific, conservative-annotation, MYC-2 inclusion, and external-cohort
   reruns. Explicit nboot=10000L; saved seed, cell IDs, donor assignments,
   annotations, expression, database, and selected interactions are preserved.
2. analysis/code/run_verified_10000.R
   Pooled MYC, SHH, and TYR models under state_resolved and all_author_tumor
   definitions. Explicit nboot=10000, seed.use=42.
3. controlled_unannotated_20261012/code/prepare_controlled.R and run_controlled.R
   Controlled comparison adding Unannotated malignant nuclei only, with fixed
   patients, group sizes, glial cells, and a common interaction family.
4. analysis/code/validated_fast_cellchat.R and fast_bootstrap_means.cpp
   Accelerated computation called by these drivers. It wraps CellChat's
   computeCommunProb while preserving its permutation generation and tests.
   The validation scripts and CSVs document comparisons with native outputs.

MAIN SETTINGS
CellChat 2.1.2; human database retained in each input object; RNA data;
type='truncatedMean'; trim=0.1; raw.use=TRUE; population.size=FALSE;
10,000 permutations. See saved sessionInfo files for the original environment.
The drivers require CellChat, Matrix, future, Rcpp, and jsonlite, plus a working
C++ compiler. Python summaries require their imported packages (including
NumPy and pandas). The archived Makevars and /private/tmp Rcpp cache location
are machine-specific; configure these for the collaborator's computer.

WHAT IS PERMUTED
CellChat generates random cell-group label assignments for its null scores.
Original expression, saved cell annotations, and donor/sample IDs are not
rewritten. Donor-specific models contain only their own donor; pooled subtype
models use their pooled cell inputs. This is not a donor-label permutation test.

P VALUES AND BH CORRECTION
The accelerated implementation preserves CellChat 2.1.2's strict null-score
exceedance rule. b is the count of null scores exceeding the observed score.
The final reporting convention is finite P=(b+1)/10001, with P=1 when the
observed score is zero. BH is applied across the full tested source-target-
interaction family within each model, before selecting highlighted candidates.
See unified_permutations_20261012/code/rerun.R and summarize.py for the final
reporting logic, and validate_saved.R for serialized-object/CSV checks.
Historical pooled scripts retain older native-nominal compatibility columns;
summarize.py harmonizes the final finite-P summaries. Do not use a historical
nominal_interactions.csv alone as the final revised reporting convention.
Pathway aggregation retains native CellChat P values; explicit pval_finite,
qval_BH, and final interaction tables are the finite-P/BH inference outputs.

RUNNING THE ORIGINAL DRIVERS
Restore the original relative project layout beneath PROJECT and provide the
saved inputs. INPUT_LOCATIONS.csv identifies the donor/external source folders;
each needs CellChat.rds and input_cells.csv. The original rerun_manifest.csv
contains absolute paths from the original computer: replace its source_dir
prefix with the collaborator's PROJECT location before running. The Rcpp cache
path may also need changing. Python summary scripts have original paths and
must be configured before use. Use a separate working copy; completion checks
intentionally prevent replacing completed runs.

Example commands (replace PROJECT with the actual project-root path):
Rscript PROJECT/Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000/unified_permutations_20261012/code/rerun.R PROJECT 0 1
The last two arguments mean worker 0 of 1 (all manifest rows, sequentially).

Rscript PROJECT/Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000/analysis/code/run_verified_10000.R PROJECT MYC state_resolved
Repeat for SHH/TYR and all_author_tumor as appropriate. These pooled drivers
require original 100-permutation CellChat objects under Revision_20261007/
analysis/cellchat and the supplied acceleration-equivalence validation table.

Rscript PROJECT/Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000/controlled_unannotated_20261012/code/prepare_controlled.R PROJECT
Rscript PROJECT/Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000/controlled_unannotated_20261012/code/run_controlled.R PROJECT MYC
Repeat the second command for SHH and TYR. Preparation requires the discovery
count matrix, genes, metadata, and existing pooled objects referenced in that
script. run_controlled.R needs prepared.rds in each model folder and completed
pooled state-resolved models for its comparison assertions.

VALIDATION AND PROVENANCE
validate_saved.R checks the 10,000-permutation setting, RDS/CSV agreement,
full-family BH calculation, and external group eligibility. Original validation
tables are provided for reference. Original script/table bytes are unchanged.
SHA256SUMS.txt covers all package files except itself. Historical figure/table
numbers in archived filenames precede manuscript renumbering; current main
Figure 6 corresponds to archived Figure 5, and external Figure S8 to old S9.
