args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]); worker <- as.integer(args[2]); workers <- as.integer(args[3])
rev <- file.path(root,'Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000')
out <- file.path(rev,'unified_permutations_20261012')
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(future)})
stopifnot(as.character(packageVersion('CellChat'))=='2.1.2')
plan(sequential);options(future.globals.maxSize=12*1024^3)
source(file.path(rev,'analysis/code/validated_fast_cellchat.R'))
fn <- make_fast_cellchat(file.path(rev,'analysis/code'))
manifest <- read.csv(file.path(out,'tables/rerun_manifest.csv'),stringsAsFactors=FALSE)
apply_filter <- function(obj,row) {
 if(row$scope=='external') return(filterCommunication(obj,min.cells=row$minimum-1L))
 sizes<-table(obj@idents)
 for(g in names(sizes)) {
  cutoff<-if(g=='Astrocytes')5 else if(g=='Tumor')10 else if(row$minimum==5)5 else 10
  if(sizes[[g]]<cutoff) {obj@net$prob[g,,]<-0;obj@net$prob[,g,]<-0;obj@net$pval[g,,]<-1;obj@net$pval[,g,]<-1}
 }
 obj
}
for(i in seq_len(nrow(manifest))) {
 if((i-1L)%%workers!=worker)next
 row<-manifest[i,];name<-row$run;dest<-file.path(out,'cellchat',name)
 dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 if(file.exists(file.path(dest,'COMPLETE.json')))next
 message('START ',name,' ',format(Sys.time()))
 started<-Sys.time();old<-readRDS(file.path(row$source_dir,'CellChat.rds'))
 seed<-old@options$parameter$seed.use;oldB<-old@options$parameter$nboot
 stopifnot(oldB %in% c(100,1000),old@options$parameter$trim==.1,!old@options$parameter$population.size)
 input<-read.csv(file.path(row$source_dir,'input_cells.csv'),check.names=FALSE,stringsAsFactors=FALSE)
 cellkey<-if(row$scope=='external')'cell_id' else 'Cell'
 stopifnot(identical(colnames(old@data),input[[cellkey]]),identical(rownames(old@meta),input[[cellkey]]))
 # Validate acceleration against the full saved native result for each model.
 check<-apply_filter(fn(old,type='truncatedMean',trim=.1,raw.use=TRUE,population.size=FALSE,nboot=oldB,seed.use=seed),row)
 dp<-max(abs(old@net$prob-check@net$prob));dv<-max(abs(old@net$pval-check@net$pval))
 stopifnot(dp<1e-12,dv==0)
 write.csv(data.frame(run=name,original_permutations=oldB,seed=seed,probability_max_difference=dp,p_value_max_difference=dv,all_p_values_identical=TRUE),file.path(dest,'acceleration_validation.csv'),row.names=FALSE)
 rm(check);gc()
 obj<-apply_filter(fn(old,type='truncatedMean',trim=.1,raw.use=TRUE,population.size=FALSE,nboot=10000L,seed.use=seed),row)
 stopifnot(identical(old@meta,obj@meta),identical(old@idents,obj@idents),identical(old@data,obj@data),identical(old@data.signaling,obj@data.signaling),identical(old@DB,obj@DB),identical(old@LR,obj@LR),identical(dimnames(old@net$prob),dimnames(obj@net$prob)))
 delta<-max(abs(old@net$prob-obj@net$prob));stopifnot(delta<1e-12)
 pr<-as.data.frame.table(obj@net$prob,responseName='probability',stringsAsFactors=FALSE);names(pr)[1:3]<-c('source','target','interaction_name')
 pr$p_native<-as.vector(obj@net$pval);pr$permutation_exceedances<-as.integer(round(pr$p_native*10000))
 pr$p_finite<-(pr$permutation_exceedances+1)/10001;pr$p_finite[pr$probability==0]<-1
 pr$q_BH<-p.adjust(pr$p_finite,'BH')
 pr$native_nominal<-pr$p_native<.05 & pr$probability>0
 pr$finite_nominal<-pr$p_finite<.05 & pr$probability>0
 pr$BH_supported<-pr$q_BH<.05 & pr$probability>0
 pr$n_permutations<-10000L;pr$seed<-seed;pr$run<-name
 # Compatibility aliases retain native P while all primary nominal summaries use finite P.
 pr$p_value<-pr$p_native;pr$p_value_plus1<-pr$p_finite;pr$p_plus1<-pr$p_finite
 pr$q_BH_plus1_within_object<-pr$q_BH;pr$q_BH_within_object<-pr$q_BH
 pr$nominal_p_lt_05<-pr$finite_nominal;pr$BH_q_lt_05<-pr$BH_supported
 info<-obj@DB$interaction;info$interaction_name<-rownames(info);pr<-merge(pr,info,by='interaction_name',all.x=TRUE,sort=FALSE)
 write.csv(pr,file.path(dest,'all_tested_interactions.csv'),row.names=FALSE)
 write.csv(pr[pr$finite_nominal,],file.path(dest,'nominal_interactions.csv'),row.names=FALSE)
 write.csv(pr[pr$BH_supported,],file.path(dest,'BH_supported_interactions.csv'),row.names=FALSE)
 obj@options$null_calculation_implementation<-'Validated computational acceleration; native RNG and strict exceedance rule preserved'
 obj@options$reporting<-'Finite P=(b+1)/10001; zero observed score P=1; BH over full within-object family; nominal means finite P<0.05'
 # Keep native pval in its original slot for provenance; use explicit finite/BH arrays for inference.
 obj@net$pval_finite<-array((round(obj@net$pval*10000)+1)/10001,dim=dim(obj@net$pval),dimnames=dimnames(obj@net$pval))
 obj@net$pval_finite[obj@net$prob==0]<-1
 obj@net$qval_BH<-array(p.adjust(as.vector(obj@net$pval_finite),'BH'),dim=dim(obj@net$pval),dimnames=dimnames(obj@net$pval))
 obj<-computeCommunProbPathway(obj);obj<-aggregateNet(obj)
 saveRDS(obj,file.path(dest,'CellChat.rds'),compress=TRUE)
 file.copy(file.path(row$source_dir,'input_cells.csv'),file.path(dest,'input_cells.csv'),overwrite=FALSE)
 for(f in c('group_eligibility.csv','status.csv'))if(file.exists(file.path(row$source_dir,f)))file.copy(file.path(row$source_dir,f),file.path(dest,f),overwrite=FALSE)
 audit<-data.frame(run=name,scope=row$scope,n_cells=nrow(input),n_tests=nrow(pr),n_permutations=10000,seed=seed,all_input_slots_preserved=TRUE,observed_probability_max_difference=delta,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
 write.csv(audit,file.path(dest,'input_and_model_audit.csv'),row.names=FALSE)
 writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
 jsonlite::write_json(list(status='Complete',run=name,n_permutations=10000,n_cells=nrow(input),completed=as.character(Sys.time())),file.path(dest,'COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
 message('COMPLETE ',name,' ',format(Sys.time()));rm(old,obj,pr,input);gc()
}
