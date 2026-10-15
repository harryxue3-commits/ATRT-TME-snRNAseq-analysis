#!/usr/bin/env Rscript
# CellChat 2.1.2 filterCommunication uses <= min.cells. Restore the intended >=N
# boundary for two exactly-ten endothelial groups without rerunning or changing
# the 1,000-permutation null. The function does not alter pval during filtering.
args<-commandArgs(trailingOnly=TRUE);out<-normalizePath(args[1])
suppressPackageStartupMessages({library(CellChat);library(future)})
plan(sequential)
stopifnot(as.character(packageVersion('CellChat'))=='2.1.2')
audit<-list()
libs<-if(length(args)>=2)strsplit(args[2],',')[[1]] else c('SCPCL000811','SCPCL000812')
for(lib in libs) {
 name<-paste0('ALSF__',lib,'__conservative__min10');d<-file.path(out,name)
 stopifnot(read.csv(file.path(d,'status.csv'))$status=='Complete')
 obj<-readRDS(file.path(d,'CellChat.rds'))
 if(file.exists(file.path(d,'inclusive_filter_correction.csv'))) {
  audit[[lib]]<-read.csv(file.path(d,'inclusive_filter_correction.csv'));next
 }
 counts<-table(obj@idents);stopifnot(any(counts==10))
 origP<-obj@net$pval;origProb<-obj@net$prob;options.original<-obj@options
 # Observed probabilities are deterministic; bootstrap count only affects pval.
 # Keep the stored 1,000-permutation pval array and restore it after this call.
 observed<-computeCommunProb(obj,type='truncatedMean',trim=.1,raw.use=TRUE,population.size=FALSE,nboot=1,seed.use=20261008)
 stopifnot(identical(dimnames(observed@net$prob),dimnames(origProb)))
 keep<-which(counts>10)
 discrepancy<-max(abs(observed@net$prob[keep,keep,,drop=FALSE]-origProb[keep,keep,,drop=FALSE]))
 stopifnot(discrepancy<1e-12)
 obj@net$prob<-observed@net$prob;obj@net$pval<-origP;obj@options<-options.original
 obj<-filterCommunication(obj,min.cells=9)
 stopifnot(identical(obj@net$pval,origP))
 pr<-as.data.frame.table(obj@net$prob,responseName='probability',stringsAsFactors=FALSE);names(pr)[1:3]<-c('source','target','interaction_name')
 pr$p_native<-as.vector(obj@net$pval);pr$permutation_exceedances<-as.integer(round(pr$p_native*1000))
 pr$p_plus1<-(pr$permutation_exceedances+1)/1001;pr$q_BH_within_object<-p.adjust(pr$p_plus1,'BH')
 pr$nominal_p_lt_05<-pr$p_native<.05 & pr$probability>0;pr$BH_q_lt_05<-pr$q_BH_within_object<.05 & pr$probability>0
 pr$run<-name;info<-obj@DB$interaction;info$interaction_name<-rownames(info);pr<-merge(pr,info,by='interaction_name',all.x=TRUE,sort=FALSE)
 write.csv(pr,file.path(d,'all_tested_interactions.csv'),row.names=FALSE)
 write.csv(pr[pr$nominal_p_lt_05,],file.path(d,'nominal_interactions.csv'),row.names=FALSE)
 obj<-computeCommunProbPathway(obj);obj<-aggregateNet(obj)
 saveRDS(obj,file.path(d,'CellChat.rds'),compress=TRUE)
 row<-data.frame(run=name,corrected_group=paste(names(counts)[counts==10],collapse=';'),group_count=10,reported_inclusive_minimum=10,CellChat_internal_min_cells=9,unaffected_probability_max_difference=discrepancy,original_1000_permutation_p_values_preserved=TRUE,method='Restore deterministic observed probabilities; retain original 1000-permutation pval array unchanged')
 write.csv(row,file.path(d,'inclusive_filter_correction.csv'),row.names=FALSE);audit[[lib]]<-row
}
write.csv(do.call(rbind,audit),file.path(dirname(out),'tables/S18_inclusive_filter_correction.csv'),row.names=FALSE)
writeLines(capture.output(filterCommunication),file.path(dirname(out),'code/filterCommunication_source.R'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo_inclusive_filter_correction.txt'))
