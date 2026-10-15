#!/usr/bin/env Rscript
# Independent sample networks. Paired specimens remain linked by patient_id.
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(args[1]); out <- normalizePath(args[2]); parity <- as.integer(args[3])
input <- file.path(project,'ALSF-ATRT-snRNAseq/analysis_20261006/cellchat')
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(future)})
plan(sequential);options(future.globals.maxSize=12*1024^3)
stopifnot(as.character(packageVersion('CellChat'))=='2.1.2')
data(CellChatDB.human,package='CellChat');db<-CellChatDB.human
if(parity==0) {
 saveRDS(db,file.path(out,'CellChatDB_human_v2.1.2.rds'))
 write.csv(db$interaction,file.path(out,'database_interactions.csv'),row.names=TRUE)
 write.csv(db$complex,file.path(out,'database_complexes.csv'),row.names=TRUE)
}
manifest<-read.csv(file.path(input,'input_manifest.csv'),stringsAsFactors=FALSE)
specs<-data.frame(tier=c('conservative','broad','conservative'),minimum=c(20,20,10))
# Optional execution-only selector; the statistical settings are identical.
if(length(args)>=5)specs<-specs[paste0(specs$tier,specs$minimum)%in%strsplit(args[5],',')[[1]],,drop=FALSE]
for(j in seq_len(nrow(manifest))) {
 if(length(args)>=4) { if(!manifest$library_id[j]%in%strsplit(args[4],',')[[1]])next } else { if((j%%2)!=parity)next }
 lib<-manifest$library_id[j];p<-file.path(input,lib)
 x<-as(readMM(gzfile(file.path(p,'counts.mtx.gz'))),'dgCMatrix')
 genes<-readLines(file.path(p,'genes.tsv'));cells<-readLines(file.path(p,'barcodes.tsv'))
 meta<-read.csv(file.path(p,'metadata.csv'),stringsAsFactors=FALSE,check.names=FALSE)
 stopifnot(identical(meta$cell_id,cells),!anyDuplicated(cells),all(x@x==round(x@x)),sum(x)==manifest$total_UMI_export[j])
 dimnames(x)<-list(genes,cells);rownames(meta)<-meta$cell_id
 norm<-x%*%Diagonal(x=1e4/colSums(x));norm@x<-log1p(norm@x);dimnames(norm)<-dimnames(x)
 for(i in seq_len(nrow(specs))) {
  tier<-specs$tier[i];minimum<-specs$minimum[i]
  name<-paste('ALSF',lib,tier,paste0('min',minimum),sep='__');d<-file.path(out,name);dir.create(d,recursive=TRUE,showWarnings=FALSE)
  m<-meta;m$group<-m[[paste0('cellchat_',tier)]];counts<-table(m$group[m$group!='Exclude'])
  eligible<-names(counts)[counts>=minimum]
  write.csv(data.frame(group=names(counts),n_cells=as.integer(counts),eligible=names(counts)%in%eligible),file.path(d,'group_eligibility.csv'),row.names=FALSE)
  m<-m[m$group%in%eligible,,drop=FALSE];m$group<-droplevels(factor(m$group));m$samples<-factor(m$sample_id)
  write.csv(m,file.path(d,'input_cells.csv'),row.names=FALSE)
  if(nlevels(m$group)<2 || !'Tumor_CNV_supported'%in%eligible) {
   write.csv(data.frame(run=name,status='Skipped: insufficient groups or no CNV-supported tumor',n_cells=nrow(m),n_groups=nlevels(m$group)),file.path(d,'status.csv'),row.names=FALSE);next
  }
  if(file.exists(file.path(d,'CellChat.rds')))next
  message('START ',name,' ',format(Sys.time()),' cells=',nrow(m),' groups=',nlevels(m$group))
  tryCatch({
   set.seed(20261008)
   obj<-createCellChat(object=norm[,m$cell_id,drop=FALSE],meta=m,group.by='group');obj@DB<-db
   obj<-subsetData(obj);obj<-identifyOverExpressedGenes(obj,do.fast=FALSE);obj<-identifyOverExpressedInteractions(obj)
   obj<-computeCommunProb(obj,type='truncatedMean',trim=.1,raw.use=TRUE,population.size=FALSE,nboot=1000,seed.use=20261008)
   obj<-filterCommunication(obj,min.cells=minimum-1)
   pr<-as.data.frame.table(obj@net$prob,responseName='probability',stringsAsFactors=FALSE);names(pr)[1:3]<-c('source','target','interaction_name')
   pr$p_native<-as.vector(obj@net$pval)
   pr$permutation_exceedances<-as.integer(round(pr$p_native*1000))
   pr$p_plus1<-(pr$permutation_exceedances+1)/1001
   pr$q_BH_within_object<-p.adjust(pr$p_plus1,'BH')
   pr$nominal_p_lt_05<-pr$p_native<.05 & pr$probability>0
   pr$BH_q_lt_05<-pr$q_BH_within_object<.05 & pr$probability>0
   pr$run<-name;info<-db$interaction;info$interaction_name<-rownames(info)
   pr<-merge(pr,info,by='interaction_name',all.x=TRUE,sort=FALSE)
   write.csv(pr,file.path(d,'all_tested_interactions.csv'),row.names=FALSE)
   write.csv(pr[pr$nominal_p_lt_05,],file.path(d,'nominal_interactions.csv'),row.names=FALSE)
   obj<-computeCommunProbPathway(obj);obj<-aggregateNet(obj);saveRDS(obj,file.path(d,'CellChat.rds'),compress=TRUE)
   write.csv(data.frame(run=name,status='Complete',n_cells=nrow(m),n_groups=nlevels(m$group)),file.path(d,'status.csv'),row.names=FALSE)
   message('COMPLETE ',name,' ',format(Sys.time()));rm(obj);gc()
  },error=function(e){writeLines(conditionMessage(e),file.path(d,'ERROR.txt'));write.csv(data.frame(run=name,status=paste('ERROR',conditionMessage(e)),n_cells=nrow(m),n_groups=nlevels(m$group)),file.path(d,'status.csv'),row.names=FALSE);message('ERROR ',name,' ',conditionMessage(e))})
 }
 rm(x,norm);gc()
}
writeLines(capture.output(sessionInfo()),file.path(out,paste0('sessionInfo_worker',parity,'.txt')))
