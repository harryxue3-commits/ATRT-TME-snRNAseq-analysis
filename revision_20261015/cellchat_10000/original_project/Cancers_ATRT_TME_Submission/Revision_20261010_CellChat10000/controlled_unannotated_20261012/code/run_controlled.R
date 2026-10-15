args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(args[1]); sub <- args[2]
stopifnot(sub %in% c('MYC','SHH','TYR'))
rev <- file.path(project,'Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000')
out <- file.path(rev,'controlled_unannotated_20261012')
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(future)})
plan(sequential); options(future.globals.maxSize=8*1024^3)
source(file.path(rev,'analysis/code/validated_fast_cellchat.R'))
fn <- make_fast_cellchat(file.path(rev,'analysis/code'))
for(scope in c('restricted','plus_unannotated')) {
  name <- paste(sub,scope,sep='__');path <- file.path(out,'models',name)
  stopifnot(!file.exists(file.path(path,'COMPLETE.json')))
  obj <- readRDS(file.path(path,'prepared.rds')); original <- obj
  started <- Sys.time();message('START ',name,' ',started)
  obj <- fn(obj,type='truncatedMean',trim=.1,raw.use=TRUE,population.size=FALSE,nboot=10000,seed.use=42)
  stopifnot(identical(original@meta,obj@meta),identical(original@idents,obj@idents),identical(original@data,obj@data),identical(original@DB,obj@DB),identical(original@LR,obj@LR))
  if(scope=='restricted') {
    old <- readRDS(file.path(rev,'analysis/cellchat',paste0('subtype__ATRT-',sub,'__state_resolved__exclude_MYC2'),'CellChat.rds'))
    common <- dimnames(old@net$prob)[[3]]
    stopifnot(identical(obj@net$prob[,,common],old@net$prob),identical(obj@net$pval[,,common],old@net$pval))
  }
  d <- as.data.frame.table(obj@net$prob,responseName='probability',stringsAsFactors=FALSE)
  names(d)[1:3] <- c('source','target','interaction_name')
  d$p_native <- as.vector(obj@net$pval)
  d$permutation_exceedances <- as.integer(round(d$p_native*10000))
  d$p_finite <- (d$permutation_exceedances+1)/10001
  d$p_finite[d$probability==0] <- 1
  d$q_BH <- p.adjust(d$p_finite,'BH')
  d$nominal_finite <- d$p_finite<.05 & d$probability>0
  d$nominal_native <- d$p_native<.05 & d$probability>0
  d$BH_supported <- d$q_BH<.05 & d$probability>0
  d$subtype <- sub; d$scope <- scope;d$n_permutations <- 10000L
  info <- obj@DB$interaction;info$interaction_name <- rownames(info)
  d <- merge(d,info,by='interaction_name',all.x=TRUE,sort=FALSE)
  write.csv(d,file.path(path,'all_tested_interactions.csv'),row.names=FALSE)
  write.csv(d[d$BH_supported,],file.path(path,'BH_supported_interactions.csv'),row.names=FALSE)
  obj@options$controlled_comparison <- 'Fixed original donors, group sizes, glial cells and complete represented-gene DB family; plus Unannotated only; no Immune-like cells'
  obj@options$null_calculation_implementation <- 'Validated acceleration retaining native CellChat shuffle indices and score equations'
  obj <- computeCommunProbPathway(obj);obj <- aggregateNet(obj)
  saveRDS(obj,file.path(path,'CellChat.rds'),compress=TRUE)
  jsonlite::write_json(list(status='Complete',run=name,n_cells=ncol(obj@data),n_tests=nrow(d),n_LR=nrow(obj@LR$LRsig),n_permutations=10000,seed=42,input_metadata_preserved=TRUE,baseline_common_entries_exact=if(scope=='restricted')TRUE else NA,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs'))),file.path(path,'COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
  message('COMPLETE ',name,' ',Sys.time())
  rm(obj,original);gc()
}
