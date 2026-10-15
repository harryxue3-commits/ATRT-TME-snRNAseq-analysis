#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]); out <- file.path(root,'cellchat'); dir.create(out,showWarnings=FALSE)
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(future)})
plan(sequential);options(future.globals.maxSize=8*1024^3)
set.seed(42)
stopifnot(as.character(packageVersion('CellChat'))=='2.1.2')
counts <- as(readMM(file.path(root,'data/GSE283839_counts.mtx')),'dgCMatrix')
genes <- readLines(file.path(root,'data/GSE283839_genes.tsv'))
meta <- read.csv(file.path(root,'data/GSE283839_metadata.csv'),check.names=FALSE,stringsAsFactors=FALSE)
rownames(counts)<-genes;colnames(counts)<-meta$Cell;rownames(meta)<-meta$Cell
stopifnot(nrow(counts)==36601,ncol(counts)==17564,!anyDuplicated(meta$Cell))
norm <- counts %*% Diagonal(x=1e4/Matrix::colSums(counts));norm@x<-log1p(norm@x);dimnames(norm)<-dimnames(counts)
data(CellChatDB.human,package='CellChat');db<-CellChatDB.human
saveRDS(db,file.path(out,'CellChatDB_human_v2.1.2.rds'))
write.csv(db$interaction,file.path(out,'database_interactions.csv'),row.names=TRUE)
write.csv(db$complex,file.path(out,'database_complexes.csv'),row.names=TRUE)
states<-c('CP-like','Cilia-like','Hypoxic','IPC-like','Mesenchymal-like','NPC-like','OPC-like','RG-like')
meta$group<-ifelse(meta$Final_Annotation %in% states,'Tumor',ifelse(meta$Final_Annotation=='Astrocytes','Astrocytes',ifelse(meta$Final_Annotation=='Microglia','Microglia',NA)))
contributions<-list();statuses<-list();allcells<-list()
run_one<-function(name,m,population.size=FALSE,min_glia=10) {
  message('START ',name,' ',format(Sys.time())); path<-file.path(out,name);dir.create(path,showWarnings=FALSE)
  m$group<-droplevels(factor(m$group));rownames(m)<-m$Cell
  m$samples<-factor(m$ID)
  write.csv(m,file.path(path,'input_cells.csv'),row.names=FALSE)
  if(!nrow(m)) {
    contributions[[name]]<<-data.frame(ID=character(),group=character(),n_cells=integer(),run=character())
    statuses[[name]]<<-data.frame(run=name,status='Skipped: no eligible groups',n_cells=0)
    return(invisible(NULL))
  }
  tab<-as.data.frame(table(m$ID,m$group),stringsAsFactors=FALSE);names(tab)<-c('ID','group','n_cells');tab<-tab[tab$n_cells>0,];tab$run<-name
  contributions[[name]]<<-tab
  if (nlevels(m$group)<2) {statuses[[name]]<<-data.frame(run=name,status='Skipped: fewer than two groups',n_cells=nrow(m));return(invisible(NULL))}
  if(file.exists(file.path(path,'CellChat.rds'))) {statuses[[name]]<<-data.frame(run=name,status='Complete (existing saved run)',n_cells=nrow(m));return(invisible(NULL))}
  tryCatch({
    set.seed(42);obj<-createCellChat(object=norm[,m$Cell,drop=FALSE],meta=m,group.by='group');obj@DB<-db
    obj<-subsetData(obj);obj<-identifyOverExpressedGenes(obj,do.fast=FALSE);obj<-identifyOverExpressedInteractions(obj)
    obj<-computeCommunProb(obj,type='truncatedMean',trim=.1,raw.use=TRUE,population.size=population.size,nboot=100,seed.use=42)
    # The original five-nucleus astrocyte exception is explicit; other glia default to ten.
    sizes<-table(obj@idents)
    for(g in names(sizes)) {
      threshold<-if(g=='Astrocytes') 5 else if(g=='Tumor') 10 else min_glia
      if(sizes[[g]]<threshold) {obj@net$prob[g,,]<-0;obj@net$prob[,g,]<-0;obj@net$pval[g,,]<-1;obj@net$pval[,g,]<-1}
    }
    prob<-as.data.frame.table(obj@net$prob,responseName='probability',stringsAsFactors=FALSE);names(prob)[1:3]<-c('source','target','interaction_name')
    prob$p_value<-as.vector(obj@net$pval);prob$nominal_p_lt_05<-prob$p_value<.05 & prob$probability>0
    prob$q_BH_within_object<-p.adjust(prob$p_value,'BH');prob$run<-name
    info<-db$interaction;info$interaction_name<-rownames(info)
    prob<-merge(prob,info,by='interaction_name',all.x=TRUE,sort=FALSE)
    write.csv(prob,file.path(path,'all_tested_interactions.csv'),row.names=FALSE)
    write.csv(prob[prob$nominal_p_lt_05,],file.path(path,'nominal_interactions.csv'),row.names=FALSE)
    obj<-computeCommunProbPathway(obj);obj<-aggregateNet(obj)
    saveRDS(obj,file.path(path,'CellChat.rds'),compress=TRUE)
    statuses[[name]]<<-data.frame(run=name,status='Complete',n_cells=nrow(m))
    rm(obj);gc()
  },error=function(e){writeLines(conditionMessage(e),file.path(path,'ERROR.txt'));statuses[[name]]<<-data.frame(run=name,status=paste('ERROR',conditionMessage(e)),n_cells=nrow(m));message(conditionMessage(e))})
  write.csv(do.call(rbind,statuses),file.path(out,'run_status.csv'),row.names=FALSE)
  write.csv(do.call(rbind,contributions),file.path(out,'donor_contributions.csv'),row.names=FALSE)
}
select_subtype<-function(sub,scope,include_MYC2=FALSE) {
  m<-meta[meta$subtype==sub & (include_MYC2|meta$ID!='ATRT-MYC-2'),]
  if(scope=='all_author_tumor') m$group[m$Final_Annotation_Focus!='TME']<-'Tumor'
  selected<-list()
  for(g in c('Tumor','Microglia','Astrocytes')) {
    d<-m[!is.na(m$group)&m$group==g,,drop=FALSE];tab<-table(d$ID)
    threshold<-switch(g,Tumor=150,Microglia=10,Astrocytes=5)
    eligible<-names(tab)[tab>=threshold]
    if(!length(eligible))next
    target<-min(tab[eligible]);
    for(patient in eligible) {
      v<-d[d$ID==patient,,drop=FALSE]
      if(g!='Astrocytes') {set.seed(42);v<-v[sample.int(nrow(v),target),,drop=FALSE]}
      selected[[paste(g,patient)]]<-v
    }
  }
  do.call(rbind,selected)
}
for(scope in c('state_resolved','all_author_tumor')) {
  for(sub in c('ATRT-MYC','ATRT-SHH','ATRT-TYR')) {
    run_one(paste('subtype',sub,scope,'exclude_MYC2',sep='__'),select_subtype(sub,scope,FALSE))
  }
  run_one(paste('subtype','ATRT-MYC',scope,'include_MYC2',sep='__'),select_subtype('ATRT-MYC',scope,TRUE))
}
# Same-donor networks retain every qualifying nucleus, without between-donor pooling.
for(patient in sort(unique(meta$ID))) {
  m<-meta[meta$ID==patient & !is.na(meta$group),,drop=FALSE]
  tab<-table(m$group);allowed<-names(tab)[tab>=ifelse(names(tab)=='Astrocytes',5,ifelse(names(tab)=='Tumor',20,10))]
  m<-m[m$group%in%allowed,,drop=FALSE]
  run_one(paste('donor',patient,'state_resolved',sep='__'),m)
}
# Conservative identity sensitivity, restricted to source TME annotations only.
prior<-file.path(dirname(dirname(root)),'Revision_20261005/analysis/data/analysis_sets.json')
sets<-jsonlite::fromJSON(prior)
for(patient in sort(unique(meta$ID))) {
  m<-meta[meta$ID==patient,,drop=FALSE];m$group<-NA_character_
  m$group[m$Final_Annotation%in%states]<-'Tumor'
  m$group[m$Cell%in%sets[['Conservative astrocytes']]]<-'Astrocytes'
  m$group[m$Cell%in%sets[['Conservative microglia']]]<-'Microglia'
  m$group[m$Cell%in%sets[['Macrophage-like including BAM']]]<-'Macrophage_BAM_like'
  m<-m[!is.na(m$group),,drop=FALSE];tab<-table(m$group)
  allowed<-names(tab)[tab>=ifelse(names(tab)=='Tumor',20,5)]
  run_one(paste('conservative',patient,'min5',sep='__'),m[m$group%in%allowed,,drop=FALSE],min_glia=5)
}
write.csv(do.call(rbind,statuses),file.path(out,'run_status.csv'),row.names=FALSE)
write.csv(do.call(rbind,contributions),file.path(out,'donor_contributions.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(c('Fresh sensitivity fits, not exact reproductions of unavailable original objects.','CellChat 2.1.2, bundled full human database; natural log(1+CP10k); truncated mean trim=0.1; nboot=100; seed=42; population.size=FALSE.','The original molecular subtype analysis pooled unequal donor sets across sender and receiver groups. Exact input IDs and donor counts are exported.','Nominal P<0.05 is shown for comparability, with separate BH-adjusted values; neither is evidence of replicated subtype specificity.','All source malignant assignments are preserved. Expanded tumor includes source Unannotated and Immune-like tumor nuclei.'),file.path(out,'METHODS.txt'))
if(any(grepl('ERROR',do.call(rbind,statuses)$status)))quit(status=1)
