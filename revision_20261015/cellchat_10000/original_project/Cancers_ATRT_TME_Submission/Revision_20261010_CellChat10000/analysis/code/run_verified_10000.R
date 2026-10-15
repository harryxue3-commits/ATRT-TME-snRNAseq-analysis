#!/usr/bin/env Rscript
# Repeat only computeCommunProb on frozen CellChat 2.1.2 subtype objects.
# Original cells, labels, input expression, LR selection, DB and seed are retained.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3, args[2] %in% c('MYC','SHH','TYR'), args[3] %in% c('state_resolved','all_author_tumor'))
project <- normalizePath(args[1]); subtype <- args[2]
rev <- file.path(project,'Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000')
old <- file.path(project,'Cancers_ATRT_TME_Submission/Revision_20261007/analysis/cellchat')
name <- paste0('subtype__ATRT-',subtype,'__',args[3],'__exclude_MYC2')
out <- file.path(rev,'analysis/cellchat',name);dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(future)})
stopifnot(as.character(packageVersion('CellChat'))=='2.1.2')
plan(sequential);options(future.globals.maxSize=8*1024^3)
if(file.exists(file.path(out,'COMPLETE.json')))stop('Already complete: refusing to replace verified results')
source <- readRDS(file.path(old,name,'CellChat.rds'))
input <- read.csv(file.path(old,name,'input_cells.csv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(identical(colnames(source@data.signaling),input$Cell),identical(rownames(source@meta),input$Cell))
stopifnot(identical(as.character(source@idents),input$group),all(source@meta$ID==input$ID))
stopifnot(source@options$datatype=='RNA',source@options$parameter$nboot==100,
          source@options$parameter$seed.use==42,source@options$parameter$trim==0.1,
          identical(source@options$parameter$population.size,FALSE))
file.copy(file.path(old,name,'input_cells.csv'),file.path(out,'input_cells.csv'),overwrite=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(deparse(CellChat::computeCommunProb),file.path(out,'computeCommunProb_source.R'))
writeLines(as.character(Sys.getpid()),file.path(out,'worker_pid.txt'))
message('START ',name,' ',format(Sys.time()),' nuclei=',nrow(input),' permutations=10000')
started <- Sys.time()
eq <- read.csv(file.path(rev,'analysis/tables/S13_acceleration_equivalence_validation.csv'))
stopifnot(nrow(eq)==7,all(eq$all_entries_identical))
source(file.path(rev,'analysis/code/validated_fast_cellchat.R'))
fn <- make_fast_cellchat(file.path(rev,'analysis/code'))
writeLines(attr(fn,'acceleration_source'),file.path(out,'computeCommunProb_accelerated_source.R'))
obj <- fn(source,type='truncatedMean',trim=.1,raw.use=TRUE,
                         population.size=FALSE,nboot=10000,seed.use=42)
stopifnot(identical(source@meta,obj@meta),identical(source@idents,obj@idents),
          identical(source@data,obj@data),identical(source@data.signaling,obj@data.signaling),
          identical(source@DB,obj@DB),identical(source@LR,obj@LR),
          identical(dimnames(source@net$prob),dimnames(obj@net$prob)))
delta <- max(abs(source@net$prob-obj@net$prob))
stopifnot(delta < 1e-12)
prob <- as.data.frame.table(obj@net$prob,responseName='probability',stringsAsFactors=FALSE)
names(prob)[1:3]<-c('source','target','interaction_name')
prob$p_value <- as.vector(obj@net$pval)
prob$permutation_exceedances <- as.integer(round(prob$p_value*10000))
prob$p_value_plus1 <- (prob$permutation_exceedances+1)/10001
prob$q_BH_plus1_within_object <- p.adjust(prob$p_value_plus1,'BH')
prob$nominal_p_lt_05 <- prob$p_value < .05 & prob$probability > 0
prob$BH_q_lt_05 <- prob$q_BH_plus1_within_object < .05 & prob$probability > 0
prob$n_permutations <- 10000L;prob$seed <- 42L;prob$run<-name
info<-obj@DB$interaction;info$interaction_name<-rownames(info)
prob<-merge(prob,info,by='interaction_name',all.x=TRUE,sort=FALSE)
write.csv(prob,file.path(out,'all_tested_interactions.csv'),row.names=FALSE)
write.csv(prob[prob$nominal_p_lt_05,],file.path(out,'nominal_interactions.csv'),row.names=FALSE)
write.csv(prob[prob$BH_q_lt_05,],file.path(out,'BH_supported_interactions.csv'),row.names=FALSE)
obj@options$null_calculation_implementation <- 'Validated computational acceleration; native CellChat RNG, equations, helpers and strict exceedance rule unchanged'
obj <- computeCommunProbPathway(obj);obj<-aggregateNet(obj)
saveRDS(obj,file.path(out,'CellChat.rds'),compress=TRUE)
audit<-data.frame(run=name,n_cells=nrow(input),n_tests=nrow(prob),n_permutations=10000,seed=42,
                 cell_order_preserved=TRUE,annotations_preserved=TRUE,patient_IDs_preserved=TRUE,
                 expression_preserved=TRUE,database_preserved=TRUE,LR_selection_preserved=TRUE,
                 observed_probability_max_difference=delta,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
write.csv(audit,file.path(out,'input_and_model_audit.csv'),row.names=FALSE)
jsonlite::write_json(list(status='Complete',run=name,n_cells=nrow(input),n_tests=nrow(prob),
                         n_permutations=10000,seed=42,completed=as.character(Sys.time())),
                    file.path(out,'COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
message('COMPLETE ',name,' ',format(Sys.time()))
