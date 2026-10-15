args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(args[1])
rev <- file.path(project,'Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000')
out <- file.path(rev,'controlled_unannotated_20261012')
old <- file.path(project,'Cancers_ATRT_TME_Submission/Revision_20261007/analysis/data')
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(future)})
plan(sequential)
stopifnot(as.character(packageVersion('CellChat'))=='2.1.2')
meta <- read.csv(file.path(old,'GSE283839_metadata.csv'),check.names=FALSE,stringsAsFactors=FALSE)
counts <- as(readMM(file.path(old,'GSE283839_counts.mtx')),'dgCMatrix')
rownames(counts) <- readLines(file.path(old,'GSE283839_genes.tsv'))
colnames(counts) <- meta$Cell
stopifnot(nrow(counts)==36601,ncol(counts)==17564,!anyDuplicated(meta$Cell))
norm <- counts %*% Diagonal(x=1e4/Matrix::colSums(counts)); norm@x <- log1p(norm@x)
dimnames(norm) <- dimnames(counts)
rm(counts); gc()
states <- c('CP-like','Cilia-like','Hypoxic','IPC-like','Mesenchymal-like','NPC-like','OPC-like','RG-like')
candidates <- c('APP_CD74','NCAM1_FGFR1','SPP1_ITGAV_ITGB1','SPP1_CD44','PTN_PTPRZ1')
contributions <- list(); audits <- list(); screens <- list(); families <- list()
screen <- function(obj,subtype,scope) {
  ans <- list()
  for(gene in c('APP','CD74')) for(g in levels(obj@idents)) {
    a <- as.numeric(obj@data.signaling[gene,obj@idents==g]); b <- as.numeric(obj@data.signaling[gene,obj@idents!=g])
    fc <- log(mean(expm1(a))+1)-log(mean(expm1(b))+1)
    p <- suppressWarnings(wilcox.test(a,b)$p.value)
    pass <- max(round(mean(a>0),3),round(mean(b>0),3))>0 && fc>0 && p<.05
    ans[[length(ans)+1]] <- data.frame(subtype=subtype,scope=scope,gene=gene,group=g,n_group=length(a),n_other=length(b),pct_group=mean(a>0),pct_other=mean(b>0),logFC=fc,screening_P=p,screening_pass=pass)
  }
  do.call(rbind,ans)
}
for(sub in c('MYC','SHH','TYR')) {
  original <- file.path(rev,'analysis/cellchat',paste0('subtype__ATRT-',sub,'__state_resolved__exclude_MYC2'))
  base <- readRDS(file.path(original,'CellChat.rds'))
  m0 <- read.csv(file.path(original,'input_cells.csv'),check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(identical(m0$Cell,colnames(base@data)),identical(m0$Cell,rownames(base@meta)))
  stopifnot(max(abs(base@data-norm[,m0$Cell]))==0)
  m1 <- m0
  tumor_donors <- unique(m0$ID[m0$group=='Tumor'])
  for(donor in tumor_donors) {
    ix <- which(m0$ID==donor & m0$group=='Tumor'); target <- length(ix)
    pool <- meta[meta$ID==donor & meta$Final_Annotation %in% c(states,'Unannotated'),,drop=FALSE]
    stopifnot(nrow(pool)>=target,donor!='ATRT-MYC-2')
    set.seed(42); chosen <- pool[sample.int(nrow(pool),target),,drop=FALSE]
    for(col in names(meta)) m1[ix,col] <- chosen[[col]]
    m1$group[ix] <- 'Tumor'; m1$samples[ix] <- donor
  }
  stopifnot(!anyDuplicated(m1$Cell),identical(m0$ID,m1$ID),identical(m0$group,m1$group))
  glia <- which(m0$group!='Tumor')
  stopifnot(identical(m0[glia,],m1[glia,]),!any(m1$Final_Annotation=='Immune-like'))
  stopifnot(all(m1$Final_Annotation[m1$group=='Tumor'] %in% c(states,'Unannotated')))
  for(scope in c('restricted','plus_unannotated')) {
    m <- if(scope=='restricted')m0 else m1
    name <- paste(sub,scope,sep='__'); path <- file.path(out,'models',name)
    dir.create(path,recursive=TRUE,showWarnings=FALSE)
    stopifnot(!file.exists(file.path(path,'prepared.rds')))
    rownames(m) <- m$Cell; m$group <- factor(m$group,levels=levels(base@idents));m$samples <- factor(m$ID)
    obj <- createCellChat(object=norm[,m$Cell,drop=FALSE],meta=m,group.by='group')
    obj@DB <- base@DB; obj <- subsetData(obj)
    stopifnot(identical(rownames(obj@data.signaling),rownames(base@data.signaling)))
    # All database interactions with complete gene components represented in the
    # matrix are fixed a priori; differential-expression significance is not a gate.
    obj <- identifyOverExpressedInteractions(obj,features=rownames(obj@data.signaling))
    families[[name]] <- rownames(obj@LR$LRsig)
    stopifnot(all(candidates %in% families[[name]]))
    screens[[name]] <- screen(obj,sub,scope)
    write.csv(m,file.path(path,'input_cells.csv'),row.names=FALSE)
    write.csv(obj@LR$LRsig,file.path(path,'testing_family.csv'),row.names=FALSE)
    saveRDS(obj,file.path(path,'prepared.rds'),compress=TRUE)
    tab <- aggregate(list(n_cells=rep(1,nrow(m))),list(subtype=rep(sub,nrow(m)),scope=rep(scope,nrow(m)),patient=m$ID,group=m$group),sum)
    contributions[[name]] <- tab
    audits[[name]] <- data.frame(run=name,n_cells=nrow(m),n_tumor=sum(m$group=='Tumor'),n_unannotated=sum(m$Final_Annotation=='Unannotated'),n_immune_like=sum(m$Final_Annotation=='Immune-like'),n_LR=nrow(obj@LR$LRsig),n_tests=nrow(obj@LR$LRsig)*nlevels(obj@idents)^2,patient_group_order_fixed=identical(m0$ID,as.character(m$ID))&&identical(m0$group,as.character(m$group)),glial_cells_fixed=identical(m0$Cell[glia],m$Cell[glia]),source_annotations_preserved=all(m$Final_Annotation==meta$Final_Annotation[match(m$Cell,meta$Cell)]))
    message('PREPARED ',name,' nuclei=',nrow(m),' LR=',nrow(obj@LR$LRsig))
  }
  # Diagnose the historical gene-selection gate without changing old outputs.
  for(scope in c('state_resolved','all_author_tumor')) {
    h <- readRDS(file.path(rev,'analysis/cellchat',paste0('subtype__ATRT-',sub,'__',scope,'__exclude_MYC2'),'CellChat.rds'))
    s <- screen(h,sub,paste0('historical_',scope))
    info <- h@var.features$features.info
    expected <- vapply(seq_len(nrow(s)),function(i)any(info$features==s$gene[i]&info$clusters==s$group[i]),logical(1))
    stopifnot(identical(expected,s$screening_pass))
    screens[[paste(sub,scope,sep='__')]] <- s
  }
}
stopifnot(all(vapply(families,function(x)identical(x,families[[1]]),logical(1))))
write.csv(do.call(rbind,contributions),file.path(out,'tables/S13_controlled_donor_contributions.csv'),row.names=FALSE)
write.csv(do.call(rbind,audits),file.path(out,'tables/S13_controlled_input_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,screens),file.path(out,'tables/S13_APP_CD74_gene_screen_audit.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'data/sessionInfo.txt'))
writeLines(deparse(CellChat::identifyOverExpressedGenes),file.path(out,'code/identifyOverExpressedGenes_installed.R'))
writeLines(deparse(CellChat::identifyOverExpressedInteractions),file.path(out,'code/identifyOverExpressedInteractions_installed.R'))
message('ALL SIX MODELS PREPARED; common LR family=',length(families[[1]]))
