# Computational acceleration only; original permutation generation, observed
# scores, LR/complex/cofactor helpers, strict exceedance test, and outputs kept.
make_fast_cellchat <- function(code_dir) {
  Rcpp::sourceCpp(file.path(code_dir,'fast_bootstrap_means.cpp'),cacheDir='/private/tmp/atrt_cellchat_rcpp_cache',showOutput=FALSE)
  src<-deparse(CellChat::computeCommunProb)
  start<-grep('^    data.use.avg.boot <-',src)
  end<-grep('^    pb <- txtProgressBar',src)-1L
  stopifnot(length(start)==1,length(end)==1)
  src<-c(src[seq_len(start-1)],
         '    stopifnot(object@options$datatype == "RNA", !population.size, type == "truncatedMean")',
         '    data.use.avg.wide <- fast_bootstrap_means(data.use, as.integer(group), permutation, trim)',
         src[seq(end+1,length(src))])
  start<-grep('^            Pboot <- sapply',src)
  end<-grep('^            nReject <-',src)-1L
  stopifnot(length(start)==1,length(end)==1)
  block<-c(
    '            dataLavgB <- as.vector(computeExpr_LR(geneL[i], data.use.avg.wide, complex_input))',
    '            dataRavgB <- as.vector(computeExpr_LR(geneR[i], data.use.avg.wide, complex_input))',
    '            coA <- as.vector(computeExpr_coreceptor(cofactor_input, data.use.avg.wide, pairLRsig[i,,drop=FALSE], type="A"))',
    '            coI <- as.vector(computeExpr_coreceptor(cofactor_input, data.use.avg.wide, pairLRsig[i,,drop=FALSE], type="I"))',
    '            dataRavgB <- dataRavgB * coA / coI',
    '            agon <- if(i %in% index.agonist) as.vector(computeExpr_agonist(data.use.avg.wide, pairLRsig, cofactor_input, index.agonist=i, Kh=Kh, n=n)) else rep(1, numCluster*nboot)',
    '            antag <- if(i %in% index.antagonist) as.vector(computeExpr_antagonist(data.use.avg.wide, pairLRsig, cofactor_input, index.antagonist=i, Kh=Kh, n=n)) else rep(1, numCluster*nboot)',
    '            Pboot <- matrix(0,numCluster*numCluster,nboot)',
    '            for(tt in seq_len(numCluster)) for(ss in seq_len(numCluster)) {',
    '              si <- seq.int(ss,by=numCluster,length.out=nboot); ti <- seq.int(tt,by=numCluster,length.out=nboot)',
    '              lr <- dataLavgB[si]*dataRavgB[ti]',
    '              p1 <- lr^n/(Kh^n+lr^n)',
    '              p2 <- agon[si]*agon[ti]; p3 <- antag[si]*antag[ti]',
    '              Pboot[ss+(tt-1L)*numCluster,] <- p1*p2*p3*P.spatial[ss,tt]',
    '            }',
    # Re-evaluate any comparison close to the observed score with the original
    # R aggregate/mean and per-permutation helper calls. The 1e-12 threshold only
    # selects calculations to redo; the FINAL test remains strictly > 0.
    '            near <- which(abs(Pboot - Pnull) <= 1e-12 & Pnull > 0, arr.ind=TRUE)',
    '            if(nrow(near)) {',
    '              keys <- c(geneL[i],geneR[i]); needed <- intersect(keys,rownames(data.use))',
    '              cx <- setdiff(keys,rownames(data.use)); if(length(cx)) needed <- c(needed,unlist(complex_input[cx,grepl("subunit",colnames(complex_input)),drop=FALSE],use.names=FALSE))',
    '              cf <- unlist(pairLRsig[i,c("co_A_receptor","co_I_receptor","agonist","antagonist")],use.names=FALSE); cf<-unique(cf[!is.na(cf)&cf!=""])',
    '              if(length(cf))needed<-c(needed,unlist(cofactor_input[cf,grepl("cofactor",colnames(cofactor_input)),drop=FALSE],use.names=FALSE))',
    '              needed<-intersect(unique(needed),rownames(data.use))',
    '              for(nE in unique(near[,2])) {',
    '                avg <- aggregate(t(data.use[needed,,drop=FALSE]),list(group[permutation[,nE]]),FUN=FunMean); avg<-t(avg[,-1,drop=FALSE])',
    '                l <- computeExpr_LR(geneL[i],avg,complex_input); r<-computeExpr_LR(geneR[i],avg,complex_input)',
    '                ca<-computeExpr_coreceptor(cofactor_input,avg,pairLRsig[i,,drop=FALSE],type="A");ci<-computeExpr_coreceptor(cofactor_input,avg,pairLRsig[i,,drop=FALSE],type="I")',
    '                r<-r*ca/ci; lr<-Matrix::crossprod(l,r); p1<-lr^n/(Kh^n+lr^n)',
    '                p2<-if(i %in% index.agonist)Matrix::crossprod(matrix(computeExpr_agonist(avg,pairLRsig,cofactor_input,index.agonist=i,Kh=Kh,n=n),nrow=1)) else matrix(1,numCluster,numCluster)',
    '                p3<-if(i %in% index.antagonist)Matrix::crossprod(matrix(computeExpr_antagonist(avg,pairLRsig,cofactor_input,index.antagonist=i,Kh=Kh,n=n),nrow=1)) else matrix(1,numCluster,numCluster)',
    '                Pboot[,nE]<-as.vector(p1*p2*p3*P.spatial)',
    '              }',
    '            }')
  src<-c(src[seq_len(start-1)],block,src[seq(end+1,length(src))])
  e<-new.env(parent=asNamespace('CellChat'));e$fast_bootstrap_means<-fast_bootstrap_means
  fn<-eval(parse(text=src),envir=e)
  attr(fn,'acceleration_source')<-src
  fn
}
