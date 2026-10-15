args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
suppressPackageStartupMessages(library(CellChat))
base <- file.path(args[1],'Cancers_ATRT_TME_Submission/Revision_20261010_CellChat10000/analysis')
rows <- list()
for (folder in list.dirs(file.path(base,'cellchat'),recursive=FALSE)) {
  stopifnot(file.exists(file.path(folder,'COMPLETE.json')))
  x <- readRDS(file.path(folder,'CellChat.rds'))
  input <- read.csv(file.path(folder,'input_cells.csv'),stringsAsFactors=FALSE)
  exported <- read.csv(file.path(folder,'all_tested_interactions.csv'),stringsAsFactors=FALSE)
  stopifnot(identical(colnames(x@data.signaling),input$Cell),
            identical(rownames(x@meta),input$Cell),
            identical(as.character(x@idents),input$group),
            identical(as.character(x@meta$ID),input$ID),
            x@options$parameter$nboot==10000,x@options$parameter$seed.use==42)
  z <- as.data.frame.table(x@net$prob,responseName='probability',stringsAsFactors=FALSE)
  names(z)[1:3] <- c('source','target','interaction_name')
  z$p_value <- as.vector(x@net$pval)
  z$p_plus1 <- (round(z$p_value*10000)+1)/10001
  z$q <- p.adjust(z$p_plus1,'BH')
  m <- merge(z,exported,by=c('source','target','interaction_name'),suffixes=c('_rds','_csv'))
  stopifnot(nrow(m)==nrow(z),max(abs(m$probability_rds-m$probability_csv))<1e-12,
            max(abs(m$p_value_rds-m$p_value_csv))<1e-12,
            max(abs(m$q-m$q_BH_plus1_within_object))<1e-12,
            all(is.finite(m$probability_rds)),all(m$p_plus1>0 & m$p_plus1<=1))
  rows[[length(rows)+1]] <- data.frame(run=basename(folder),n_cells=nrow(input),
    n_tests=nrow(z),n_permutations=10000,serialized_metadata_verified=TRUE,
    serialized_probabilities_verified=TRUE,serialized_p_values_verified=TRUE,
    full_family_BH_verified=TRUE)
}
stopifnot(length(rows)==6)
write.csv(do.call(rbind,rows),file.path(base,'tables/S13_saved_object_validation.csv'),row.names=FALSE)
cat('All six serialized objects and full statistical exports verified.\n')
