function (object, features.name = "features", variable.both = TRUE, 
    features = NULL, return.object = TRUE) 
{
    gene.use <- row.names(object@data.signaling)
    DB <- object@DB
    if (is.null(features)) {
        if (is.list(object@var.features)) {
            features.sig <- object@var.features[[features.name]]
        }
        else {
            stop("Please update your CellChat object via `updateCellChat()`")
        }
    }
    else {
        features.sig <- features
    }
    interaction_input <- DB$interaction
    complex_input <- DB$complex
    my.sapply <- ifelse(test = future::nbrOfWorkers() == 1, yes = pbapply::pbsapply, 
        no = future.apply::future_sapply)
    complexSubunits <- complex_input[, grepl("subunit", colnames(complex_input))]
    index.sig <- unlist(x = my.sapply(X = 1:nrow(complexSubunits), 
        FUN = function(x) {
            complexsubunitsV <- unlist(complexSubunits[x, ], 
                use.names = F)
            complexsubunitsV <- complexsubunitsV[complexsubunitsV != 
                ""]
            if (length(intersect(complexsubunitsV, features.sig)) > 
                0 & all(complexsubunitsV %in% gene.use)) {
                return(x)
            }
        }))
    complexSubunits.sig <- complexSubunits[index.sig, ]
    index.use <- unlist(x = my.sapply(X = 1:nrow(complexSubunits), 
        FUN = function(x) {
            complexsubunitsV <- unlist(complexSubunits[x, ], 
                use.names = F)
            complexsubunitsV <- complexsubunitsV[complexsubunitsV != 
                ""]
            if (all(complexsubunitsV %in% gene.use)) {
                return(x)
            }
        }))
    complexSubunits.use <- complexSubunits[index.use, ]
    pairLR <- select(interaction_input, ligand, receptor)
    if (variable.both) {
        index.sig <- unlist(x = my.sapply(X = 1:nrow(pairLR), 
            FUN = function(x) {
                if (all(unlist(pairLR[x, ], use.names = F) %in% 
                  c(features.sig, rownames(complexSubunits.sig)))) {
                  return(x)
                }
            }))
    }
    else {
        index.sig <- unlist(x = my.sapply(X = 1:nrow(pairLR), 
            FUN = function(x) {
                if (all(unlist(pairLR[x, ], use.names = F) %in% 
                  c(gene.use, rownames(complexSubunits.use))) & 
                  (length(intersect(unlist(pairLR[x, ], use.names = F), 
                    c(features.sig, rownames(complexSubunits.sig)))) > 
                    0)) {
                  return(x)
                }
            }))
    }
    pairLRsig <- interaction_input[index.sig, ]
    object@LR$LRsig <- pairLRsig
    cat("The number of highly variable ligand-receptor pairs used for signaling inference is", 
        nrow(pairLRsig), "\n")
    if (return.object) {
        return(object)
    }
    else {
        return(pairLRsig)
    }
}
