function (object, data.use = NULL, group.by = NULL, idents.use = NULL, 
    invert = FALSE, group.dataset = NULL, pos.dataset = NULL, 
    group.DE.combined = FALSE, features.name = "features", only.pos = TRUE, 
    features = NULL, return.object = TRUE, thresh.pc = 0, thresh.fc = 0, 
    thresh.p = 0.05, do.DE = TRUE, do.fast = TRUE, min.cells = 10) 
{
    if (!is.list(object@var.features)) {
        stop("Please update your CellChat object via `updateCellChat()`")
    }
    if (is.null(data.use)) {
        X <- object@data.signaling
        if (nrow(X) < 3) {
            stop("Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks OK.")
        }
    }
    else {
        X <- data.use
    }
    if (is.null(features)) {
        features.use <- row.names(X)
    }
    else {
        features.use <- intersect(features, row.names(X))
    }
    data.use <- X[features.use, ]
    if (do.DE) {
        data.use <- as.matrix(data.use)
        if (is.null(group.by)) {
            labels <- object@idents
            if (!is.factor(labels)) {
                message("Use the joint cell labels from the merged CellChat object")
                labels <- object@idents$joint
            }
        }
        else {
            labels <- object@meta[[group.by]]
        }
        if (!is.factor(labels)) {
            labels <- factor(labels)
        }
        level.use <- levels(labels)[levels(labels) %in% unique(labels)]
        if (!is.null(idents.use)) {
            if (invert) {
                level.use <- level.use[!(level.use %in% idents.use)]
            }
            else {
                level.use <- level.use[level.use %in% idents.use]
            }
        }
        numCluster <- length(level.use)
        if (!is.null(group.dataset)) {
            labels.dataset <- as.character(object@meta[[group.dataset]])
            if (!(pos.dataset %in% unique(labels.dataset))) {
                cat("Please set pos.dataset to be one of the following dataset names: ", 
                  unique(as.character(labels.dataset)))
                stop()
            }
            labels.dataset[labels.dataset != pos.dataset] <- toString(setdiff(unique(labels.dataset), 
                pos.dataset))
            labels.dataset <- factor(labels.dataset, levels = c(pos.dataset, 
                setdiff(unique(labels.dataset), pos.dataset)))
        }
        if (do.fast) {
            presto.check <- rlang::is_installed(c("presto"))
            if (!presto.check) {
                stop("For a faster implementation of the Wilcoxon Test, please install the presto package", 
                  "\n--------------------------------------------", 
                  "\n devtools::install_github('immunogenomics/presto')", 
                  "\n--------------------------------------------", 
                  "\n Otherwise, plase set `do.fast = FALSE` for running the standard Wilcoxon Test!\n")
            }
            if (is.null(group.dataset)) {
                genes.de <- presto::wilcoxauc(data.use, labels, 
                  groups_use = level.use)
                colnames(genes.de) <- plyr::mapvalues(colnames(genes.de), 
                  from = c("group", "feature", "pval", "logFC", 
                    "pct_in", "pct_out", "padj"), to = c("clusters", 
                    "features", "pvalues", "logFC", "pct.1", 
                    "pct.2", "pvalues.adj"), warn_missing = TRUE)
                genes.de$logFC_abs <- abs(genes.de$logFC)
                pct.max <- apply(genes.de[, c("pct.1", "pct.2"), 
                  drop = FALSE], MARGIN = 1, FUN = max)
                genes.de$pct.max <- pct.max
                markers.all <- dplyr::filter(genes.de, pvalues < 
                  thresh.p, logFC_abs >= thresh.fc, pct.max > 
                  thresh.pc * 100) %>% arrange(pvalues)
            }
            else if ((!is.null(group.dataset)) & (group.DE.combined == 
                FALSE)) {
                genes.de <- data.frame()
                for (i in 1:numCluster) {
                  idx <- which(labels == level.use[i])
                  data.use.i <- data.use[, idx]
                  labels.i <- labels.dataset[idx]
                  genes.de.i <- presto::wilcoxauc(data.use.i, 
                    labels.i)
                  genes.de.i$clusters <- level.use[i]
                  genes.de <- rbind(genes.de, genes.de.i)
                }
                colnames(genes.de) <- plyr::mapvalues(colnames(genes.de), 
                  from = c("group", "feature", "pval", "logFC", 
                    "pct_in", "pct_out", "padj"), to = c("datasets", 
                    "features", "pvalues", "logFC", "pct.1", 
                    "pct.2", "pvalues.adj"), warn_missing = TRUE)
                genes.de$logFC_abs <- abs(genes.de$logFC)
                pct.max <- apply(genes.de[, c("pct.1", "pct.2"), 
                  drop = FALSE], MARGIN = 1, FUN = max)
                genes.de$pct.max <- pct.max
                markers.all <- dplyr::filter(genes.de, pvalues < 
                  thresh.p, logFC_abs >= thresh.fc, pct.max > 
                  thresh.pc * 100)
                markers.all$datasets <- factor(markers.all$datasets, 
                  levels = levels(labels.dataset))
                markers.all <- markers.all[order(markers.all$datasets, 
                  markers.all$pvalues, -markers.all$logFC), ]
            }
            else if ((!is.null(group.dataset)) & (group.DE.combined == 
                TRUE)) {
                genes.de.c <- presto::wilcoxauc(data.use, labels.dataset)
                genes.de.c <- genes.de.c[1:(nrow(genes.de.c)/2), 
                  ]
                genes.de <- data.frame()
                for (i in 1:numCluster) {
                  genes.de.c$clusters <- level.use[i]
                  genes.de <- rbind(genes.de, genes.de.c)
                }
                colnames(genes.de) <- plyr::mapvalues(colnames(genes.de), 
                  from = c("group", "feature", "pval", "logFC", 
                    "pct_in", "pct_out", "padj"), to = c("datasets", 
                    "features", "pvalues", "logFC", "pct.1", 
                    "pct.2", "pvalues.adj"), warn_missing = TRUE)
                genes.de$logFC_abs <- abs(genes.de$logFC)
                pct.max <- apply(genes.de[, c("pct.1", "pct.2"), 
                  drop = FALSE], MARGIN = 1, FUN = max)
                genes.de$pct.max <- pct.max
                markers.all <- dplyr::filter(genes.de, pvalues < 
                  thresh.p, logFC_abs >= thresh.fc, pct.max > 
                  thresh.pc * 100) %>% arrange(pvalues)
            }
            markers.all <- dplyr::select(markers.all, -c("logFC_abs", 
                "statistic", "pct.max"))
            if (only.pos & nrow(markers.all) > 0) {
                markers.all <- subset(markers.all, subset = logFC > 
                  0)
            }
            markers.all$features <- as.character(markers.all$features)
            features.sig <- markers.all$features
            object@var.features[[features.name]] <- unique(features.sig)
            features.name <- paste0(features.name, ".info")
            object@var.features[[features.name]] <- markers.all
        }
        else {
            my.sapply <- ifelse(test = future::nbrOfWorkers() == 
                1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
            mean.fxn <- function(x) {
                return(log(x = mean(x = expm1(x = x)) + 1))
            }
            labels <- as.character(labels)
            genes.de <- vector("list", length = numCluster)
            for (i in 1:numCluster) {
                features <- features.use
                if (is.null(group.dataset)) {
                  cell.use1 <- which(labels == level.use[i])
                  cell.use2 <- base::setdiff(1:length(labels), 
                    cell.use1)
                }
                else if ((!is.null(group.dataset)) & (group.DE.combined == 
                  FALSE)) {
                  cell.use1 <- which((labels == level.use[i]) & 
                    (labels.dataset == pos.dataset))
                  cell.use2 <- which((labels == level.use[i]) & 
                    (labels.dataset != pos.dataset))
                }
                else if ((!is.null(group.dataset)) & (group.DE.combined == 
                  TRUE)) {
                  cell.use1 <- which(labels.dataset == pos.dataset)
                  cell.use2 <- which(labels.dataset != pos.dataset)
                }
                thresh.min <- 0
                pct.1 <- round(x = rowSums(data.use[features, 
                  cell.use1, drop = FALSE] > thresh.min)/length(x = cell.use1), 
                  digits = 3)
                pct.2 <- round(x = rowSums(data.use[features, 
                  cell.use2, drop = FALSE] > thresh.min)/length(x = cell.use2), 
                  digits = 3)
                data.alpha <- cbind(pct.1, pct.2)
                colnames(x = data.alpha) <- c("pct.1", "pct.2")
                alpha.min <- apply(X = data.alpha, MARGIN = 1, 
                  FUN = max)
                names(x = alpha.min) <- rownames(x = data.alpha)
                features <- names(x = which(x = alpha.min > thresh.pc))
                if (length(x = features) == 0) {
                  next
                }
                data.1 <- apply(X = data.use[features, cell.use1, 
                  drop = FALSE], MARGIN = 1, FUN = mean.fxn)
                data.2 <- apply(X = data.use[features, cell.use2, 
                  drop = FALSE], MARGIN = 1, FUN = mean.fxn)
                FC <- (data.1 - data.2)
                if (only.pos) {
                  features.diff <- names(which(FC > thresh.fc))
                }
                else {
                  features.diff <- names(which(abs(FC) > thresh.fc))
                }
                features <- intersect(x = features, y = features.diff)
                if (length(x = features) == 0) {
                  next
                }
                data1 <- data.use[features, cell.use1, drop = FALSE]
                data2 <- data.use[features, cell.use2, drop = FALSE]
                pvalues <- unlist(x = my.sapply(X = 1:nrow(x = data1), 
                  FUN = function(x) {
                    return(wilcox.test(data1[x, ], data2[x, ])$p.value)
                  }))
                pval.adj = stats::p.adjust(p = pvalues, method = "bonferroni", 
                  n = nrow(X))
                genes.de[[i]] <- data.frame(clusters = level.use[i], 
                  features = as.character(rownames(data1)), pvalues = pvalues, 
                  logFC = FC[features], data.alpha[features, 
                    , drop = F], pvalues.adj = pval.adj, stringsAsFactors = FALSE)
            }
            markers.all <- data.frame()
            for (i in 1:numCluster) {
                gde <- genes.de[[i]]
                if (!is.null(gde)) {
                  gde <- gde[order(gde$pvalues, -gde$logFC), 
                    ]
                  gde <- subset(gde, subset = pvalues < thresh.p)
                  if (nrow(gde) > 0) {
                    markers.all <- rbind(markers.all, gde)
                  }
                }
            }
            if (only.pos & nrow(markers.all) > 0) {
                markers.all <- subset(markers.all, subset = logFC > 
                  0)
            }
            if (!is.null(group.dataset)) {
                markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
                markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), 
                  pos.dataset)
                markers.all$datasets <- factor(markers.all$datasets, 
                  levels = levels(labels.dataset))
                markers.all <- markers.all[order(markers.all$datasets, 
                  markers.all$pvalues, -markers.all$logFC), ]
            }
            markers.all$features <- as.character(markers.all$features)
            features.sig <- markers.all$features
            object@var.features[[features.name]] <- features.sig
            features.name <- paste0(features.name, ".info")
            object@var.features[[features.name]] <- markers.all
        }
    }
    else {
        markers.all <- data.frame(features = as.character(rownames(data.use)), 
            nCells = rowSums(data.use > 0))
        markers.all <- dplyr::filter(markers.all, nCells >= min.cells)
        features.sig <- markers.all$features
        object@var.features[[features.name]] <- unique(features.sig)
        features.name <- paste0(features.name, ".info")
        object@var.features[[features.name]] <- markers.all
    }
    if (return.object) {
        return(object)
    }
    else {
        return(markers.all)
    }
}
