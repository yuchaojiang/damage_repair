plotCounts.RABMP=function (ddsRA, ddsBMP, gene, intgroup = "condition", normalized = TRUE, 
                         transform = TRUE, main, xlab = "group", returnData = FALSE, 
                         replaced = FALSE, pc, ...) 
{
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(ddsRA)))))
  if (!all(intgroup %in% names(colData(ddsRA)))) 
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(ddsRA)[[v]], 
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (transform) 
      0.5
    else 0
  }
  if (is.null(sizeFactors(ddsRA)) & is.null(normalizationFactors(ddsRA))) {
    ddsRA <- estimateSizeFactors(ddsRA)
  }
  cnts <- counts(ddsRA, normalized = normalized, replaced = replaced)[gene, 
                                                                    ]
  group <- if (length(intgroup) == 1) {
    colData(ddsRA)[[intgroup]]
  } else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(ddsRA)[[intgroup[1]]]), 
                              levels(colData(ddsRA)[[intgroup[2]]]), function(x, 
                                                                            y) paste(x, y, sep = ":"))))
    droplevels(factor(apply(as.data.frame(colData(ddsRA)[, 
                                                       intgroup, drop = FALSE]), 1, paste, collapse = ":"), 
                      levels = lvls))
  } else {
    factor(apply(as.data.frame(colData(ddsRA)[, intgroup, drop = FALSE]), 
                 1, paste, collapse = ":"))
  }
  dataRA <- data.frame(count = cnts + pc, group = as.integer(group))
  
  
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                          (gene >= 1 & gene <= nrow(ddsBMP)))))
  if (!all(intgroup %in% names(colData(ddsBMP)))) 
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(ddsBMP)[[v]], 
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (tBMPnsform) 
      0.5
    else 0
  }
  if (is.null(sizeFactors(ddsBMP)) & is.null(normalizationFactors(ddsBMP))) {
    ddsBMP <- estimateSizeFactors(ddsBMP)
  }
  cnts <- counts(ddsBMP, normalized = normalized, replaced = replaced)[gene, 
                                                                       ]
  group <- if (length(intgroup) == 1) {
    colData(ddsBMP)[[intgroup]]
  } else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(ddsBMP)[[intgroup[1]]]), 
                              levels(colData(ddsBMP)[[intgroup[2]]]), function(x, 
                                                                               y) paste(x, y, sep = ":"))))
    droplevels(factor(apply(as.data.fBMPme(colData(ddsBMP)[, 
                                                           intgroup, drop = FALSE]), 1, paste, collapse = ":"), 
                      levels = lvls))
  } else {
    factor(apply(as.data.fBMPme(colData(ddsBMP)[, intgroup, drop = FALSE]), 
                 1, paste, collapse = ":"))
  }
  dataBMP <- data.frame(count = cnts + pc, group = as.integer(group))
  
  
  
  
  
  logxy <- if (transform) 
    "y" else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(ddsRA)[gene]
    } else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    return(data.frame(count = data$count, colData(ddsRA)[intgroup]))
  plot(dataRA$group + runif(ncol(ddsRA), -0.05, 0.05), dataRA$count, 
       xlim = c(0.5, max(dataRA$group) + 0.5),
       ylim = c(min(c(dataRA$count, dataBMP$count)),max(c(dataRA$count, dataBMP$count))),
       log = logxy, xaxt = "n", 
       xlab = xlab, ylab = ylab, main = main, col="#E69F00",...)
  data.avg=aggregate(dataRA$count, list(dataRA$group), mean)
  points(data.avg$Group.1, data.avg$x,type='l',col="#E69F00")
  
  points(dataBMP$group + runif(ncol(ddsBMP), -0.05, 0.05), dataBMP$count, 
       xlim = c(0.5, max(dataBMP$group) + 0.5),
       ylim = c(min(c(dataRA$count, dataBMP$count)),max(c(dataRA$count, dataBMP$count))),
       log = logxy, xaxt = "n", 
       xlab = xlab, ylab = ylab, main = main,col="#009E73",...)
  data.avg=aggregate(dataBMP$count, list(dataBMP$group), mean)
  points(data.avg$Group.1, data.avg$x,type='l',col="#009E73")
  axis(1, at = seq_along(levels(group)), levels(group))
  
}