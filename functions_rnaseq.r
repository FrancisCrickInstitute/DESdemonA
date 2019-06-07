applySubset <- function(dds, subs) {
  tmp <- dds[,subs]
  colData(tmp) <- droplevels(colData(tmp))
  tmp
}
applyLRT <- function(dds, mdl) {
  design(dds) <- mdl$full
  dds <- DESeq(dds, test="LRT", full=mdl$full, reduced=mdl$reduced)
  cols <- resultsNames(dds)
  inBoth <- intersect(attr(terms(mdl$full), "term.labels"),
                     attr(terms(mdl$reduced), "term.labels"))
  testTerm <- setdiff(attr(terms(mdl$full), "term.labels"),
                     attr(terms(mdl$reduced), "term.labels"))
  if (length(inBoth)) {
    cols <- cols[!grepl(inBoth, cols)]
  }
  toOrd <- mcols(dds)[cols]
  toOrd[,1] <- 0
  ords <-  apply(toOrd, 1, order)
  str <- matrix(levels(dds[[testTerm]])[ords], nrow=nrow(ords), ncol=ncol(ords))
  mcols(dds)$class <- apply(str, 2, paste, collapse="<")
  dds
}
applyResults <- function(dds, ...) {
  r <- results(dds, ...)
  r$symbol <- mcols(dds)$symbol
  r$tcr <- mcols(dds)$tcr
  r$treg <- mcols(dds)$treg
  if ("LRTPvalue" %in% names(mcols(dds))) {
    r$class <- mcols(dds)$class
    r$class[is.na(r$padj) | is.na(r$pvalue) | r$baseMean==0] <- NA
    r$log2FoldChange <- NULL
    r$lfcSE <- NULL
  }  else {
    r$class <- ifelse(r$log2FoldChange >0, "Up", "Down")
  }
  ind <- which(r$padj<metadata(r)$alpha)
  r$class[ind] <- paste0(r$class[ind], "*")
  r$class[is.na(r$pvalue)] <- "Outlier"
  r$class[is.na(r$padj)] <- "Low Count"
  r$class[r$baseMean==0] <- "Zero Count"
  r
}
sigNot <- function(r) ifelse(grepl("\\*", r$class), "Sig","-")
