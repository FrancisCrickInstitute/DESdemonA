## *** Useful Functions
## Only use certain samples
subsample <- function(dds, subs) {
  if (is_formula(subs)) {
    grps <- Reduce(interaction, colData(dds)[all.vars(subs)])
    grps <- names(table(grps))[table(grps)!=0]
    out <- lapply(setNames(grps,grps),  function(x) subsample(dds, grps==x))
  } else {
    out <- dds[,subs]
    colData(out) <- droplevels(colData(out))
    out
  }
  out
}


fit_models <- function(dds, ...) {
  lapply(metadata(dds)$model,
         function(mdl)  {
           out <- list()
           is_lrt <- sapply(mdl$comparison, is_formula)
           if (any(!is_lrt)) {
             wald <- dds
             design(wald) <- mdl$design
             wald <- DESeq(wald, test="Wald", ...)
             metadata(wald)$model <- mdl$design
             out <- lapply(mdl$comparison[!is_lrt], function(cntr) {
               metadata(wald)$comparison <- cntr
               wald})
           }
           if (any(is_lrt)) {
             lrt <- lapply(mdl$comparison[is_lrt],
                          function(reduced) {
                            dds_lrt <- fitLRT(dds, full=mdl$design, reduced=reduced, ...)
                            metadata(dds_lrt)$model <- mdl$design
                            metadata(dds_lrt)$comparison <- reduced
                            dds_lrt
                          })
             out <- c(out, lrt)
           }
           out
         }
         )
}

## Do likelihood ratio test, and classify in order of effect size
## Doesn't work for interactions, obviously

fitLRT <- function(dds, full, reduced, ...) {
  design(dds) <- full
  dds <- DESeq(dds, test="LRT", full=full, reduced=reduced, ...)
  ## cols <- resultsNames(dds)
  ## inBoth <- intersect(attr(terms(full), "term.labels"),
  ##                    attr(terms(reduced), "term.labels"))
  ## testTerm <- setdiff(attr(terms(full), "term.labels"),
  ##                    attr(terms(reduced), "term.labels"))
  ## if (length(inBoth)) {
  ##   cols <- cols[!grepl(inBoth, cols)]
  ## }
  ## toOrd <- mcols(dds)[cols]
  ## toOrd[,1] <- 0
  ## ords <-  apply(toOrd, 1, order)
  ## str <- matrix(levels(dds[[testTerm]])[ords], nrow=nrow(ords), ncol=ncol(ords))
  ## mcols(dds)$class <- apply(str, 2, paste, collapse="<")
  dds
}

## apply contrast, and transfer across interesting mcols from the dds

get_result <- function(dds, mcols=c("symbol", "entrez"), filterFun=IHW::ihw, ...) {
  comp <- metadata(dds)$comparison
  if (!is_formula(comp)) {
    if (is.character(comp) && length(comp)==1) { #  it's a name
      r <- results(dds, filterFun=filterFun, name=metadata(dds)$comparison, ...)
    } else { # it's a contrast
      if (is.list(comp) && "listValues" %in% names(comp)) {
        r <- results(dds, filterFun=filterFun, contrast=metadata(dds)$comparison[names(comp) != "listValues"], listValues=comp$listValues, ...)
      } else {
        r <- results(dds, filterFun=filterFun, contrast=metadata(dds)$comparison, ...)
      }
    }
  } else { # it's LRT
    r <- results(dds, filterFun=filterFun, ...)
  }
  r[mcols] <- mcols(dds)[mcols]
  if ("LRTPvalue" %in% names(mcols(dds))) {
    r$class <- mcols(dds)$class
    r$class[is.na(r$padj) | is.na(r$pvalue) | r$baseMean==0] <- NA
    term <- setdiff(.resNames(dds,  metadata(dds)$model),
                   .resNames(dds,  metadata(dds)$comparison))
    convertNames <- DESeq2:::renameModelMatrixColumns(colData(dds), metadata(dds)$model)
    convertNames <- subset(convertNames, from %in% term)
    term[match(convertNames$from, term)] <- convertNames$to
    # take the biggest fold-change vs baseline, for MA and reporting?
    effect_matrix <- as.matrix(mcols(dds)[,term])
    ind <- apply(abs(effect_matrix), 1, which.max)
    ind <- sapply(ind, function(x) ifelse(length(x)==0, NA, x))
    r$log2FoldChange <- effect_matrix[cbind(1:length(ind), ind)]
    r$lfcSE <- as.matrix(mcols(dds)[,paste0("SE_", term)])[cbind(1:length(ind), ind)]
    fit <- ashr::ash(r$log2FoldChange, r$lfcSE, mixcompdist = "normal", 
                    method = "shrink")
    r$shrunkLFC <- fit$result$PosteriorMean
    r$whichLFC <- term[ind]
    r$lfcSE <- fit$result$PosteriorSD
    effect_matrix <- cbind(0, effect_matrix)
    colnames(effect_matrix)[1] <- "0"
    ind <- apply(effect_matrix, 1, order)
    r$class <- apply(ind, 2, function(x) paste(colnames(effect_matrix)[x], collapse="<"))
  }  else {
    r$shrunkLFC <- lfcShrink(dds, res=r, type="ashr", quiet=TRUE)$log2FoldChange
    r$class <- ifelse(r$log2FoldChange >0, "Up", "Down")
  }
  ind <- which(r$padj<metadata(r)$alpha)
  r$class[ind] <- paste0(r$class[ind], "*")
  r$class[is.na(r$pvalue)] <- "Outlier"
  r$class[is.na(r$padj)] <- "Low Count"
  r$class[r$baseMean==0] <- "Zero Count"
  mcols(dds)$results <- r
  dds
}

.resNames <- function(dds, mdl) {
  make.names(
  colnames(stats::model.matrix.default(
    mdl, 
    data = as.data.frame(colData(dds))
  ))
  )
  }


summarise_results <- function(dds) {
  res <- mcols(dds)$results
  as.data.frame(table(Group=sub("\\*$","",res$class),
                      Significant=factor(ifelse(grepl("\\*$",res$class), "Significant", "not"), levels=c("Significant","not"))
                      )) %>%
    spread(Significant, Freq) %>%
    mutate(Total=not+Significant) %>%
    dplyr::select(-not) %>%
    arrange(desc(Significant/Total))
}    


# sigNot <- function(r) ifelse(grepl("\\*", r$class), "Sig","-")


qc_heatmap <- function(dds, pc_x=1, pc_y=2, batch=~1, title="QC Visualisation", header="\n\n##", n=500, caption=print) {
  cat(header, " ", title, "\n", sep="")
  qc_vis <- list()
  
  var_stab <- assay(vst(dds))
  if (batch != ~1) {
    var_stab <- residuals(limma::lmFit(var_stab, model.matrix(batch, as.data.frame(colData(dds)))), var_stab)
  }
  top <- order(apply(var_stab, 1, sd), decreasing=TRUE)[1:n]

  ### Heatmap
  cat(header, "# Heatmap of variable genes", "\n", sep="") 
  plotDat <- var_stab[top,]
  plotDat <- plotDat-rowMeans(plotDat)
  vars <- metadata(dds)$labels
  colDat <- as.data.frame(colData(dds))[vars]
  colnames(plotDat) <- rownames(colDat)
  qc_vis$heatmap <- Heatmap(plotDat, name="Mean Centred", column_title="Samples", row_title="Genes", 
                           col=colorspace::diverging_hcl(5, palette="Blue-Red"),
                           top_annotation=HeatmapAnnotation(df=colDat,
                                                            col = df2colorspace(colDat)),
                           show_row_names=FALSE, show_column_names=TRUE)
  print(qc_vis$heatmap)
  caption("Heatmap of variable genes")

  cat(header, "# Heatmap of sample distances", "\n", sep="") 
  poisd <- PoiClaClu::PoissonDistance(t(counts(dds, normalized=TRUE)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <-unite(colDat, col="exp_group", !!vars, remove=FALSE)[[1]]
  colnames(samplePoisDistMatrix) <- row.names(colDat)
  qc_vis$sample_dist <- Heatmap(
    samplePoisDistMatrix,
    clustering_distance_rows = function(x) {poisd$dd},
    clustering_distance_columns = function(x) {poisd$dd},
    top_annotation=HeatmapAnnotation(df=colDat,
                                     col = df2colorspace(colDat))
    )
  print(qc_vis$sample_dist)
  caption("Heatmap of sample distances")
  
  ### PCA
  pc <- prcomp(t(var_stab[top,]), scale=FALSE)
  percentVar <- round(100 * pc$sdev^2 / sum( pc$sdev^2 ))
  fml <- as.formula(paste0("~", paste(metadata(dds)$labels, collapse="+")))
  covvar_PC <- matrix(NA_real_, length(all.vars(fml)), ncol(pc$x),
                     dimnames=list(all.vars(fml), paste0("", 1:ncol(pc$x))))
  for (ipc in 1:ncol(pc$x)) {
    fit0 <- lm(pc$x[,ipc]  ~ 1, data=colData(dds))
    fit1 <- add1(fit0, fml, test="Chisq")[-1,"Pr(>Chi)", drop=FALSE]
    covvar_PC[rownames(fit1),ipc] <- as.vector(fit1[,1])
  }
  plotFrame <- expand.grid(Covariate=dimnames(covvar_PC)[[1]], PC=dimnames(covvar_PC)[[2]])
  plotFrame$p <- as.vector(covvar_PC)
  qc_vis$PC_phenotype <- ggplot(plotFrame, aes(x=PC, y=Covariate, fill=-log10(p))) +
    geom_raster() +
    scale_fill_gradient(low="grey90", high="red") +
    theme_classic() + coord_fixed()

  print(qc_vis$PC_phenotype)
  caption("Covariate-PC association")
  qc_vis$PC <- list()
  cat(header, "# Visualisation of PCs ", pc_x, " and ", pc_y, " coloured by covariate", "\n", sep="") 
  for (j in vars) {
    pc.df <- data.frame(PC1=pc$x[,pc_x], PC2=pc$x[,pc_y], col=colDat[[j]], sample=rownames(colDat))
    qc_vis$PC[[j]] <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=col))  + geom_point(size=3) +
      geom_text_repel(aes(label=sample)) +
      xlab(paste0("PC ", pc_x, ": ", percentVar[pc_x], "% variance")) +
      ylab(paste0("PC ", pc_y, ": ", percentVar[pc_y], "% variance")) +
      labs(colour=j) + coord_fixed()
    print(qc_vis$PC[[j]])
    caption(paste0("Coloured by ", j))
  }
  invisible(qc_vis)
}



enrichment <- function(ddsList, fun, showCategory, max_width=30) {
  genes <- lapply(ddsList, function(dds) {
    res <- mcols(dds)$results
    na.omit(res$entrez[grepl("\\*", res$class)])
  })
  genes <- genes[sapply(genes, length)!=0]
  if (length(genes)<=1) {
    return()
  }
  if (fun=="enrichGO") {
    reactome <- compareCluster(genes, fun=fun, OrgDb=metadata(ddsList[[1]])$organism$org, universe=na.omit(metadata(ddsList[[1]])$entrez))
  } else {
    orgs <- c(anopheles = "org.Ag.eg.db", arabidopsis = "org.At.tair.db", 
             bovine = "org.Bt.eg.db", canine = "org.Cf.eg.db", celegans = "org.Ce.eg.db", 
             chicken = "org.Gg.eg.db", chimp = "org.Pt.eg.db", coelicolor = "org.Sco.eg.db", 
             ecolik12 = "org.EcK12.eg.db", ecsakai = "org.EcSakai.eg.db", 
             fly = "org.Dm.eg.db", gondii = "org.Tgondii.eg.db", human = "org.Hs.eg.db", 
             malaria = "org.Pf.plasmo.db", mouse = "org.Mm.eg.db", 
             pig = "org.Ss.eg.db", rat = "org.Rn.eg.db", rhesus = "org.Mmu.eg.db", 
             xenopus = "org.Xl.eg.db", yeast = "org.Sc.sgd.db", zebrafish = "org.Dr.eg.db"
             )
    reactome_org <- names(orgs[orgs==metadata(ddsList[[1]])$organism$org])
    reactome <- compareCluster(genes, fun=fun, organism=reactome_org, universe=na.omit(metadata(ddsList[[1]])$entrez))
  }
  enrich_table <- as.data.frame(reactome)[c("Cluster", "ID", "Description","GeneRatio","BgRatio")]
  reactome@compareClusterResult$Description <- ifelse(
    nchar(reactome@compareClusterResult$Description)>max_width,
    reactome@compareClusterResult$ID,
    reactome@compareClusterResult$Description)
  pl <- dotplot(reactome, showCategory=showCategory) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6),
          axis.text.y = element_text(size=8)
          )
  list(plot=pl, table=enrich_table)
}

tidy_significant_dds <- function(dds, res, tidy_fn=NULL) {
  ind <- grepl("\\*$", res$class)
  mat <- assay(vst(dds))[ind,]
  tidy_dat <- tidy_per_gene(mat, as.data.frame(colData(dds)), tidy_fn)
  return(tidy_dat)
}

tidy_per_gene <- function(mat, pdat,  tidy_fn) {
  if (is.null(tidy_fn)) {
    return(list(mat=mat, pdat=pdat))
  }
  pdat$.value <- mat[1,]
  tidy_mat <- apply(mat, 1, function(x) {pdat$.value=x; tidy_fn(pdat)$.value})
  tidy_pdat <- tidy_fn(pdat) %>% dplyr::select(-.value)
  list(mat=t(tidy_mat), pdat=tidy_pdat)
}

## tidy_per_gene <- function(mat, pdat,  tidy_fn) {
##   long_mat <- tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(mat), var=".gene"), -.gene, names_to=".sample", values_to=".value")
##   long_df <- dplyr::left_join(long_mat, dplyr::mutate(pdat, .sample=long_mat$.sample[1:nrow(pdat)]), by=".sample")
##   long_tidied <- long_df %>% group_by(.gene) %>% tidy_fn()
##   tidied_mat <- tidyr::pivot_wider(ungroup(long_tidied), .gene, values_from=.value, names_from=.sample) %>%
##     tibble::column_to_rownames(".gene")
##   pdat_tidied <- tidy_fn(dplyr::mutate(pdat, .value=mat[1,])) %>% dplyr::select(-.value)
##   list(mat=tidied_mat[row.names(mat),], pdat=pdat_tidied)
## }



differential_heatmap <- function(ddsList, tidy_fn=NULL, title="Differential heatmap", header="\n\n##") {
  cat(header, " ", title, "\n", sep="")
  pal <- RColorBrewer::brewer.pal(12, "Set3")
  for (i in names(ddsList)) {
    tidied_data <- tidy_significant_dds(ddsList[[i]], mcols(ddsList[[i]])$results, tidy_fn)
    if (i == names(ddsList)[1]) {
      pdat <- tidied_data$pdat
      if (length(group_vars(pdat))) {
        column_split=pdat[group_vars(pdat)]
      } else {
        column_split=NULL
      }
      pdat <- mutate_if(as.data.frame(pdat), is.character, as.factor)
      colList <- df2colorspace(pdat)
    }
    name <- sub(".*\\t", "", i)
    ha <- HeatmapAnnotation(df=pdat, col=colList)
    pl <- Heatmap(tidied_data$mat,
                 name=sub(".*\\t", "", i),
                 cluster_columns = FALSE, show_column_names = FALSE,
                 column_split = column_split,
                 top_annotation = ha,
                 show_row_names = FALSE)
    draw(pl)
    caption(paste0("Heatmap of differential genes ", name))
  }
}



differential_MA <- function(ddsList, caption) {
  for (i in names(ddsList)) {
    res <- as.data.frame(mcols(ddsList[[i]])$results)
    md <- metadata(ddsList[[i]])
    pal <- RColorBrewer::brewer.pal(12, "Set3")
    if (is_formula(md$comparison)) {
      term <- sort(unique(mcols(ddsList[[i]])$results$whichLFC))
      res$class <- ifelse(grepl("\\*", res$class),res$whichLFC, "None")
      cols <- setNames(c(pal[1:length(term)], "grey"), c(term, "None"))
      yax <- "Largest fold change"
      lpos <- "bottom"
    } else {
      res$class <- ifelse(grepl("\\*", res$class), "Sig", "Not")
      cols <- setNames(c("blue","grey"), c("Sig","Not"))
      yax <- "Log fold change (shrunk)"
      lpos <- "none"
    }
    pl <- ggplot(res,
                aes(x=baseMean,
                    y=shrunkLFC,
                    colour=class,
                    size=ifelse(class %in% c("Not","None"),.1 , .5),
                    label=ifelse(class %in% c("Not","None"), "" , symbol)
                    )
                ) +
      geom_point() +
      scale_x_log10() + scale_colour_manual(values=cols) + scale_size_identity() +
      theme(legend.position=lpos) +
      labs(x="Mean of normalised counts", y=yax)
    if (sum(!(res$class %in% c("Not","None")))<30) {
      pl <- pl+ggrepel::geom_text_repel(size=4)
    }
    print(pl)
    caption(paste0("MA plot of differential genes ", sub(".*\\t", "", i)))
  }
}

write_results <- function(ddsList, param) {
  si <- session_info()
  crick_colours <-list(
    primary=list(red="#e3001a",yellow="#ffe608",blue="#4066aa",green="#7ab51d",purple="#bb90bd"),
    secondary=list(pink="#fadbd9", yellow="#fff6a7", green="#adcf82", orange="#ffe7ab", blue="#bee2e6"),
    spare=list(blue="#95ceed"))
  hs1 <- createStyle(fgFill = crick_colours$secondary$orange, textDecoration = "italic",
                    border = "Bottom")
  hs2 <- createStyle(fgFill = crick_colours$secondary$blue, textDecoration = "italic",
                    border = "Bottom")
  summaries <- map_depth(ddsList, 3, babsRNASeq::summarise_results)
  out <- lapply(ddsList, function(x) "")
  for (dataset in names(ddsList)) {
    wb <- openxlsx::createWorkbook(title=param$get("title"),
                                  creator="Gavin Kelly")
    tmp <- param$describe()
    dframe <- data.frame(id=names(tmp), description=unlist(tmp))
    sn <- "Parameters"
    addWorksheet(wb, sn)
    writeData(wb, sn, dframe,rowNames=FALSE, colNames=TRUE)
    ## Design
    samples_used <- as.data.frame(colData(ddsList[[dataset]][[1]][[1]]))
    sn <- "Design"
    addWorksheet(wb, sn)
    writeData(wb, sn, samples_used, headerStyle=hs2)
    sn <- "Class Sizes"
    addWorksheet(wb, sn)
    dframe <- babsRNASeq::rbind_summary(
      summaries[[dataset]],
      levels=c("Design","Comparison")
    )
    writeData(wb, sn, dframe, headerStyle=hs2)
    ## Differential gene-lists
    for (design_ind in 1:length(ddsList[[dataset]])) {
      for (contrast_name in names(ddsList[[dataset]][[design_ind]])) {
        dframe <- as.data.frame(mcols(ddsList[[dataset]][[design_ind]][[contrast_name]])$results) %>%
          tibble::rownames_to_column("id") %>%
          dplyr::filter(padj<param$get("alpha")) %>%
          dplyr::arrange(desc(shrunkLFC)) %>%
          dplyr::select(-pvalue, -padj)
        sn <- paste0(contrast_name, ", ", names(ddsList[[dataset]])[design_ind])
        addWorksheet(wb, sn, tabColour=crick_colours$secondary[[design_ind]])
        writeData(wb, sn, dframe, headerStyle=hs1)
      }
    }
    ## sn <- "GO terms"
    ## addWorksheet(wb, sn)
    ## writeData(wb, sn, go_df, headerStyle=hs1)
    sn <- "R Packages"
    addWorksheet(wb, sn)
    writeData(wb, sn, as.data.frame(si$packages), headerStyle=hs2)
    sn <- "R Details"
    addWorksheet(wb, sn)
    writeData(wb, sn, data.frame(setting = names(si$platform),
                                 value = unlist(si$platform),
                                 stringsAsFactors = FALSE),
              headerStyle=hs2)
    out[[dataset]] <- vDevice(file= paste0("differential_genelists_", dataset, ".xlsx"))
    saveWorkbook(wb, out[[dataset]], overwrite=TRUE)
  }
  out
}



write_all_results <- function(ddsList) {
  for (i in names(ddsList)) {
    for (j in names(ddsList[[i]])) { 
      for (k in names(ddsList[[i]][[k]])) { 
        vDevice(txt, file=paste(i,j,k, sep="-"), x=as.data.frame(mcols(ddsList[[i]][[j]][[k]])$results) %>% dplyr::select(log2FoldChange, stat, symbol, class))
      }
    }
  }
}
