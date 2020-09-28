## *** Useful Functions


##' Adjust a nested batch effect
##'
##' When one covariate entirely predicts another, then including both in a
##' linear model will result in an unspecified model.  This can happen if we
##' have e.g. groups of cell-lines that represent a genotype.  If each of those
##' cell-lines receives a treatment, then we need both the group-level information to
##' test the interesting treatment x genotype interaction, but also the cell-line information
##' to remove any batch effect they represent.  To achieve this, we can recode the batch effect
##' so that each genotype has a common set of levels, and then we can include that batch effect
##' in interaction with the grouping factor.
##'
##' If one group has a larger number of sub-treatments, then it will be necessary to remove
##' the columns of the design matrix that are all zero separately
##' @title Adjust a confounded batch effect
##' @param inner A factor representing the variable to be recoded, e.g. the cell-line
##' @param within The parent factor that 'inner' is inside, e.g. the genotype of the cell-line
##' @return A recoded factor that can now replace the inner (cell-line) variable in your model
##' @author Gavin Kelly
nest_batch <- function(inner, within, set_to=".") {
  n_instances <- apply(table(inner, within)!=0, 1, sum)
  if (any(n_instances>1)) {
    stop(paste(names(n_instances)[n_instances>1], collapse=", "), " appear in multiple parents")
  }
  levels_in <- levels(inner)
  # each inner level has in a unique group, find it
  corresponding_group <-within[match(levels_in, inner)]
  # for each group, find index of first inner level
  which_inner <- match(levels(within), corresponding_group)
  # and set it to the common value
  if (any(levels(inner)[-which_inner]==set_to )) { # we're going to duplicate an existing level so best stop
    stop("One of the batches is already called ", set_to, ". Please use a non-existing level")
  }
  levels(inner)[which_inner] <- set_to
  inner
  }

recode_within <- function(inner, within) {
  tab <- table(inner, within)!=0 # which batches are in which nest
  if (any(rowSums(tab)>1)) {
    stop("Some inner levels appear in multiple outer groups.")
  }
  factor(apply(tab, 2, cumsum)[cbind(as.character(inner),as.character(within))]) # cumsum to get incrementing index within group.
  }

clean_nested_imbalance <- function(dds) {
  mm <- model.matrix(design(dds), colData(dds))
  mm <- mm[,!apply(mm==0, 2, all)] # remove zeroed columns
  design(dds) <- mm
  dds
}


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

mult_comp <- function(spec, ...) {
  obj <- list(spec=spec,...)
  class(obj) <- "post_hoc"
  obj
}

emcontrasts <- function(dds, spec, ...) {
  df <- as.data.frame(colData(dds))
  df$.x <- counts(dds, norm=TRUE)[1,]
  if (is_formula(design(dds))) {
    fml <- as.formula(paste0(".x ~ ", as.character(design(dds)[2])))
    fit <- lm(fml, data=df)
    dds_ind <- seq_along(resultsNames(dds))
  } else {
    mm <- design(dds)
    ind <- colnames(mm) == "(Intercept)"
    colnames(mm)[ind] <- "Intercept"
    colnames(mm) <- make.names(colnames(mm))
    fit <- lm(df$.x ~ . -1, data.frame(mm))
    ddsNames <- match(resultsNames(dds), names(coef(fit)))
  }
  emfit <- emmeans::emmeans(fit, spec,...)
  contr_frame <- as.data.frame(summary(emfit$contrasts))
  contr_frame <- contr_frame[1:(which(names(contr_frame)=="estimate")-1)]
  contr <- lapply(seq_len(nrow(contr_frame)), function(i) emfit$contrast@linfct[i,])
  names(contr) <- do.call(paste, c(contr_frame,sep= "|"))
  contr
}


fit_models <- function(dds, ...) {
  lapply(metadata(dds)$model,
         function(mdl)  {
           out <- list()
           is_lrt <- sapply(mdl$comparison, is_formula)
           if (any(!is_lrt)) {
             wald <- dds
             design(wald) <- mdl$design
             comps <- mdl$comparison[!is_lrt]
             is_post_hoc <- sapply(comps, class)=="post_hoc"
             if (any(is_post_hoc)) {
               comps[is_post_hoc] <- lapply(comps[is_post_hoc], function(ph) {
                 do.call(emcontrasts, list(dds=wald, spec=ph$spec, ph[-1]))
               })
               comps[!is_post_hoc] <- lapply(comps[!is_post_hoc], list) # protect existing lists from unlist
               comps <- unlist(comps, recursive=FALSE)
             }
             wald <- DESeq(wald, test="Wald", ...)
             metadata(wald)$model <- mdl$design
             out <- lapply(comps, function(cntr) {
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
  if (is.null(filterFun)) filterFun <- rlang::missing_arg()
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
    colnames(effect_matrix) <- sub(paste0("(", paste0( all.vars(metadata(dds)$model), "_", collapse="|"),")"),  "", colnames(effect_matrix))
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
                           heatmap_legend_param = list(direction = "horizontal" ),
                           col=colorspace::diverging_hcl(5, palette="Blue-Red"),
                           top_annotation=HeatmapAnnotation(df=colDat,
                                                            col = df2colorspace(colDat)),
                           show_row_names=FALSE, show_column_names=TRUE)
  draw(qc_vis$heatmap, heatmap_legend_side="top")
  caption("Heatmap of variable genes")

  cat(header, "# Heatmap of sample distances", "\n", sep="") 
  poisd <- PoiClaClu::PoissonDistance(t(counts(dds, normalized=TRUE)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <-unite(colDat, col="exp_group", !!vars, remove=FALSE)[[1]]
  colnames(samplePoisDistMatrix) <- row.names(colDat)
  qc_vis$sample_dist <- Heatmap(
    samplePoisDistMatrix,
    name="Poisson Distance", 
    clustering_distance_rows = function(x) {poisd$dd},
    clustering_distance_columns = function(x) {poisd$dd},
    heatmap_legend_param = list(direction = "horizontal" ),
    top_annotation=HeatmapAnnotation(df=colDat,
                                     col = df2colorspace(colDat))
    )
  draw(qc_vis$sample_dist, heatmap_legend_side="top")
  caption("Heatmap of sample distances")
  
  ### PCA
  pc <- prcomp(t(var_stab), scale=FALSE)
  percentVar <- round(100 * pc$sdev^2 / sum( pc$sdev^2 ))
  fml <- as.formula(paste0("~", paste(metadata(dds)$labels, collapse="+")))
  covvar_PC <- matrix(NA_real_, length(all.vars(fml)), ncol(pc$x),
                     dimnames=list(all.vars(fml), paste0("", 1:ncol(pc$x))))
  for (ipc in 1:ncol(pc$x)) {
    fit0 <- lm(pc$x[,ipc]  ~ 1, data=colData(dds))
    fit1 <- add1(fit0, fml, test="Chisq")
#    covvar_PC[rownames(fit1)[-1],ipc] <- -log10(fit1$`Pr(>Chi)`[-1])
    covvar_PC[rownames(fit1)[-1],ipc] <- 1-fit1$RSS[2]/fit1$RSS[1]
  }
  plotFrame <- expand.grid(Covariate=dimnames(covvar_PC)[[1]], PC=dimnames(covvar_PC)[[2]])
  plotFrame$Assoc <- as.vector(covvar_PC)
  qc_vis$PC_phenotype <- ggplot(plotFrame, aes(x=PC, y=Covariate, fill=Assoc)) +
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
  qc_vis$PC_assoc <- list()
  cat(header, "# Visualisation of associated PCs coloured by covariate", "\n", sep="") 
  for (j in vars) {
    tmp <- order(covvar_PC[j,], decreasing=TRUE)
    pc_x <- tmp[1]
    pc_y <- tmp[2]
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



differential_heatmap <- function(ddsList, tidy_fn=NULL, caption) {
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
      pdat[sapply(pdat, is.character)] <- lapply(pdat[sapply(pdat, is.character)], 
                                                as.factor)
      colList <- df2colorspace(pdat)
    }
    colnames(tidied_data$mat) <- rownames(pdat)
    name <- sub(".*\\t", "", i)
    ha <- HeatmapAnnotation(df=pdat, col=colList)
    pl <- Heatmap(tidied_data$mat,
                 heatmap_legend_param = list(direction = "horizontal" ),
                 name=sub(".*\\t", "", i),
                 cluster_columns = is.null(column_split), show_column_names = TRUE,
                 column_split = column_split,
                 top_annotation = ha,
                 row_names_gp = gpar(fontsize = 6),
                 show_row_names = nrow(tidied_data$mat)<100)
    draw(pl, heatmap_legend_side="top")
    caption(paste0("Heatmap on differential genes ", name))
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
    ind <- order(res$shrunkLFC)[c(1:5, nrow(res)-(0:4))]
    pl <- pl+ggrepel::geom_text_repel(size=4, data=res[ind,])
    print(pl)
    caption(paste0("MA plot of differential genes ", sub(".*\\t", "", i)))
  }
}

write_results <- function(ddsList, param, dir=".") {
  si <- session_info()
  crick_colours <-list(
    primary=list(red="#e3001a",yellow="#ffe608",blue="#4066aa",green="#7ab51d",purple="#bb90bd"),
    secondary=list(pink="#fadbd9", yellow="#fff6a7", green="#adcf82", orange="#ffe7ab", blue="#bee2e6"),
    spare=list(blue="#95ceed"))
  hs1 <- createStyle(fgFill = crick_colours$secondary$orange, textDecoration = "italic",
                    border = "Bottom")
  hs2 <- createStyle(fgFill = crick_colours$secondary$blue, textDecoration = "italic",
                    border = "Bottom")
  summaries <- map_depth(ddsList, 3, babsrnaseq::summarise_results)
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
    dframe <- babsrnaseq::rbind_summary(
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
        if (length(ddsList[[dataset]])==1) {
          sn <- contrast_name
        } else {
          sn <- paste0(contrast_name, ", ", names(ddsList[[dataset]])[design_ind])
        }
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
    out[[dataset]] <- file.path(dir, paste0("differential_genelists_", dataset, ".xlsx"))
    saveWorkbook(wb, out[[dataset]], overwrite=TRUE)
  }
  out
}



write_all_results <- function(ddsList, dir=".") {
  for (i in names(ddsList)) {
    for (j in names(ddsList[[i]])) { 
      for (k in names(ddsList[[i]][[k]])) { 
        readr::write_excel_csv(
          path=file.path(dir, sprintf("allgenes_%s_%s_%s.csv", i, j, k)),
          x=as.data.frame(mcols(ddsList[[i]][[j]][[k]])$results) %>% dplyr::select(log2FoldChange, stat, symbol, class))
      }
    }
  }
}
