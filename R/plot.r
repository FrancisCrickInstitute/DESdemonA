##' Generate visualisations of raw data
##'
##' Plot heatmaps of the raw data, along with plots of the samples
##' projected onto various slices of the PCA space.
##' @title Data Visualisation
##' @param dds The [DESeq2::DESeqDataSet-class] object
##' @param pc_x The default choice for which principal component to
##'   plot on the x axis
##' @param pc_y The default choice for which principal component to
##'   plot on the y axis
##' @param batch A null model to take the residuals against, if you
##'   want to visualise the plots having removed any nuisance factors.
##' @param family Choice of PCA method
##' @param title Data Visualisation
##' @param header The markdown prefix that should precede any section
##' @param n The number of most-variable-genes to use for clustering
##' @param caption The function to use to create captions
##' @return A list (invisible) of ggplot2 objects
##' @author Gavin Kelly
#' @export
qc_heatmap <- function(dds, pc_x=1, pc_y=2, batch=~1, family="norm", title="QC Visualisation", header="\n\n##", n=500, caption=print) {
  cat(header, " ", title, "\n", sep="")
  qc_vis <- list()
  
  var_stab <- assay(dds, "vst")
  top <- order(apply(var_stab, 1, sd), decreasing=TRUE)[1:n] 

  ### Heatmap
  cat(header, "# Heatmap of variable genes", "\n", sep="") 
  plotDat <- var_stab[top,]
  plotDat <- plotDat-rowMeans(plotDat)
  vars <- get_terms(dds)
  colDat <- as.data.frame(colData(dds))
  colnames(plotDat) <- rownames(colDat)
  qc_vis$heatmap <- ComplexHeatmap::Heatmap(plotDat, name="Mean Centred", column_title="Samples", row_title="Genes",
                           heatmap_legend_param = list(direction = "horizontal" ),
                           col=colorspace::diverging_hcl(5, palette="Blue-Red"),
                           top_annotation=ComplexHeatmap::HeatmapAnnotation(df=colDat[vars],
                                                            col = df2colorspace(colDat[vars])),
                           show_row_names=FALSE, show_column_names=TRUE)
  draw(qc_vis$heatmap, heatmap_legend_side="top")
  caption("Heatmap of variable genes")

  cat(header, "# Heatmap of sample distances", "\n", sep="") 
  poisd <- PoiClaClu::PoissonDistance(t(counts(dds, normalized=TRUE)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <-tidyr::unite(colDat[vars], col="exp_group", !!vars, remove=FALSE)[[1]]
  colnames(samplePoisDistMatrix) <- row.names(colDat)
  qc_vis$sample_dist <- ComplexHeatmap::Heatmap(
    samplePoisDistMatrix,
    name="Poisson Distance", 
    clustering_distance_rows = function(x) {poisd$dd},
    clustering_distance_columns = function(x) {poisd$dd},
    heatmap_legend_param = list(direction = "horizontal" ),
    top_annotation=ComplexHeatmap::HeatmapAnnotation(df=colDat[vars],
                                     col = df2colorspace(colDat[vars]))
    )
  draw(qc_vis$sample_dist, heatmap_legend_side="top")
  caption("Heatmap of sample distances")
  
  ### PCA
  pc <- as.matrix(colData(dds)$.PCA)
  percentVar <- metadata(colData(dds)$.PCA)$percentVar
  is_vary <- sapply(colData(dds)[vars], function(v) length(unique(v))!=1)
  #fml <- as.formula(paste0("~", paste(metadata(dds)$labels[is_vary], collapse="+")))

  qc_vis$PC <- list(list())
  cat(header, "# Visualisation of PCs ", pc_x, " and ", pc_y, " coloured by covariate", "\n", sep="")
  do_labels <- nrow(colDat)<10
  for (j in vars[is_vary]) {
    pc.df <- data.frame(PC1=pc[,pc_x], PC2=pc[,pc_y], col=colDat[[j]], sample=rownames(colDat))
    qc_vis$PC[[1]][[j]] <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=col))  + geom_point(size=3) +
      xlab(paste0("PC ", pc_x, ": ", percentVar[pc_x], "% variance")) +
      ylab(paste0("PC ", pc_y, ": ", percentVar[pc_y], "% variance")) +
      labs(colour=j)
    if (do_labels) {qc_vis$PC[[1]][[j]] <- qc_vis$PC[[1]][[j]] + geom_text_repel(aes(label=sample))}
    print(qc_vis$PC[[1]][[j]])
    caption(paste0("Coloured by ", j))
  }
  
  qc_vis$PC_phenotype <- list()
  fitFrame <- colDat
  yvar <- make.unique(c(colnames(fitFrame), "y", sep=""))[ncol(fitFrame)+1]
  models_for_qc <- sapply(metadata(dds)$models, "[[", "plot_qc")
  for (model_name in names(metadata(dds)$models)[models_for_qc]) {
    qc_vis$PC[[model_name]] <- list()
    fml <- update(metadata(dds)$models[[model_name]]$design, paste(yvar, "~ ."))
    plotFrame <- expand.grid(Covariate=attr(terms(fml), "term.labels"),
                            PC=1:ncol(pc))
    plotFrame$Assoc <- NA
    npc <- ncol(pc)
    for (ipc in 1:npc) {
      fitFrame[[yvar]] <- pc[,ipc]
      fit1 <-  lm(fml, data=fitFrame)
      fit_selected <- MASS::stepAIC(fit1, trace=0)
      ind <- attr(terms(fit_selected),"term.labels")
      if (length(ind)) {
        fit0 <- as.data.frame(anova(fit_selected))
        rss <- fit0[ind,"Sum Sq"]/sum(fit0[,"Sum Sq"])
        ind_pc <- plotFrame$PC==ipc
        ind_fit <- match(ind, plotFrame$Covariate[ind_pc])
        plotFrame$Assoc[ind_pc][ind_fit] <- rss
      }
    }
    plotFrame$wrap <- (plotFrame$PC-1) %/% 20
    plotFrame$wrap <- paste0("PCs ", plotFrame$wrap*20+1, "-", min((plotFrame$wrap+1) * 20, npc))
    plotFrame$PC <- sprintf("%02d", plotFrame$PC)
    qc_vis$PC_phenotype[[model_name]] <- ggplot(plotFrame, aes(x=PC, y=Covariate, fill=Assoc)) +
      geom_raster() +
      facet_wrap(~wrap, scales="free_x", ncol=1) + 
      scale_fill_gradient2(low="#4575b4", mid="grey90", high="#d73027") +
      theme_classic() + theme(aspect.ratio = length(unique(plotFrame$Covariate)) / min(20, npc))
    print(qc_vis$PC_phenotype[[model_name]])
    caption("Covariate-PC association")
    
    cat(header, "# Visualisation of associated PCs coloured by covariate in ",model_name, "\n", sep="")
    qc_vis$PC[[model_name]] <- list()
    factors_1 <- attr(terms(fit1),"factors")
    for (j in unique(plotFrame$Covariate)) {
      this_covar <- subset(plotFrame, Covariate==j & !is.na(Assoc) & PC!=PC[nrow(plotFrame)])
      if (nrow(this_covar)==0) next
      this_covar <- this_covar[order(this_covar$Assoc, decreasing=TRUE),]
      pc_x <- as.integer(as.character(this_covar$PC[1]))
      if (nrow(this_covar)==1) {
        pc_y <- pc_x
      } else {
        pc_y <- as.integer(as.character(this_covar$PC[2]))
      }
      inter_vars <- rownames(factors_1)[factors_1[,j]!=0]
      pc.df <- data.frame(X=pc[,pc_x], Y=pc[,pc_y], sample=rownames(colDat))
      lab1 <- paste0("PC ", pc_x, ": ", percentVar[pc_x], "% variance")
      lab2 <- paste0("PC ", pc_y, ": ", percentVar[pc_y], "% variance")
      if (length(inter_vars)>1) {
        if (length(inter_vars)==2) {
          pc.df$Y <- pc.df$X # ie PC1
          pc.df$col <- colDat[,inter_vars[1]]
          pc.df$X <- colDat[,inter_vars[2]]
          my_lbl <- labs(colour=inter_vars[1], x=inter_vars[2], y=lab1)
        } else {
          pc.df$Y <- pc.df$X
          pc.df$X <- colDat[,inter_vars[1]]
          pc.df$col <- Reduce(interaction, colDat[,inter_vars[-1]])
          my_lbl <- labs(
            x=inter_vars[1],
            colour=paste(inter_vars[-1],collapse="x"),
            y=lab1)
        }
        my_aes <- ggplot2::aes(x=X, y=Y, colour=col)
      } else {
        pc.df$col <- colDat[[j]]
        my_aes <- ggplot2::aes(x=X, y=Y, colour=col)
        my_lbl <- labs(colour=j,x=lab1, y=lab2)
      }
      qc_vis$PC[[model_name]][[j]] <- ggplot(pc.df, my_aes)  + geom_point(size=3) +
        my_lbl
      if (do_labels) {qc_vis$PC[[model_name]][[j]] <- qc_vis$PC[[model_name]][[j]] + geom_text_repel(aes(label=sample))}
      print(qc_vis$PC[[model_name]][[j]])
      caption(paste0("Coloured by ", j))
    }
  }
  invisible(qc_vis)
}

get_terms <- function(dds) {
  if ("model" %in% metadata(dds)) {
    return(all.vars(metadata(dds)$model))
  } else {
    term_list <- lapply(metadata(dds)$models, function(mdl) {all.vars(mdl$design)})
    return(unique(unlist(term_list)))
  }
}
    

##' Heatmaps of expression for differential genes
##'
##' For each set of results, plot the heatmap of counts, limited to
##' differential genes
##' @title Heatmaps of differential genes
##' @param ddsList DESdemonA-generated list of [DESeq2::DESeqDataSet-class]s
##' @param tidy_fn A dplyr pipeline to transform the expression values
##'   in a convenient manner. Use this for example to 'normalise' the
##'   data so that all values are relative to a particular condition.
##' @param caption The function to use to create captions
##' @return
##' @author Gavin Kelly
#' @export
differential_heatmap <- function(ddsList, tidy_fn=NULL, caption) {
  pal <- RColorBrewer::brewer.pal(12, "Set3")
  first_done <- FALSE
  for (i in names(ddsList)) {
    if (!any(grepl("\\*$", mcols(ddsList[[i]])$results$class))) {
      next
    }
    tidied_data <- tidy_significant_dds(ddsList[[i]], mcols(ddsList[[i]])$results, tidy_fn)
    if (!first_done) {
      pdat <- tidied_data$pdat
      grouper <- setdiff(group_vars(pdat), ".gene")
      if (length(grouper)) {
        column_split=pdat[grouper]
      } else {
        column_split=NULL
      }
      pdat[sapply(pdat, is.character)] <- lapply(pdat[sapply(pdat, is.character)], 
                                                as.factor)
      pdat <- pdat[,get_terms(ddsList[[i]])]
      colList <- df2colorspace(pdat)
      first_done <- TRUE
    }
    colnames(tidied_data$mat) <- rownames(pdat)
    name <- sub(".*\\t", "", i)
    ha <- ComplexHeatmap::HeatmapAnnotation(df=pdat, col=colList)
    pl <- ComplexHeatmap::Heatmap(tidied_data$mat,
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
##' Generate MA Plots
##'
##' MA plots of each differential genelist
##' @title Generate MA Plots
##' @param ddsList DESdemonA-generated list of DESeqDataSets
##' @param caption The function to use to create captions
##' @return 
##' @author Gavin Kelly
#' @export
differential_MA <- function(ddsList, caption) {
  for (i in names(ddsList)) {
    res <- as.data.frame(mcols(ddsList[[i]])$results)
    md <- metadata(ddsList[[i]])
    pal <- RColorBrewer::brewer.pal(12, "Set3")
    if (is_formula(md$comparison)) {
      res$class <- ifelse(grepl("\\*", res$class),res$class, "None")
      colind <- (seq(along=unique(res$class), from=0)) %% length(pal) + 1
      yax <- "Largest fold change"
      if (length(colind)>12) {
        res$class <- ifelse(res$class=="None","None","Sig")
        cols <- setNames(c("blue","grey"), c("Sig","None"))
        lpos <- "none"
      } else {
        cols <- setNames(c(pal[colind]), unique(res$class))
        cols["None"] <- "grey"
        lpos <- "bottom"
      }
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

##' Heatmap colour-scheme generator
##'
##' For each column in a dataframe, generate a sensible colour palette
##' for each column
##' @title Heatmap colour-scheme generator
##' @param df
##' @return
##' @author Gavin Kelly
df2colorspace <- function(df) {
  pal <- RColorBrewer::brewer.pal(12, "Set3")
  df <- dplyr::mutate_if(as.data.frame(df), is.character, as.factor)
  purrr::map2(df,
              cumsum(c(0,purrr::map(df, nlevels)))[1:length(df)],
              function(column, start_level) {
                if (is.factor(column)) {
                  setNames(pal[(seq(start_level, length=nlevels(column)) %% length(pal)) + 1],
                           levels(column))
                } else {
                  circlize::colorRamp2(range(column), c("white","red"))
                }
              }
              )
}
