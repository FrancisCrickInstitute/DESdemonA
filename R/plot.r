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
  pc <- as.matrix(colData(dds)$.PCA)
  percentVar <- metadata(colData(dds)$.PCA)$percentVar
  is_vary <- sapply(colData(dds)[metadata(dds)$labels], function(v) length(unique(v))!=1)
  fml <- as.formula(paste0("~", paste(metadata(dds)$labels[is_vary], collapse="+")))
  plotFrame <- expand.grid(Covariate=all.vars(fml),
                          PC=1:ncol(pc))
  plotFrame$Assoc <- 0.0
  plotFrame$AIC <- 0.0
  for (ipc in 1:ncol(pc)) {
    fit0 <- lm(pc[,ipc]  ~ 1, data=colData(dds))
    fit1 <- add1(fit0, fml, test="Chisq")
#    covvar_PC[rownames(fit1)[-1],ipc] <- -log10(fit1$`Pr(>Chi)`[-1])
#    covvar_PC[rownames(fit1)[-1],ipc] <- 1-fit1$RSS[-1]/fit1$RSS[1]
    ind <- plotFrame$PC==ipc
    ind_fit <- match(plotFrame$Covariate[ind], row.names(fit1)[-1])
    aic <-  ifelse(fit1$RSS[1] > fit1$RSS[-1], 1, -1)[ind_fit]
    plotFrame$Assoc[ind] <- (1-(fit1$RSS[-1])/(fit1$RSS[1]))[ind_fit] * aic
  }
  plotFrame$PC <- sprintf("%02d", plotFrame$PC)
  qc_vis$PC_phenotype <- ggplot(plotFrame, aes(x=PC, y=Covariate, fill=Assoc)) +
    geom_raster() +
    scale_fill_gradient2(low="#4575b4", mid="grey90", high="#d73027") +
    theme_classic() + coord_fixed()

  print(qc_vis$PC_phenotype)
  caption("Covariate-PC association")
  qc_vis$PC <- list()
  cat(header, "# Visualisation of PCs ", pc_x, " and ", pc_y, " coloured by covariate", "\n", sep="")
  do_labels <- nrow(colDat)<10
  for (j in vars[is_vary]) {
    pc.df <- data.frame(PC1=pc[,pc_x], PC2=pc[,pc_y], col=colDat[[j]], sample=rownames(colDat))
    qc_vis$PC[[j]] <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=col))  + geom_point(size=3) +
      xlab(paste0("PC ", pc_x, ": ", percentVar[pc_x], "% variance")) +
      ylab(paste0("PC ", pc_y, ": ", percentVar[pc_y], "% variance")) +
      labs(colour=j)
    if (do_labels) {qc_vis$PC[[j]] <- qc_vis$PC[[j]] + geom_text_repel(aes(label=sample))}
    print(qc_vis$PC[[j]])
    caption(paste0("Coloured by ", j))
  }
  qc_vis$PC_assoc <- list()
  cat(header, "# Visualisation of associated PCs coloured by covariate", "\n", sep="") 
  for (j in vars[is_vary]) {
    tmp <- subset(plotFrame, Covariate==j)
    tmp <- tmp[order(tmp$Assoc, decreasing=TRUE),]
    pc_x <- as.integer(as.character(tmp$PC[1]))
    pc_y <- as.integer(as.character(tmp$PC[2]))
    pc.df <- data.frame(PC1=pc[,pc_x], PC2=pc[,pc_y], col=colDat[[j]], sample=rownames(colDat))
    qc_vis$PC[[j]] <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=col))  + geom_point(size=3) +
      xlab(paste0("PC ", pc_x, ": ", percentVar[pc_x], "% variance")) +
      ylab(paste0("PC ", pc_y, ": ", percentVar[pc_y], "% variance")) +
      labs(colour=j)
    if (do_labels) {qc_vis$PC[[j]] <- qc_vis$PC[[j]] + geom_text_repel(aes(label=sample))}

    print(qc_vis$PC[[j]])
    caption(paste0("Coloured by ", j))
  }

  invisible(qc_vis)
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
      pdat <- pdat[,metadata(ddsList[[i]])$labels]
      colList <- df2colorspace(pdat)
      first_done <- TRUE
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
              cumsum(c(0,map(df, nlevels)))[1:length(df)],
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
