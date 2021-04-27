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
#  dend <- dendro_all(plotDat, colDat[[c(vars$groups, vars$fixed)[1]]])
  qc_vis$heatmap <- ComplexHeatmap::Heatmap(
    plotDat, name="Mean Centred", column_title="Samples", row_title="Genes",
#    cluster_columns=dend,
    heatmap_legend_param = list(direction = "horizontal" ),
    col=colorspace::diverging_hcl(5, palette="Blue-Red"),
    top_annotation=ComplexHeatmap::HeatmapAnnotation(df=colDat[vars$fixed],
                                                     col = df2colorspace(colDat[vars$fixed])),
    show_row_names=FALSE, show_column_names=TRUE)
  draw(qc_vis$heatmap, heatmap_legend_side="top")
  caption("Heatmap of variable genes")

  cat(header, "# Heatmap of sample distances", "\n", sep="") 
  poisd <- PoiClaClu::PoissonDistance(t(counts(dds, normalized=TRUE)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  if (length(vars$fixed)>1) {
    rownames(samplePoisDistMatrix) <-Reduce(function(...) paste(..., sep="_"), colDat[vars$fixed])
  } else {
    rownames(samplePoisDistMatrix) <-colDat[[vars$fixed]]
  }
  if (is.null(vars$groups)) {
    colnames(samplePoisDistMatrix) <- row.names(colDat)
  } else {
    if (length(vars$groups)>1) {
      rownames(samplePoisDistMatrix) <-Reduce(function(...) paste(..., sep="_"), colDat[vars$groups])
    } else {
      rownames(samplePoisDistMatrix) <-colDat[[vars$groups]]
    }
  }
  qc_vis$sample_dist <- ComplexHeatmap::Heatmap(
    samplePoisDistMatrix,
    name="Poisson Distance", 
    clustering_distance_rows = function(x) {poisd$dd},
    clustering_distance_columns = function(x) {poisd$dd},
    heatmap_legend_param = list(direction = "horizontal" ),
    top_annotation=ComplexHeatmap::HeatmapAnnotation(df=colDat[vars$fixed],
                                     col = df2colorspace(colDat[vars$fixed]))
    )
  draw(qc_vis$sample_dist, heatmap_legend_side="top")
  caption("Heatmap of sample distances")
  
  ### PCA
  pc <- as.matrix(colData(dds)$.PCA)
  percentVar <- metadata(colData(dds)$.PCA)$percentVar
  is_vary <- sapply(colData(dds)[vars$fixed], function(v) length(unique(v))!=1)
  #fml <- as.formula(paste0("~", paste(metadata(dds)$labels[is_vary], collapse="+")))

  qc_vis$PC <- list(list())
  cat(header, "# Visualisation of PCs ", pc_x, " and ", pc_y, " coloured by covariate", "\n", sep="")
  do_labels <- nrow(colDat)<10
  pc_frame <- expand.grid(sample=rownames(colDat), PC=1:ncol(pc))
  pc_frame$coord <- as.vector(pc)
  pc_frame <- cbind(pc_frame, colDat)
                        
  for (j in vars$fixed[is_vary]) {
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
  models_for_qc <- sapply(metadata(dds)$models, "[[", "plot_qc") # TODO Make this more sensible
  for (model_name in names(metadata(dds)$models)[models_for_qc]) {
    qc_vis$PC[[model_name]] <- list()
    fml <- update(metadata(dds)$models[[model_name]]$design, paste(yvar, "~ ."))
    plotFrame <- expand.grid(Covariate=attr(terms(fml), "term.labels"),
                            PC=1:ncol(pc))
    plotFrame$Assoc <- NA
    npc <- ncol(pc)
    fit_selected <- list()
    for (ipc in 1:npc) {
      fitFrame[[yvar]] <- pc[,ipc]
      fit1 <-  lm(fml, data=fitFrame)
      fit_selected[[ipc]] <- MASS::stepAIC(fit1, trace=0)
      ind <- attr(terms(fit_selected[[ipc]]),"term.labels")
      if (length(ind)) {
        fit0 <- as.data.frame(anova(fit_selected[[ipc]]))
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
    
    cat(header, "# Partial Residuals of Associated PCs ",model_name, "\n", sep="")
    qc_vis$PC[[model_name]] <- list()
    factors_1 <- attr(terms(fit1),"factors")
    many_levels <- sapply(colDat, function(x) length(unique(x)))
    for (j in unique(plotFrame$Covariate)) {
      this_covar <- subset(plotFrame, Covariate==j & !is.na(Assoc) & PC!=PC[nrow(plotFrame)])
      if (nrow(this_covar)==0) next
      pc_first <- as.integer(as.character(this_covar$PC[1]))
      this_covar <- this_covar[order(this_covar$Assoc, decreasing=TRUE),]
      pc_strong <- as.integer(as.character(this_covar$PC[1]))
      pc.df <- subset(pc_frame, PC %in% c(pc_first, pc_strong))
      pc.df$label <- as.character(pc.df$PC)
      pc.df$label[pc.df$PC==pc_first] <- paste0(pc.df$label[pc.df$PC==pc_first], "+First")
      pc.df$label[pc.df$PC==pc_strong] <- paste0(pc.df$label[pc.df$PC==pc_strong], "+Best")
      pc.df$label <- paste0(pc.df$label, " (", percentVar[pc.df$PC], "%)")
      inter_vars <- rownames(factors_1)[factors_1[,j]!=0]
      inter_vars <- inter_vars[order(many_levels[inter_vars], decreasing=TRUE)]
      pc.df$X <- colDat[,inter_vars[1]]
      for (ipc in unique(c(pc_first, pc_strong))) {
        pc.df$coord[pc.df$PC==ipc] <- part.resid(fit_selected[[ipc]])[,j]
      }
      if (length(inter_vars)>1) {
        if (length(inter_vars)==2) {
          pc.df$col <-  colDat[,inter_vars[2]]
        } else {
          pc.df$col <- Reduce(interaction, colDat[,inter_vars[-1]])
        }
        my_aes <- ggplot2::aes(x=X, y=coord, colour=col, group=col)
        my_lbl <- labs(
          x=inter_vars[1],
          colour=paste(inter_vars[-1],collapse="x"),
          y="Residual")
      } else {
        my_aes <- ggplot2::aes(x=X, y=coord, group=1)
        my_lbl <- labs(colour=j,x=inter_vars[1], y="Residual")
      }
      qc_vis$PC[[model_name]][[j]] <- ggplot(pc.df, my_aes)  + geom_point(size=2) +
        stat_summary(fun = mean, geom="line") + 
        my_lbl + facet_wrap(~label, scales="free_y") 
      if (do_labels) {qc_vis$PC[[model_name]][[j]] <- qc_vis$PC[[model_name]][[j]] + geom_text_repel(aes(label=sample))}
      print(qc_vis$PC[[model_name]][[j]])
      caption(paste0("Against ", j))
    }
  }
  invisible(qc_vis)
}

get_terms <- function(dds) {
  ret <- list(fixed=NULL, groups=NULL)
  if ("full_model" %in% names(metadata(dds))) {
    ret <- classify_terms(metadata(dds)$full_model)
    return(ret)
  }
  if ("model" %in% names(metadata(dds))) {
    ret$fixed <- all.vars(metadata(dds)$model$design)
    return(ret)
  } else {
    term_list <- lapply(metadata(dds)$models, function(mdl) {all.vars(mdl$design)})
    ret$fixed <- unique(unlist(term_list))
  }
  return(ret)
}

part.resid <- function(fit) {
  pterms <- predict(fit, type="terms")
  apply(pterms,2,function(x)x+resid(fit))
}
  
  
  
  

dendro_all <- function(mat, var) {
  col_dend  <- as.dendrogram(hclust(dist(t(mat))))
  col  <-  df2colorspace(data.frame(v=var))[[1]]
  dendrapply(col_dend, function(n) {
    if (length(unique(var[unlist(n)]))==1) {
      attr(n, "edgePar")$col <- col[unlist(n)[1]]
      if (!is.leaf(n)) {
        attr(n, "nodePar")$lwd <- 2
      }
    }
    n
  })
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
differential_heatmap <- function(ddsList, tidy_fn=NULL, caption, colList=df2colorspace(colData(ddsList[[1]]))) {
  first_done <- FALSE
  for (i in names(ddsList)) {
    if (!any(grepl("\\*$", mcols(ddsList[[i]])$results$class))) {
      next
    }
    comp <- metadata(ddsList[[i]])$comparison
    if ("spec" %in% names(attributes(comp))) {
      tidy_fn <- emmeans:::.parse.by.formula(attr(comp, "spec"))
      tidy_fn$rhs <- c(tidy_fn$rhs, setdiff(all.vars(metadata(ddsList[[i]])$model$design), unlist(tidy_fn)))
    } else {
      tidy_fn <- list(lhs="", rhs=all.vars(metadata(ddsList[[i]])$model$design), by=NULL)
    }
    tidied_data <- tidy_significant_dds(ddsList[[i]], mcols(ddsList[[i]])$results, tidy_fn)
    pdat <- tidied_data$pdat
    pdat[sapply(pdat, is.character)] <- lapply(pdat[sapply(pdat, is.character)], 
                                              as.factor)
    colnames(tidied_data$mat) <- rownames(pdat)
    name <- sub(".*\\t", "", i)
    if (length(tidy_fn$by)==0) {
      col_split <- NULL
    } else {
      col_split <- apply(sapply(pdat[, tidy_fn$by, drop=FALSE], as.character), 1, paste, collapse=" ")
    }
    ha <- ComplexHeatmap::HeatmapAnnotation(df=pdat, col=colList)
    pl <- ComplexHeatmap::Heatmap(tidied_data$mat,
                 heatmap_legend_param = list(direction = "horizontal" ),
                 name=sub(".*\\t", "", i),
                 cluster_columns = FALSE,
                 show_column_names = TRUE,
                 column_split = col_split,
                 top_annotation = ha,
                 row_names_gp = gpar(fontsize = 6),
                 show_row_names = nrow(tidied_data$mat)<100)
    draw(pl, heatmap_legend_side="top")
    caption(paste0("Heatmap on differential genes ", name))
    part_resid <- residual_heatmap_transform(tidied_data$mat, pdat, metadata(ddsList[[i]])$model$design)
    for (i_term in 1:(dim(part_resid)[3])) {
      term_name <- dimnames(part_resid)[i_term]
          pl <- ComplexHeatmap::Heatmap(part_resid[,,i_terms],
                 heatmap_legend_param = list(direction = "horizontal" ),
                 name=paste(sub(".*\\t", "", i),term_name),
                 cluster_columns = FALSE,
                 show_column_names = TRUE,
                 column_split = col_split,
                 top_annotation = ha,
                 row_names_gp = gpar(fontsize = 6),
                 show_row_names = nrow(tidied_data$mat)<100)
          draw(pl, heatmap_legend_side="top")
          caption(paste0("Heatmap on differential genes ", name, ", ", term_name, "-focussed"))
    }
  }
}

residual_heatmap_transform <- function(mat, cdata, fml) {
  tmat <- t(mat)
  fml <- update(fml, tmat ~ .)
  fit <- lm(fml, data=cdata)
  fit1 <- fit
  class(fit1) <- "lm"
  ind <- c("coefficients","residuals","effects","fitted.values")
  fit1[ind] <- lapply(fit[ind], function(x) x[,1])
  p1 <- predict(fit1, type="terms")
  out <- array(0, c(dim(fit$residuals), ncol(p1)), dimnames=c(dimnames(fit$residuals), list(colnames(p1))))
  for (i in 1:(dim(out)[2])) {
    fit1[ind] <- lapply(fit[ind], function(x) x[,i])
    out[,i,] <- predict(fit1, type="terms")
  }
  for (j in 1:(dim(out)[3])) {
    out[,,j] <- out[,,j]+fit$residuals
  }
  out
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
