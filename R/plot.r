sym_colour <- function(dat, lo="blue",zero="white", hi="red") {
  mx <- max(abs(dat))
  circlize::colorRamp2(c(-mx, 0, mx), colors=c(lo, zero, hi))
}


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
  var_stab <- assay(dds, "vst")
  top <- order(apply(var_stab, 1, sd), decreasing=TRUE)[1:n] 
  models_for_qc <- metadata(dds)$models[sapply(metadata(dds)$models, "[[", "plot_qc")] # TODO Make this more sensible
  ### Heatmap
  cat(header, "# Heatmap of variable genes", "\n", sep="") 
  plotDat <- var_stab[top,]
  plotDat <- plotDat-rowMeans(plotDat)
  vars <- get_terms(dds)
  colDat <- as.data.frame(colData(dds))
  colnames(plotDat) <- rownames(colDat)
  #  dend <- dendro_all(plotDat, colDat[[c(vars$groups, vars$fixed)[1]]])
  pl <- ComplexHeatmap::Heatmap(
    plotDat, name="Mean Centred", column_title="Samples", row_title="Genes",
    col=sym_colour(plotDat),
    #    cluster_columns=dend,
    heatmap_legend_param = list(direction = "horizontal" ),
#    col=colorspace::diverging_hcl(5, palette="Blue-Red"),
#    col = circlize::colorRamp2(sym_colour(plotDat), colors=c("blue", "white", "red")),
    top_annotation=ComplexHeatmap::HeatmapAnnotation(df=colDat[vars$fixed],
                                                     col = metadata(colData(dds))$palette[vars$fixed]),
    show_row_names=FALSE, show_column_names=TRUE)
  draw(pl, heatmap_legend_side="top")
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
  pl <- ComplexHeatmap::Heatmap(
    samplePoisDistMatrix,
    col = circlize::colorRamp2(c(0,max(samplePoisDistMatrix)), colors=c("white", "red")),
    name="Poisson Distance", 
    clustering_distance_rows = function(x) {poisd$dd},
    clustering_distance_columns = function(x) {poisd$dd},
    heatmap_legend_param = list(direction = "horizontal" ),
    top_annotation=ComplexHeatmap::HeatmapAnnotation(df=colDat[vars$fixed],
                                                     col = metadata(colData(dds))$palette[vars$fixed]
                                                     )
  )
  draw(pl, heatmap_legend_side="top")
  caption("Heatmap of sample distances")
  
  ### PCA
  pc <- as.matrix(colData(dds)$.PCA)
  percentVar <- metadata(colData(dds)$.PCA)$percentVar
  is_vary <- sapply(colData(dds)[vars$fixed], function(v) length(unique(v))!=1)
  #fml <- as.formula(paste0("~", paste(metadata(dds)$labels[is_vary], collapse="+")))
  
  
  do_focus <- length(vars$fixed)>1 # ie also remove effect of every other covariate
  do_batch <- FALSE                # ie just remove effect of specified batch covariates
  if (do_focus) {
    sample_gene_factor <- residual_heatmap_transform(
      assay(dds, "vst"),
      colData(dds),
      models_for_qc[[1]]$design)
    if ("batch" %in% names(models_for_qc[[1]])) {
      batch_vars <- all.vars(models_for_qc[[1]]$batch)
      pc_batch <- prcomp(t(assay(dds, "vst")) - apply(sample_gene_factor$terms[,,batch_vars, drop=FALSE], 1:2,sum),
                  scale=FALSE)
      pc_resid <- list(pc=pc_batch$x, percent=round(100 * pc_batch$sdev^2 / sum( pc_batch$sdev^2 )))
      do_focus <- FALSE
      do_batch <- TRUE
    } else {
      pc_resid <- lapply(
        dimnames(sample_gene_factor$terms)[[3]],
        function(fac) {
          pc <- prcomp(sample_gene_factor$terms[,,fac,drop=TRUE] + sample_gene_factor$resid,
                      scale=FALSE) 
          list(pc=pc$x, percent=round(100 * pc$sdev^2 / sum( pc$sdev^2 )))
        }
      )
      names(pc_resid) <- dimnames(sample_gene_factor$terms)[[3]]
    }
  }
  
  cat(header, "# Visualisation of PCs ", pc_x, " and ", pc_y, " coloured by covariate", "\n", sep="")
  do_labels <- nrow(colDat)<10
  
  for (j in vars$fixed[is_vary]) {
    pc.df <- data.frame(PC1=pc[,pc_x], PC2=pc[,pc_y], col=colDat[[j]], sample=rownames(colDat))
    pl <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=col))  + geom_point(size=3) +
      xlab(paste0("PC ", pc_x, ": ", percentVar[pc_x], "% variance")) +
      ylab(paste0("PC ", pc_y, ": ", percentVar[pc_y], "% variance")) +
      scale_colour_manual(values=metadata(colData(dds))$palette[[j]]) + 
      labs(colour=j)
    if (do_labels) {pl <- pl + geom_text_repel(aes(label=sample))}
    print(pl)
    caption(paste0("Coloured by ", j))
    if (do_batch) {
      pc.df <- data.frame(PC1=pc_resid$pc[,pc_x], PC2=pc_resid$pc[,pc_y], col=colDat[[j]], sample=rownames(colDat))
      pl <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=col))  + geom_point(size=3) +
        xlab(paste0("PC ", pc_x, ": ", pc_resid$percent[pc_x], "% variance")) +
        ylab(paste0("PC ", pc_y, ": ", pc_resid$percent[pc_y], "% variance")) +
        scale_colour_manual(values=metadata(colData(dds))$palette[[j]]) + 
        labs(colour=j)
      if (do_labels) {pl <- pl + geom_text_repel(aes(label=sample))}
      print(pl)
      caption(paste0("Coloured by ", j, ", correcting for ", paste(batch_vars, collapse=", ")))
    }
    if (do_focus && j %in% names(pc_resid)) {
      pc.df$PC1 <- pc_resid[[j]]$pc[,pc_x]
      pc.df$PC2 <- pc_resid[[j]]$pc[,pc_y]
      pl <- ggplot(pc.df, aes(x=PC1, y=PC2, colour=col))  +
        annotate("text", x = Inf, y = -Inf, label = "PROOF ONLY",
                 hjust=1.1, vjust=-1.1, col="white", cex=6,
                 fontface = "bold", alpha = 0.8) + 
        geom_point(size=3) +
        scale_colour_manual(values=metadata(colData(dds))$palette[[j]]) + 
        xlab(paste0("PC ", pc_x, ": ", pc_resid[[j]]$percent[pc_x], "% variance")) +
        ylab(paste0("PC ", pc_y, ": ", pc_resid[[j]]$percent[pc_y], "% variance")) +
        labs(colour=j)
      if (do_labels) {pl <- pl + geom_text_repel(aes(label=sample))}
      print(pl)
      caption(paste0("Focussed on ", j))
    }
  }
  
  pc_frame <- expand.grid(sample=rownames(colDat), PC=1:ncol(pc))
  pc_frame$coord <- as.vector(pc)
  pc_frame <- cbind(pc_frame, colDat)
  fitFrame <- colDat
  yvar <- make.unique(c(colnames(fitFrame), "y", sep=""))[ncol(fitFrame)+1]
  for (model_name in names(models_for_qc)) {
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
    pl <- ggplot(plotFrame, aes(x=PC, y=Covariate, fill=Assoc)) +
      geom_raster() +
      facet_wrap(~wrap, scales="free_x", ncol=1) + 
      scale_fill_gradient2(low="#4575b4", mid="grey90", high="#d73027") +
      theme_classic() + theme(aspect.ratio = length(unique(plotFrame$Covariate)) / min(20, npc))
    print(pl)
    caption("Covariate-PC association")
    
    cat(header, "# Partial Residuals of Associated PCs ",model_name, "\n", sep="")
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
        pc.df$col <-  colDat[,inter_vars[2]]
        if (length(inter_vars)==2) {
          pl <- ggplot(pc.df, aes(x=X, y=coord, colour=col, group=col)) +
            scale_colour_manual(values=metadata(colData(dds))$palette[[inter_vars[2]]]) +
            stat_summary(fun = mean, geom="line") + 
            labs(
              x=inter_vars[1],
              colour=paste(inter_vars[-1],collapse="x"),
              y="Residual")
        } else {
          pc.df$grp <- Reduce(interaction, colDat[,inter_vars[-1]])
          pl <- ggplot(pc.df, aes(x=X, y=coord, colour=col, group=grp)) +
            scale_colour_manual(values=metadata(colData(dds))$palette[[inter_vars[2]]]) + 
            stat_summary(fun = mean, geom="line") + 
            labs(
              x=inter_vars[1],
              colour=paste(inter_vars[-1],collapse="x"),
              y="Residual")
        }
      } else {
        pl <- ggplot(pc.df, aes(x=X, y=coord, colour=X, group=1)) + 
          scale_colour_manual(values=metadata(colData(dds))$palette[[inter_vars[1]]]) +
          stat_summary(fun = mean, geom="line", colour="black") + 
          labs(colour=j,x=inter_vars[1], y="Residual")
      }
      pl <- pl  + geom_point(size=2) +
        facet_wrap(~label, scales="free_y") 
      if (do_labels) {pl <- pl + geom_text_repel(aes(label=sample))}
      print(pl)
      caption(paste0("Against ", j))
    }
  }
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
differential_heatmap <- function(ddsList, tidy_fn=NULL, caption) {
  first_done <- FALSE
  for (i in names(ddsList)) {
    if (!any(grepl("\\*$", mcols(ddsList[[i]])$results$class))) {
      next
    }
    comp <- metadata(ddsList[[i]])$comparison
    fml <- metadata(ddsList[[i]])$model$design
    if ("spec" %in% names(attributes(comp))) {
      var_roles <- emmeans:::.parse.by.formula(attr(comp, "spec"))
      var_roles$all <- c(var_roles$by, var_roles$rhs, setdiff(all.vars(fml), unlist(var_roles)))
      mdl <- metadata(ddsList[[i]])$model
      if ("mat" %in% names(mdl)) {
        mmat <- mdl$mat
      } else {
        mmat <- model.matrix(mdl$design, colData(ddsList[[i]]))
      }
      weights <- t(metadata(ddsList[[i]])$comparison %*% MASS::ginv(mmat))
      is_denom <- weights < -sqrt(.Machine$double.eps)
      weights[!is_denom] <- 0
    } else {
      var_roles <- list(lhs="", rhs=all.vars(fml), by=NULL, all=all.vars(fml))
      weights <- NULL
    }
    tidied_data <- tidy_significant_dds(ddsList[[i]], mcols(ddsList[[i]])$results, var_roles, weights=weights)
    pdat <- tidied_data$pdat
    pdat[sapply(pdat, is.character)] <- lapply(pdat[sapply(pdat, is.character)], 
                                              as.factor)
    colnames(tidied_data$mat) <- rownames(pdat)
    name <- sub(".*\\t", "", i)
    if (length(var_roles$by)==0) {
      col_split <- NULL
    } else {
      col_split <- apply(sapply(pdat[, var_roles$by, drop=FALSE], as.character), 1, paste, collapse=" ")
    }
    
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df=pdat,
      col=metadata(colData(ddsList[[i]]))$palette[names(pdat)])
    pl <- ComplexHeatmap::Heatmap(
      tidied_data$mat,
      col=sym_colour(tidied_data$mat),
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
    if (length(all.vars(fml))>1) {
      part_resid <- residual_heatmap_transform(tidied_data$mat, pdat, fml)
      term_names <- intersect(dimnames(part_resid$terms)[[3]], var_roles$rhs)
      for (term_name in term_names) {
        tdat <- t(part_resid$terms[,,term_name, drop=TRUE]) +part_resid$const + t(part_resid$resid)
        pl <- ComplexHeatmap::Heatmap(
          tdat,
          col=sym_colour(tdat),
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
}

residual_heatmap_transform <- function(mat, cdata, fml) {
  assign("tmat", t(mat), envir=environment(fml))
  fml <- update(fml, tmat ~ .)
  fit <- lm(fml, data=cdata)
  fit1 <- fit
  class(fit1) <- "lm"
  ind <- c("coefficients","residuals","effects","fitted.values")
  fit1[ind] <- lapply(fit[ind], function(x) x[,1])
  p1 <- predict(fit1, type="terms")
  out <- array(0, c(dim(fit$residuals), ncol(p1)), dimnames=c(dimnames(fit$residuals), list(colnames(p1))))
  const <- numeric(dim(out)[2])
  for (i in 1:(dim(out)[2])) {
    fit1[ind] <- lapply(fit[ind], function(x) x[,i])
    pred <- predict(fit1, type="terms")
    out[,i,] <- predict(fit1, type="terms")
    const[i] <- attr(pred, "constant")
  }
  list(terms=out, const=const, resid=fit$residuals)
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
df2colorspace <- function(df, palette) {
  pal <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[palette, "maxcolors"], palette)
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
