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
##' @param ... The parent factors that 'inner' is inside, e.g. the genotype of the cell-line
##' @return A recoded factor that can now replace the inner (cell-line) variable in your model
##' @author Gavin Kelly
#' @export
load_specs <- function(file="", context) {
  if (file.exists(file)) {
    e <- as.environment(as.data.frame(colData(context)))
    parent.env(e) <- environment()
    assign("sample_set", list, envir=e)
    assign("model", list, envir=e)
    assign("specification", list, envir=e)
    assign("settings",
           function(...) {
             as.list( substitute(alist(...)))[-1]
           },
           envir=e
           )
    assign("mutate",
           function(...) {
             substitute(alist(...))
           },
           envir=e
           )
    specs <- source(file, local=e)$value
    assign("sample_set", expression, envir=e) # avoid evaluating any examples sample_sets.
    pkg_defaults <- source(system.file("templates/example.spec", package="DESdemonA"), local=e)$value
    new_settings <- setdiff(names(pkg_defaults$settings), names(specs$settings))
    if (length(new_settings)>0) {
      string_rep <- lapply(pkg_defaults$settings[new_settings], deparse)
      warning("New settings (", paste(new_settings), ") can be set in ", file, ", so please update it. The default values that will be used are:\n", paste(names(string_rep), string_rep, sep=": ", collapse="\n"))
      specs$settings[new_settings] <- pkg_defaults$settings[new_settings]
    }
    rm(list=ls(envir=e), envir=e)
  } else {
    fml <- paste("~", names(colData(dds))[ncol(colData(dds))])
    specs <- list(
      sample_sets = list(all=TRUE),
      models=list(
        "Naive" = list(
          design = as.formula(fml),
          comparisons = mult_comp(as.formula(paste("pairwise", fml)))
        )
      ),
      plot_scale = function(y) {
        y/mean(y)
      }
    )
  }
  specs
}

##' Recode nested factors to avoid matrix-rank problems
##'
##' Following the suggestion in the DESeq2 vignette, recode nested factors so that
##' they take common values in different clusters to avoid rank problems
##' @title Recode nested factors
##' @param inner A factor representing the variable to be recoded, e.g. the cell-line
##' @param ... The parent factors that 'inner' is inside, e.g. the genotype of the cell-line
##' @return A factor with recoded levels
##' @author Gavin Kelly
##' @export
recode_within <- function(inner, ...) {
  within <- do.call(interaction, alist(...))
  tab <- table(inner, within)!=0 # which batches are in which nest
  if (any(rowSums(tab)>1)) {
    stop("Some inner levels appear in multiple outer groups.")
  }
  factor(apply(tab, 2, cumsum)[cbind(as.character(inner),as.character(within))]) # cumsum to get incrementing index within group.
}



##' Expand an analysis specification into its corresponding subset list
##'
##' Generate a list of DESeq2 objects corresponding to the different
##' subsets specified
##' @title Generate subsets of DESeq2 object
##' @param dds The original DESeq2 object containing all samples
##' @param spec The analysis specificiation
##' @return A list of DESeq2 objects
##' @author Gavin Kelly
##' @export
build_dds_list <- function(dds, spec) {
  flat <- unlist(spec$sample_sets)
  flat <- flat[sapply(flat, is_formula)]
  modelled_terms <-  unlist(lapply(flat, all.vars))
  if (!"palette" %in% names(spec$settings)) {
    spec$settings$palette="Set1"
  }
  if (is.list(spec$settings$palette)) {
    default_palette <- spec$settings$palette
  } else {
    default_palette<- DESdemonA:::df2colorspace(
      colData(dds)[, intersect(modelled_terms, colnames(colData(dds))),drop=FALSE],
      spec$settings$palette
    )
  }
  metadata(colData(dds))$palette <- default_palette
  ddsList <- lapply(spec$sample_sets, function(set) {
    mdlList <- spec$models
    if (is.list(set)) {
      ind <- set$subset
      mdlList <- c(mdlList, set$models)
    } else {
      ind <- set
    }
    obj <- dds[,ind]
    metadata(obj)$full_model <- spec$full_model
    colData(obj) <- droplevels(colData(obj))
    mdlList <- lapply(mdlList, function(x) modifyList(list(plot_qc=FALSE), x))
    if (!any(sapply(mdlList, "[[", "plot_qc"))) {
      mdlList[[1]]$plot_qc <- TRUE
    }
    metadata(obj)$models <- mdlList
    if ("transform" %in% names(set)) {
      .mu <- purrr::partial(mutate, .data=as.data.frame(colData(obj)))
      tr <- set$transform
      tr[[1]] <- .mu
      cnames <- colnames(obj)
      colData(obj) <- S4Vectors::DataFrame(eval(tr))
      metadata(colData(obj)) <- metadata(colData(dds))
      colnames(obj) <- cnames
      new_cols <- intersect(modelled_terms, setdiff(names(colData(obj)), names(default_palette$Heatmap)))
      old_cols <- setdiff(modelled_terms, new_cols)
      if (length(old_cols)>0) {
        is_modified <- sapply(old_cols,
                               function(x) {
                                 if (class(colData(obj)[[x]]) != class(colData(dds)[[x]])) return(TRUE)
                                 if (is.factor(colData(obj)[[x]])) return(!all(levels(colData(obj)[[x]]) %in%  levels(colData(dds)[[x]])))
                                 if (is.character(colData(obj)[[x]])) return(!all(unique(colData(obj)[[x]]) %in%  levels(unique(dds)[[x]])))
                                 return(!all(range(colData(obj)[[x]])==range(colData(obj)[[x]])))
                               })
        new_cols <- c(new_cols, old_cols[is_modified])
      }
      if (length(new_cols)>0) {
        new_meta <- DESdemonA:::df2colorspace(
          colData(obj)[, new_cols, drop=FALSE],
          spec$settings$palette
        )
        metadata(colData(obj))$palette$Heatmap[new_cols] <- new_meta$Heatmap[new_cols]
        metadata(colData(obj))$palette$ggplot[new_cols] <- new_meta$ggplot[new_cols]
      }
    }
    if ("collapse" %in% names(set)) {
      mf <- model.frame(set$collapse, data.frame(colData(obj)))
      ind <- match(do.call(paste, c(mf, sep="\r")),
                  do.call(paste, c(unique(mf), sep="\r")))
      obj <- collapseReplicates(obj, groupby=factor(ind), renameCols=FALSE)
    }
    obj
  })
  ddsList <- imap(ddsList,
                   function(obj, dname) {
                     metadata(obj)$dmc <- list(dataset=dname)
                     obj
                   }
                   )
}

##' Calculate dimension reduction 
##'
##' Add a vst transformed assay, and a projection of the samples ont PCA space
##' @title Store dimension-reduction results in DESeq2 object
##' @param dds The original DESeq2 object containing all samples
##' @param n 
##' @param family 
##' @param batch 
##' @param spec The analysis specificiation
##' @return 
##' @author Gavin Kelly
##' @export
add_dim_reduct  <-  function(dds, n=Inf, family="norm", batch=~1) {
  var_stab <- assay(vst(dds))
  if (batch != ~1) {
    var_stab <- residuals(limma::lmFit(var_stab, model.matrix(batch, as.data.frame(colData(dds)))), var_stab)
  }
  colnames(var_stab) <- colnames(dds)
  assay(dds, "vst") <- var_stab
  if (family=="norm") {
    pc <- prcomp(t(var_stab), scale=FALSE)
    percentVar <- round(100 * pc$sdev^2 / sum( pc$sdev^2 ))
    colData(dds)$.PCA <- DataFrame(pc$x)
    metadata(colData(dds)$.PCA)$percentVar <- setNames(percentVar, colnames(pc$x))
    mcols(dds)$PCA <-DataFrame(pc$rotation)
  } else {
    co <- counts(dds, norm=FALSE)
    pc_glm <- glmpca::glmpca(Y=co[rowSums(co)!=0,],
                            L=ncol(co),
                            fam=family,
                            X=if(batch == ~1) 
                              NULL
                            else
                              model.matrix(batch, as.data.frame(colData(dds)))
                            )
    colData(dds)$.PCA <- DataFrame(pc_glm$factors)
    metadata(colData(dds)$.PCA)$percentVar <- setNames(rep(0, ncol(co)), colnames(pc$x))
  }
  dds
}

##' Fit the models of expression
##'
##' Iterate through each model (stored in the 'models' metadata of a
##' DESeqDataSet) and expand the contrasts so each contrast gets a
##' separate nested level.
##' @title Fit the DESeq2 models
##' @param dds The original DESeq2 object containing all samples
##' @param ...
##' @return
##' @author Gavin Kelly
##' @export
fit_models <- function(dds, ...) {
  model_comp <- lapply(
    metadata(dds)$models,
    function(mdl) fit_model(mdl, dds)
  )
  model_comp <- model_comp[sapply(model_comp, length)!=0]
  model_comp <- imap(model_comp, function(obj, mname) {
    lapply(obj, function(y) {metadata(y)$dmc$model <- mname; y})
  })
  model_comp
}

fit_model <- function(mdl, dds) {
  this_dds <- dds
  design(this_dds) <- mdl$design
  metadata(this_dds)$model <- mdl
  this_dds <- check_model(this_dds) 
  out <- list()
  is_lrt <- sapply(mdl$comparisons, is_formula)
  if (any(!is_lrt)) {
    comps <- mdl$comparisons[!is_lrt]
    is_post_hoc <- sapply(comps, class)=="post_hoc"
    if (any(is_post_hoc)) {
      comps[is_post_hoc] <- lapply(
        comps[is_post_hoc],
        function(ph) {emcontrasts(dds=this_dds, spec=ph$spec, extra=ph[-1])}
      )
      comps[!is_post_hoc] <- lapply(comps[!is_post_hoc], list) # protect existing lists from unlist
      comps <- unlist(comps, recursive=FALSE)
    }
    if (any(metadata(this_dds)$model$dropped)) {
      design(this_dds) <- metadata(this_dds)$model$mat
    }
    this_dds <- DESeq2::DESeq(this_dds, test="Wald", ...)
    metadata(this_dds)$models <- NULL
    metadata(this_dds)$comparisons <- NULL
    out <- lapply(comps, function(cntr) {
      metadata(this_dds)$comparison <- cntr
      this_dds})
  }
  if (any(is_lrt)) {
    lrt <- lapply(mdl$comparisons[is_lrt],
                 function(reduced) {DESdemonA:::fitLRT(this_dds, mdl=mdl, reduced=reduced, ...)}
                 )
    out <- c(out, lrt)
  }
  out <- imap(out, function(obj, cname) {metadata(obj)$dmc$comparison <- cname; obj})
  out
}

##' Check model
##'
##' Run formula through an lm to check it
##' @title Check model
##' @param mdl 
##' @param coldat 
##' @param dds The original DESeq2 object containing all samples
##' @return 
##' @author Gavin Kelly
##' @export
check_model <- function(dds) {
  mdl <- metadata(dds)$model
  mdl$dropped <- FALSE
  if (is_formula(mdl$design) ) {
    df <- as.data.frame(colData(dds))
    df$.x <- counts(dds, norm=TRUE)[1,]
    fml <- as.formula(paste0(".x ~ ", as.character(design(dds)[2])))
    fit <- lm(fml, data=df)
    mdl$lm <- fit
    if ("drop_unsupported_combinations" %in% names(mdl) && mdl$drop_unsupported_combinations==TRUE) {
      mdl$dropped <- is.na(coef(fit))
    } else {
      if (any(is.na(coef(fit)))) {
        warning("Can't estimate some coefficients in ", mdl$design, ".\n In unbalanced nested designs, use the option 'drop_unsupported_combinations=TRUE,' in the problematic model. Other causes are complete confounding or conditions with no observations.")
      }
    }
  }
  if (any(mdl$dropped)) {
    mm <- model.matrix(mdl$design, as.data.frame(colData(dds)))[,!mdl$dropped]
    colnames(mm) <- DESdemonA:::.resNames(colnames(mm))
    mdl$mat <- mm
  }
  metadata(dds)$model <- mdl
  dds
}


##' Post-hoc generator
##'
##' Wrap a formula so that emmeans can auto-expand it
##' @title Mark a formula as a multiple comparison
##' @param spec 
##' @param ... 
##' @param dds The original DESeq2 object containing all samples
##' @return 
##' @author Gavin Kelly
##' @export
mult_comp <- function(spec, ...) {
  obj <- list(spec=spec,...)
  class(obj) <- "post_hoc"
  obj
}

##' Expand post-hoc comparisons
##'
##' Use emmeans to expand keywords
##' @title Expand multiple comparisons into their contrasts
##' @param dds The original DESeq2 object containing all samples
##' @param spec 
##' @return 
##' @author Gavin Kelly
##' @export
emcontrasts <- function(dds, spec, extra=NULL) {
  if ("keep" %in% names(extra)) {
    keep <- extra$keep
    extra$keep <- NULL
  } else {
    keep <- NA
  }
  mdl <- metadata(dds)$model
  emfit <- do.call(emmeans::emmeans, c(list(object=mdl$lm, specs= spec),extra))
  contr_frame <- as.data.frame(summary(emfit$contrasts))
  ind_est <- !is.na(contr_frame$estimate)
  if (!is.na(keep[1])) {
    ind_est  <- ind_est & contr_frame$contrast %in% keep
  }
  contr_frame <- contr_frame[ind_est,1:(which(names(contr_frame)=="estimate")-1), drop=FALSE]
  contr_mat <- emfit$contrast@linfct[ind_est, !mdl$dropped, drop=FALSE]
  colnames(contr_mat) <- DESdemonA:::.resNames(colnames(contr_mat))
  contr <- lapply(seq_len(nrow(contr_frame)), function(i) contr_mat[i,,drop=TRUE])
  contr <- lapply(contr, function(vect) {attr(vect, "spec") <- spec; vect})
  names(contr) <- do.call(paste, c(contr_frame,sep= "|"))
  contr
}


##' Fit an LRT model
##'
##' Insert the 'comparison' formula into the reduced slot
##' @title Fit LRT
##' @param dds The original DESeq2 object containing all samples
##' @param reduced 
##' @param ... 
##' @return 
##' @author Gavin Kelly
fitLRT <- function(dds, mdl, reduced, ...) {
  mdl <- metadata(dds)$model
  metadata(dds)$comparison <- reduced
  if (any(mdl$dropped)) {
    full <- mdl$mat
    reduced <- model.matrix(reduced, colData(dds))
    ## unsupported_ind <- apply(reduced==0, 2, all)
    ## reduced <- reduced[, !unsupported_ind]
    ## colnames(reduced) <- DESdemonA:::.resNames(colnames(reduced))
    reduced <- reduced[,colnames(reduced) %in% colnames(full)]
    metadata(dds)$reduced_mat <- reduced
  } else {
    full <- mdl$design
  }
  design(dds) <- full
  dds <- DESeq2::DESeq(dds, test="LRT", full=full, reduced=reduced, ...)
  metadata(dds)$LRTterms=setdiff(
    colnames(attr(dds, "modelMatrix")),
    colnames(attr(dds, "reducedModelMatrix"))
  )
  metadata(dds)$models <- NULL
  metadata(dds)$comparisons <- NULL
  dds
}

## apply contrast, and transfer across interesting mcols from the dds
##' Generate results object
##'
##' Insert results columns into mcols
##' @title Generate the results for a model and comparison
##' @param dds The original DESeq2 object containing all samples
##' @param mcols 
##' @param filterFun 
##' @param lfcThreshold 
##' @param ... 
##' @return 
##' @author Gavin Kelly
##' @export
get_result <- function(dds, mcols=c("symbol", "entrez"), filterFun=IHW::ihw, lfcThreshold=0, alpha=0.1, ...) {
  if (is.null(filterFun)) filterFun <- rlang::missing_arg()
  comp <- metadata(dds)$comparison
  if (length(alpha)>1) {
    alpha <- sort(alpha)
    alpha1 <- alpha[1]
  } else {
    alpha1 <- alpha
  }
  if (!is_formula(comp)) {
    if (is.character(comp) && length(comp)==1) { #  it's a name
      r <- DESeq2::results(dds, filterFun=filterFun, lfcThreshold=lfcThreshold, name=metadata(dds)$comparison, alpha=alpha1, ...)
    } else { # it's a contrast
      if (is.list(comp) && "listValues" %in% names(comp)) {
        r <- DESeq2::results(dds, filterFun=filterFun, lfcThreshold=lfcThreshold, contrast=metadata(dds)$comparison[names(comp) != "listValues"], listValues=comp$listValues, alpha=alpha1, ...)
      } else {
        r <- DESeq2::results(dds, filterFun=filterFun, lfcThreshold=lfcThreshold, contrast=metadata(dds)$comparison, alpha=alpha1, ...)
      }
    }
  } else { # it's LRT
    r <- results(dds, filterFun=filterFun, alpha=alpha1, ...)
  }
  sigs <- rep("NS", nrow(r))
  sigs[r$padj <= alpha1] <- paste("<=", alpha1)
  prev_alpha <- alpha1
  for (a in alpha[-1]) {
    if (missing(filterFun)) {
      res_alpha <- DESeq2:::pvalueAdjustment(r, independentFiltering=TRUE, alpha=a, pAdjustMethod="BH")$padj
    }
    else {
      res_alpha <- filterFun(r, alpha=a)$padj
    }
    ind <- sigs=="NS" & res_alpha <= a
    sigs[ind] <- paste0(prev_alpha, "-", a)
    prev_alpha <- a
  }
  my_mcols <- intersect(mcols, colnames(mcols(dds)))
  if (length(my_mcols)>0) {
    r[my_mcols] <- mcols(dds)[my_mcols]
  }
  r$sig <- sigs
  if ("LRTPvalue" %in% names(mcols(dds))) {
    r$class <- mcols(dds)$class
    r$class[is.na(r$padj) | is.na(r$pvalue) | r$baseMean==0] <- NA
    term <-  metadata(dds)$LRTterms
    # take the biggest fold-change vs baseline, for MA and reporting?
    if (all(term %in% names(mcols(dds)))) {
      effect_matrix <- cbind(I=rep(0, nrow(dds)),as.matrix(mcols(dds)[,term,drop=FALSE]))
      split_effects <- strsplit(term, "_vs_")
      # if a main effect is dropped, then all the terms are probably named A vs (intercept)
      # we can make the class a bit more interpretable by renaming the case where e.g.
      # B vs A is the largest and C vs A the smallest as B v C. So change effect columns to
      # enable this.
      if (all(sapply(split_effects, length)==2)) {
        if (length(unique(sapply(split_effects, "[", 2)))==1) {
          colnames(effect_matrix) <- c(split_effects[[1]][2], sapply(split_effects, "[", 1))
        }
      }
      maxmin <- cbind(
        apply(effect_matrix, 1, which.max),
        apply(effect_matrix, 1, which.min))
      #imax imin between them locate the max and min. imin is the 'earlier' term, to allow for negative and positive fc's
      imax <- apply(maxmin, 1, max)
      imin <- apply(maxmin, 1, min)
      r$maxLog2FoldChange <- effect_matrix[cbind(1:length(imax), imax)] -
        effect_matrix[cbind(1:length(imin), imin)]
      maxlfcSE <- sqrt(
      (as.matrix(mcols(dds)[,paste0("SE_", c("Intercept", term))])[cbind(1:length(imax), imax)])^2 +
        (as.matrix(mcols(dds)[,paste0("SE_", c("Intercept", term))])[cbind(1:length(imin), imin)])^2
      )
      fit <- ashr::ash(r$maxLog2FoldChange, maxlfcSE, mixcompdist = "normal", 
                      method = "shrink")
      r$shrunkLFC <- fit$result$PosteriorMean
      r$shrunkSE <- fit$result$PosteriorSD
      r$class <- paste(colnames(effect_matrix)[imax], "V", colnames(effect_matrix)[imin])
    } else {
      warning("Couldn't work out relevant group ordering in LRT")
      r$shrunkLFC <- lfcShrink(dds, res=r, type="ashr", quiet=TRUE)$log2FoldChange
      r$class <- ""
    }
  }  else {
    r$shrunkLFC <- lfcShrink(dds, res=r, type="ashr", quiet=TRUE)$log2FoldChange
    r$class <- ifelse(r$log2FoldChange >0, "Up", "Down")
  }
  ind <- which(r$padj<metadata(r)$alpha)
  r$class[ind] <- paste0(r$class[ind], "*")
  r$class[is.na(r$padj)] <- "Low Count"
  r$class[is.na(r$pvalue)] <- "Outlier"
  r$class[r$baseMean==0] <- "Zero Count"
  mcols(dds)$results <- r
  dds
}

.resNames <- function(names) {
 names[names == "(Intercept)"] <- "Intercept"
 make.names(names)
}

##' Tabulate genelists
##'
##' Get genelist sizes - up, down and classed
##' @title Tabulate the size of the differential lists
##' @param dds 
##' @return 
##' @author Gavin Kelly
##' @export
summarise_results <- function(dds) {
  res <- mcols(dds)$results
  as.data.frame(table(
    Group=sub("\\*$","",res$class),
    Significant=factor(ifelse(grepl("\\*$",res$class), "Significant", "not"), levels=c("Significant","not"))
  )) %>%
    tidyr::spread(Significant, Freq) %>%
    dplyr::mutate(Total=not+Significant) %>%
    dplyr::select(-not) %>%
    dplyr::arrange(desc(Significant/Total))
}    


tidy_significant_dds <- function(dds, res, tidy_fn=NULL, weights=NULL) {
  ind <- grepl("\\*$", res$class)
  mat <- assay(dds, "vst")[ind,,drop=FALSE]
  if (!is.null(weights)) {
    if (is.numeric(weights)) {
      offset <- mat %*%  weights
      mat <- mat - as.vector(offset)
    }
  }
  tidy_dat <- tidy_per_gene(mat, as.data.frame(colData(dds)), tidy_fn)
  return(tidy_dat)
}

tidy_per_gene <- function(mat, pdat,  tidy_fn) {
  if (is.null(tidy_fn)) {
    return(list(mat=mat, pdat=pdat))
  }
  if (inherits(tidy_fn, "fseq")) {
    pdat_long <- dplyr::group_by(cbind(pdat,
                               .value=as.vector(t(mat)),
                               .gene=rep(rownames(mat),each=ncol(mat)),
                               .sample=colnames(mat)),
                         .gene, .add=TRUE)
    summ_long <- dplyr::ungroup(tidy_fn(pdat_long), .gene)
    tidy_pdat <- summ_long[summ_long$.gene==summ_long$.gene[1],]
    tidy_mat <- mat[, tidy_pdat$.sample,drop=FALSE]
    tidy_mat[cbind(summ_long$.gene, summ_long$.sample)] <- summ_long$.value
    tidy_pdat  <- as.data.frame(dplyr::select(tidy_pdat, -.gene, -.value, -.sample))
  } else {
    facts <- c(tidy_fn$by, tidy_fn$rhs, setdiff(tidy_fn$all, unlist(tidy_fn[c("by", "rhs")])))
    ord <- do.call(order, as.list(pdat[,facts, drop=FALSE]))
    return(list(mat=mat[,ord, drop=FALSE], pdat=pdat[ord,facts,drop=FALSE]))
  }
  list(mat=tidy_mat, pdat=tidy_pdat)
}

full_model <- function(mdlList) {
  rhs <- lapply(mdlList, function(mdl) deparse(mdl$design[[2]]))
  fml <- stats::update(as.formula(paste("~", paste(rhs, collapse=" + "))), ~ . )
}


retrieve_contrast <- function (object, expanded = FALSE, listValues=c(1,-1)) {
  comparison <- metadata(object)$comparison
  resNames <- resultsNames(object)
  resReady <- FALSE
  if (is.character(comparison)) {
    if (length(comparison)==1) {
      contrast <- ifelse(resNames==comparison, 1, 0)
      resReady <- TRUE
    } else {
      contrastFactor <- comparison[1]
      contrastNumLevel <- comparison[2]
      contrastDenomLevel <- comparison[3]
      contrastBaseLevel <- levels(colData(object)[, contrastFactor])[1]
      hasIntercept <- attr(terms(design(object)), "intercept") == 1
      firstVar <- contrastFactor == all.vars(design(object))[1]
      noInterceptPullCoef <- !hasIntercept & !firstVar & (contrastBaseLevel %in% 
                                                           c(contrastNumLevel, contrastDenomLevel))
      if (!expanded & (hasIntercept | noInterceptPullCoef)) {
        contrastNumColumn <- make.names(paste0(contrastFactor, "_", contrastNumLevel, "_vs_", contrastBaseLevel))
        contrastDenomColumn <- make.names(paste0(contrastFactor, "_", contrastDenomLevel, "_vs_", contrastBaseLevel))
        if (contrastDenomLevel == contrastBaseLevel) {
          name <- if (!noInterceptPullCoef) {
            make.names(paste0(contrastFactor, "_", contrastNumLevel, "_vs_", contrastDenomLevel))
          }
          else {
            make.names(paste0(contrastFactor, contrastNumLevel))
          }
          contrast <- ifelse(resNames==name, 1,0)
          resReady <- TRUE
        }
        else if (contrastNumLevel == contrastBaseLevel) {
          swapName <- if (!noInterceptPullCoef) {
            make.names(paste0(contrastFactor, "_", contrastDenomLevel, 
                              "_vs_", contrastNumLevel))
          }
          else {
            make.names(paste0(contrastFactor, contrastDenomLevel))
          }
          contrast <- ifelse(resNames==swapName, -1, 0)
          resReady <- TRUE
        }
      }
      else {
        contrastNumColumn <- make.names(paste0(contrastFactor, contrastNumLevel))
        contrastDenomColumn <- make.names(paste0(contrastFactor, contrastDenomLevel))
      }
    }
  }
  if (!resReady) {
    if (is.numeric(comparison)) {
      contrast <- comparison
    }
    else if (is.list(comparison)) {
      contrastNumeric <- rep(0, length(resNames))
      contrastNumeric[resNames %in% comparison[[1]]] <- listValues[1]
      contrastNumeric[resNames %in% comparison[[2]]] <- listValues[2]
      contrast <- contrastNumeric
    }
    else if (is.character(comparison)) {
      contrastNumeric <- rep(0, length(resNames))
      contrastNumeric[resNames == contrastNumColumn] <- 1
      contrastNumeric[resNames == contrastDenomColumn] <- -1
      contrast <- contrastNumeric
    }
  }
  contrast
}

##' Change a factor's reference level
##'
##' Rather than change the order of the levels, this changes the way
##' the factor is parametrised, so that the levels are in the natural order
##' but the coefficients can reflect experimental design considerations
##' @title Rebase a factor's level
##' @param x A factor to be rebased 
##' @param lev The level of the factor that is to be regarded as the 'control' to which all others will be compared
##' @return A factor with a new contrast attribute
##' @author Gavin Kelly
##' @export
rebase <- function(x, lev) {
  i <- which(levels(x)==lev)
  if (length(i)==0) {
    stop(lev, " is not a level of your factor ", deparse(substitute(x)))
  }
  contrasts(x) <- contr.treatment(nlevels(x), i)
  x
 }
