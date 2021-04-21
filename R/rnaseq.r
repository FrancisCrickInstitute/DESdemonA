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



##' Subsample
##'
##' Think this should probably be deprecated
##' @title subsample
##' @param dds Not Sure
##' @param subs Not Sure
##' @return Not Sure
##' @author Gavin Kelly
##' @export
subsample <- function(dds, subs) {
  if (is_formula(subs)) {
    grps <- Reduce(interaction, colData(dds)[all.vars(subs)])
    grps <- names(table(grps))[table(grps)!=0]
    out <- lapply(setNames(grps,grps),  function(x) DESdemonA::subsample(dds, grps==x))
  } else {
    out <- dds[,subs]
    colData(out) <- droplevels(colData(out))
    out
  }
  out
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
  lapply(spec$sample_sets, function(set) {
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
      colData(obj) <- S4Vectors::DataFrame(eval(tr))
    }
    obj
  })
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
    function(mdl)  {
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
      out
    }
  )
  model_comp[sapply(model_comp, length)!=0] 
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
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
    if ("drop_unsupported_combinations" %in% names(mdl)) {
      mdl$dropped <- is.na(coef(fit))
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


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Expand multiple comparisons into their contrasts
##' @param dds The original DESeq2 object containing all samples
##' @param spec 
##' @param ... 
##' @return 
##' @author Gavin Kelly
##' @export
emcontrasts <- function(dds, spec, extra=NULL) {
  ## df <- as.data.frame(colData(dds))
  ## df$.x <- counts(dds, norm=TRUE)[1,]
  ## if (is_formula(design(dds))) {
  ##   fml <- as.formula(paste0(".x ~ ", as.character(design(dds)[2])))
  ##   fit <- lm(fml, data=df)
  ## } else { # shouldn't happen, but
  ##   warning("Seem to be doing EM on a matrix - not sure that's great")
  ##   mm <- design(dds)
  ##   colnames(mm) <- DESdemonA:::.resNames(colnames(mm))
  ##   fit <- lm(df$.x ~ . -1, data.frame(mm))
  ##   ddsNames <- match(resultsNames(dds), names(coef(fit)))
  ## }
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
  contr_frame <- contr_frame[ind_est,1:(which(names(contr_frame)=="estimate")-1), drop=FALSE]
  if (!is.na(keep[1])) {
    contr_frame <- subset(contr_frame, contrast %in% keep)
  }
  contr_mat <- emfit$contrast@linfct[ind_est, !mdl$dropped, drop=FALSE]
  colnames(contr_mat) <- DESdemonA:::.resNames(colnames(contr_mat))
  contr <- lapply(seq_len(nrow(contr_frame)), function(i) contr_mat[i,,drop=TRUE])
  names(contr) <- do.call(paste, c(contr_frame,sep= "|"))
  contr
}


## Do likelihood ratio test, and classify in order of effect size
## Doesn't work for interactions, obviously
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Fit LRT
##' @param dds The original DESeq2 object containing all samples
##' @param reduced 
##' @param ... 
##' @return 
##' @author Gavin Kelly
fitLRT <- function(dds, reduced, ...) {
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
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Generate the results for a model and comparison
##' @param dds The original DESeq2 object containing all samples
##' @param mcols 
##' @param filterFun 
##' @param lfcThreshold 
##' @param ... 
##' @return 
##' @author Gavin Kelly
##' @export
get_result <- function(dds, mcols=c("symbol", "entrez"), filterFun=IHW::ihw, lfcThreshold=0,  ...) {
  if (is.null(filterFun)) filterFun <- rlang::missing_arg()
  comp <- metadata(dds)$comparison
  if (!is_formula(comp)) {
    if (is.character(comp) && length(comp)==1) { #  it's a name
      r <- DESeq2::results(dds, filterFun=filterFun, lfcThreshold=lfcThreshold, name=metadata(dds)$comparison, ...)
    } else { # it's a contrast
      if (is.list(comp) && "listValues" %in% names(comp)) {
        r <- DESeq2::results(dds, filterFun=filterFun, lfcThreshold=lfcThreshold, contrast=metadata(dds)$comparison[names(comp) != "listValues"], listValues=comp$listValues, ...)
      } else {
        r <- DESeq2::results(dds, filterFun=filterFun, lfcThreshold=lfcThreshold, contrast=metadata(dds)$comparison, ...)
      }
    }
  } else { # it's LRT
    r <- results(dds, filterFun=filterFun, ...)
  }
  r[mcols] <- mcols(dds)[mcols]
  if ("LRTPvalue" %in% names(mcols(dds))) {
    r$class <- mcols(dds)$class
    r$class[is.na(r$padj) | is.na(r$pvalue) | r$baseMean==0] <- NA
    term <-  metadata(dds)$LRTterms
    # take the biggest fold-change vs baseline, for MA and reporting?
    if (all(term %in% names(mcols(dds)))) {
      effect_matrix <- cbind(I=rep(0, nrow(dds)),as.matrix(mcols(dds)[,term,drop=FALSE]))
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
  r$class[is.na(r$pvalue)] <- "Outlier"
  r$class[is.na(r$padj)] <- "Low Count"
  r$class[r$baseMean==0] <- "Zero Count"
  mcols(dds)$results <- r
  dds
}

.resNames <- function(names) {
 names[names == "(Intercept)"] <- "Intercept"
 make.names(names)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Enrichment analysis
##' @param ddsList 
##' @param fun 
##' @param showCategory 
##' @param max_width 
##' @return 
##' @author Gavin Kelly
##' @export
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
  mat <- assay(vst(dds))[ind,,drop=FALSE]
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
  }
  ## pdat$.value <- mat[1,]
  ## tidy_mat <- apply(mat, 1, function(x) {pdat$.value=x; tidy_fn(pdat)$.value})
  ## tidy_pdat <- tidy_fn(pdat) %>% dplyr::select(-.value)
  list(mat=tidy_mat, pdat=tidy_pdat)
}

