

all_identical <- function(obj)  {all(sapply(obj, function(x) identical(x, obj[[1]])))}

#' Descend a result-tree
#'
#' DESdemonA creates a Dataset > Model > Comparison nested list,
#' with everything stored at the third level of this list. results_apply
#' provides a way to traverse this list, applying a function to each dds object,
#' and reporting the results of that function aggregated to a particular level.
#'
#' For instance, the counts are universal to every model within a dataset, and
#' also to every comparison done within each model in that dataset.  But the
#' actual fitted coeffecients are specific to each comparison.
#' 
#' @title Apply a function to a results list
#' @param object The list of results returned by get_result
#' @param dataset The name of the dataset you want to restrict the traversal to - optional, and if omitted all datasets will be traversed.
#' @param model The name of the model you want to restrict the traversal to - optional, and if omitted all models will be traversed.
#' @param comparison The name of the comparison you want to restrict the traversal to - optional, and if omitted all comparisons will be traversed.
#' @param fn The function that should be applied to each dds object encountered
#' @param depth The depth to which traversal should be carried out.
#' @return A (nested) list of values that represent the return value of fn on the corresponding dds
#' @author Gavin Kelly
#' @export
results_apply <- function(object, dataset, model, comparison, fn, depth="comparison") {
  if (!missing(dataset)) {
    if (!missing(model)) {
      if (!missing(comparison)) {
        return(fn(object[[dataset]][[model]][[comparison]]))
      } else {
        ret <- lapply(
          names(object[[dataset]][[model]]),
          function(comp) DESdemonA::results_apply(object, dataset, model, comp, fn, depth)
        )
        if (depth != "comparison") {
          if (all_identical(ret)) {
            ret <- ret[[1]]
          } else {
            stop("Summarising at the level of ", depth, " but some comparisons are different")
          }
        }
        return(ret)
      }
    } else {
      ret <- lapply(
        names(object[[dataset]]),
        function(mdl)  DESdemonA::results_apply(object, dataset, mdl, fn=fn, depth=depth)
      )
      if (!(depth %in% c("comparison", "model"))) {
        if (all_identical(ret)) {
          ret <- ret[[1]]
        } else {
          stop("Summarising at the level of ", depth, " but some models are different")
        }
      }
      return(ret)
    }
  } else {
      ret <- lapply(
        names(object),
        function(dset)  DESdemonA::results_apply(object, dset, fn=fn, depth=depth)
      )
      return(ret)
  }
}

#' Return all results contained in a result list
#'
#' Traverse an object returned by get_result, and extract the DESeq2 results stored at each node.
#' @title Return all results contained in a result list
#' @param object The list of results returned by get_result
#' @param dataset The name of the dataset you want to restrict the traversal to - optional, and if omitted all datasets will be traversed.
#' @param model The name of the model you want to restrict the traversal to - optional, and if omitted all models will be traversed.
#' @param comparison The name of the comparison you want to restrict the traversal to - optional, and if omitted all comparisons will be traversed.
#' @return 
#' @author Gavin Kelly
#' @export
extract_results <- function(object, dataset, model, comparison ) {
  results_apply(object, dataset, model, comparison, fn=function(dds) mcols(dds)$results, depth="comparison")
}


extract_model <- function(object, dataset) {
  results_apply(object, dataset, fn=function(dds) mcols(dds)$results, depth="model")
}

#' Return an assay from a result-list
#'
#' Traverse an object returned by get_result, and extract a specific assay
#' @title Return an assay from a result-list
#' @param dataset The name of the dataset you want to extract the assay from  - optional, and if omitted all datasets will be reported.
#' @param object The list of results returned by get_result
#' @param assay The name of the assay you want to extract from each dataset
#' @return 
#' @author Gavin Kelly
#' @export
extract_assay <- function(object,dataset,  assay) {
  results_apply(object, dataset, fn=function(dds) assay(dds, assay), depth="dataset")
}  


#' Return PCA-projected data from a result-list
#'
#' Traverse an object returned by get_result, and extract the expression projected onto PC space.
#' The models, and comparisons within those, should all have the same projections, so the results
#' are aggregated to the level of 'dataset'
#' @title Return PCA-projected data from a result-list
#' @param dataset The name of the dataset you want to extract the PC coordinates  - optional, and if omitted all datasets will be reported.
#' @param object The list of results returned by get_result
#' @return A DataFrame of PC coordinates for each dataset specified
#' @author Gavin Kelly
#' @export
extract_pca_x <- function(object, dataset) {
  results_apply(object, dataset, fn=function(dds) colData(dds)$.PCA, depth="dataset")
}

#' Return colData from a result-list
#'
#' Traverse an object returned by get_result, and extract the colData. As the
#' comparisons and models within a dataset are all based off the same colData,
#' the return value is aggregated to the level of dataset.
#' @title Return PCA-projected data from a result-list
#' @param dataset The name of the dataset you want to extract the PC coordinates  - optional, and if omitted all datasets will be reported.
#' @param object The list of results returned by get_result
#' @return A DataFrame of PC coordinates for each dataset specified
#' @author Gavin Kelly
#' @export
extract_colData <- function(object, dataset) {
  results_apply(object, dataset, fn=function(dds) colData(dds), depth="dataset")
}



#' Wrap each dataset of a  DESdemonA fitted object 
#'
#' If you want to loop through datasets that have been fitted by
#' DESdemonA, then use this function as a wrapped so that reports
#' can generate formatted output, and transformations can be applied.
#'
#' This should be the first process in a pipeline, and be followed by
#' either a per_model or per_comparison call: the primary function .f
#' will be applied to the output.  For each dataset in the dds
#' object, the 'before' function will be called before passing things
#' on (and so can be used to print a header, for example), and the
#' 'after' function will be called after (to produce captions).  All
#' functions have access to a variable '.dataset' that will take the
#' name of the current dataset
#' 
#' @title Wrap datasets from DESdemonA
#' @param .dds A three-level list in the standard DESdemonA hierarchy dataset > model > comparison
#' @param .f The function that will be called on the output of the per_model or per_comparison.
#' @param before A function that will be invoked before .f
#' @param after A function that will be invoked after .f
#' @return The return value of .f, as applied to the list of child outputs
#' @author Gavin Kelly 
#' @export
per_dataset <- function(.dds, .f=identity, before=as.null, after=as.null, deepen=FALSE, ...) {
  .data <- list(dds=.dds,
               dataset_f=.f, dataset_before=before, dataset_after=after,
               model_f=identity, model_before=as.null, model_after=as.null
               )
  class(.data) <- "DesIterator"
  if (class(.dds[[1]])=="DESeqDataSet") {
    .data$dds <- lapply(.dds, function(x) list(model=list(comparison=x)))
    .data$dataset_f=identity
    ret <- per_comparison(.data, .f, ...)
    if (deepen) {
      ret
    } else {
      lapply(ret, function(x) x$model$comparison)
    }
  } else {
    .data
  }
}

#' Wrap each model of a DESdemonA fitted object
#'
#' If you want to loop through models that have been fitted by
#' DESdemonA, then use this function as a wrapped so that reports
#' can generate formatted output, and transformations can be applied.
#'
#' This can be the first process in a pipeline, or follow a
#' 'per_dataset' process and must be followed by a per_comparison
#' call: the primary function .f will be applied to the output.  For
#' each model in the dds object, the 'before' function will be called
#' before passing things on (and so can be used to print a header, for
#' example), and the 'after' function will be called after (to produce
#' captions).  All functions have access to a variables '.dataset' and
#' '.model' that will take the name of the current dataset and model
#' 
#' @title Wrap models from DESdemonA
#' @param x Either a DESdemonA three-level list, or a previous 'per_dataset' call
#' @param .f The function that will be called on the output of per_comparison.
#' @param before A function that will be invoked before .f
#' @param after A function that will be invoked after .f
#' @return The value of .f applied to the list of child outputs from per_comparison
#' @author Gavin Kelly
#' @export
per_model <- function(x, .f, before, after, deepen) {
  UseMethod("per_model") 
}
per_model.list <- function(.dds, .f = identity, before=as.null, after=as.null, deepen=FALSE) {
  .data <- list(dds=.dds,
               dataset_f=identity, dataset_before=as.null, dataset_after=as.null
               )
  class(.data) <- "DesIterator"
  per_model(.data, .f, before, after, deepen)
}

per_model.DesIterator <- function(.data, .f = identity, before=as.null, after=as.null, deepen=FALSE) {
  .data$model_f <- .f
  .data$model_before=before
  .data$model_after=after
  if (class(.data$dds[[1]][[1]])=="DESeqDataSet") {
    .data$dds <- lapply(.data$dds, function(x) lapply(x, function(y) list(comparison=y)))
    .data$model_f=identity
    ret <- per_comparison(.data, .f)
    if (deepen) {
      ret
    } else {
      lapply(ret, function(x) lapply(x, function(y) y$comparison))
    }
  } else {
    .data
  }
}


#' Wrap each comparison of a DESdemonA fitted object
#'
#' If you want to loop through comparisons that have been fitted by
#' DESdemonA, then use this function as a wrapped so that reports
#' can generate formatted output, and transformations can be applied.
#'
#' This can be called on its own, or at the end of a pipeline of
#' per_model and per_comparison: the primary function .f will be
#' applied to each dds comparison object.  For each comparison in the
#' dds object, the 'before' function will be called before passing
#' things on (and so can be used to print a header, for example), and
#' the 'after' function will be called after (to produce captions).
#' All functions have access to a variables '.dataset', '.model' and
#' '.comparison' that will take the name of the current dataset,
#' model and comparison
#' 
#' @title Wrap models from DESdemonA
#' @param x Either a DESdemonA three-level list, or a previous 'per_dataset' or 'per_model' call
#' @param .f The function that will be called on the DESeq object.
#' @param before A function that will be invoked before .f
#' @param after A function that will be invoked after .f
#' @return A three-level list of results of .f applied to each comparison
#' @author Gavin Kelly
#' @export

per_comparison <- function(x, .f, before, after, ...) {
  UseMethod("per_comparison")
}
per_comparison.list <- function(.dds, .f=identity, before=as.null, after=as.null, ...) {
  .data <- list(dds=.dds,
               dataset_f=identity, dataset_before=as.null, dataset_after=as.null,
               model_f=identity, model_before=as.null, model_after=as.null
               )
  class(.data) <- "DesIterator"
  per_comparison.DesIterator(.data, .f=.f, before=before, after=after, ...)
}
per_comparison.DesIterator <- function(.data, .f=identity, before=as.null, after=as.null, ...) {
  dataset_ret <- list()
  for (.dataset in names(.data$dds)) {
    assign(".dataset", .dataset, envir=environment(.data$dataset_before))
    assign(".dataset", .dataset, envir=environment(.data$dataset_after))
    rlang::as_function(.data$dataset_before)()
    model_ret <- list()
    for (.model in names(.data$dds[[.dataset]])) {
      assign(".dataset", .dataset, envir=environment(.data$model_before))
      assign(".dataset", .dataset, envir=environment(.data$model_after))
      assign(".model", .model, envir=environment(.data$model_before))
      assign(".model", .model, envir=environment(.data$model_after))
      (rlang::as_function(.data$model_before))()
      comparison_ret <- list()
      for (.comparison in names(.data$dds[[.dataset]][[.model]])) {
        .dds <- .data$dds[[.dataset]][[.model]][[.comparison]]
        assign(".dataset", .dataset, envir=environment(before))
        assign(".dataset", .dataset, envir=environment(after))
        assign(".model", .model, envir=environment(before))
        assign(".model", .model, envir=environment(after))
        assign(".comparison", .comparison, envir=environment(before))
        assign(".comparison", .comparison, envir=environment(after))
        for (i in setdiff(names(.data), "dds")) {
          assign(".dds", .dds, envir=environment(.data[[i]]))
        }
        (rlang::as_function(before))()
        comparison_ret[[.comparison]] <- (rlang::as_function(.f))(.dds, ...)
        (rlang::as_function(after))()
      }
      model_ret[[.model]] <- (rlang::as_function(.data$model_f))(comparison_ret)
      (rlang::as_function(.data$model_after))()
    }
    dataset_ret[[.dataset]] <- (rlang::as_function(.data$dataset_f))(model_ret)
    (rlang::as_function(.data$dataset_after))()
  }
  dataset_ret
}

