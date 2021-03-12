

all_identical <- function(obj)  {all(sapply(obj, function(x) identical(x, obj[[1]])))}

##' Descend a result-tree
##'
##' DESdemonA creates a Dataset > Model > Comparison nested list,
##' with everything stored at the third level of this list. results_apply
##' provides a way to traverse this list, applying a function to each dds object,
##' and reporting the results of that function aggregated to a particular level.
##'
##' For instance, the counts are universal to every model within a dataset, and
##' also to every comparison done within each model in that dataset.  But the
##' actual fitted coeffecients are specific to each comparison.
##' 
##' @title Apply a function to a results list
##' @param object The list of results returned by get_result
##' @param dataset The name of the dataset you want to restrict the traversal to - optional, and if omitted all datasets will be traversed.
##' @param model The name of the model you want to restrict the traversal to - optional, and if omitted all models will be traversed.
##' @param comparison The name of the comparison you want to restrict the traversal to - optional, and if omitted all comparisons will be traversed.
##' @param fn The function that should be applied to each dds object encountered
##' @param depth The depth to which traversal should be carried out.
##' @return A (nested) list of values that represent the return value of fn on the corresponding dds
##' @author Gavin Kelly
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

##' Return all results contained in a result list
##'
##' Traverse an object returned by get_result, and extract the DESeq2 results stored at each node.
##' @title Return all results contained in a result list
##' @param object The list of results returned by get_result
##' @param dataset The name of the dataset you want to restrict the traversal to - optional, and if omitted all datasets will be traversed.
##' @param model The name of the model you want to restrict the traversal to - optional, and if omitted all models will be traversed.
##' @param comparison The name of the comparison you want to restrict the traversal to - optional, and if omitted all comparisons will be traversed.
##' @return 
##' @author Gavin Kelly
#' @export
extract_results <- function(object, dataset, model, comparison ) {
  results_apply(object, dataset, model, comparison, fn=function(dds) mcols(dds)$results, depth="comparison")
}


extract_model <- function(object, dataset) {
  results_apply(object, dataset, fn=function(dds) mcols(dds)$results, depth="model")
}

##' Return an assay from a result-list
##'
##' Traverse an object returned by get_result, and extract a specific assay
##' @title Return an assay from a result-list
##' @param dataset The name of the dataset you want to extract the assay from  - optional, and if omitted all datasets will be reported.
##' @param object The list of results returned by get_result
##' @param assay The name of the assay you want to extract from each dataset
##' @return 
##' @author Gavin Kelly
#' @export
extract_assay <- function(object,dataset,  assay) {
  results_apply(object, dataset, fn=function(dds) assay(dds, assay), depth="dataset")
}  


##' Return PCA-projected data from a result-list
##'
##' Traverse an object returned by get_result, and extract the expression projected onto PC space.
##' The models, and comparisons within those, should all have the same projections, so the results
##' are aggregated to the level of 'dataset'
##' @title Return PCA-projected data from a result-list
##' @param dataset The name of the dataset you want to extract the PC coordinates  - optional, and if omitted all datasets will be reported.
##' @param object The list of results returned by get_result
##' @return A DataFrame of PC coordinates for each dataset specified
##' @author Gavin Kelly
#' @export
extract_pca_x <- function(object, dataset) {
  results_apply(object, dataset, fn=function(dds) colData(dds)$.PCA, depth="dataset")
}

##' Return colData from a result-list
##'
##' Traverse an object returned by get_result, and extract the colData. As the
##' comparisons and models within a dataset are all based off the same colData,
##' the return value is aggregated to the level of dataset.
##' @title Return PCA-projected data from a result-list
##' @param dataset The name of the dataset you want to extract the PC coordinates  - optional, and if omitted all datasets will be reported.
##' @param object The list of results returned by get_result
##' @return A DataFrame of PC coordinates for each dataset specified
##' @author Gavin Kelly
#' @export
extract_colData <- function(object, dataset) {
  results_apply(object, dataset, fn=function(dds) colData(dds), depth="dataset")
}

