

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



#' Map over each component of a  DESdemonA object
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
#' @param .dmc A three-level list in the standard DESdemonA hierarchy dataset > model > comparison
#' @param .f The function that will be called on the output of the per_model or per_comparison.
#' @param before A function that will be invoked before .f
#' @param after A function that will be invoked after .f
#' @return The return value of .f, as applied to the list of child outputs
#' @author Gavin Kelly 

#' @export
## map_des <- function(data, f=identity, depth="comparison", ...) {
##   if (depth=="top") {
##     return((rlang::as_function(f))(data, ...))
##   }
##   dataset_ret <- list()
##   for (dataset in names(data)) {
##     if (class(data[[dataset]])=="list" && depth!="dataset") {
##       model_ret <- list()
##       for (model in names(data[[dataset]])) {
##         if (class(data[[dataset]][[model]])=="list" && depth=="comparison") {
##           comparison_ret <- list()
##           for (comparison in names(data[[dataset]][[model]])) {
##             comparison_ret[[comparison]] <- (rlang::as_function(f))(data[[dataset]][[model]][[comparison]], ...)
##           }
##           if (length(comparison_ret)>0) {
##             model_ret[[model]] <-comparison_ret
##           }
##         } else {
##           model_ret[[model]] <- (rlang::as_function(f))(data[[dataset]][[model]], ...)
##         }
##       }
##       if (length(model_ret)>0) {
##         dataset_ret[[dataset]] <- model_ret
##       }
##     } else {
##       dataset_ret[[dataset]] <- (rlang::as_function(f))(data[[dataset]], ...)
##     }
##   }
##   dataset_ret
## }


map_des <- function(data, f, depth="comparison",...) {
  if (class(data[[1]])!="list" || depth=="dataset") {
    return(map(data,f))
  }
  if (class(data[[1]][[1]])!="list" || depth=="model") {
    return(map_depth(data, 2, f))
  }
  map_depth(data, 3, f)
}

trim_map <- function(data) {
  deeper <- !sapply(data, is.null)
  if (!any(deeper)) {
    return(NULL)
  }
  lapply(data[deeper], function(x) {if (class(x)=="list") trim_map(x) else x})
}
  
  
#' @export
dataset_name <- function(dds) {
  metadata(dds)$dmc$dataset
}
#' @export
model_name <- function(dds) {
  metadata(dds)$dmc$model
}
#' @export
comparison_name <- function(dds) {
  metadata(dds)$dmc$comparison
}



#' @export
dmc2frame <- function(dmc) {
  df_rows <- map_des(dmc, function(x) {
    data.frame(dataset=dataset_name(x),
               model=model_name(x),
               comparison=comparison_name(x))
  })
  df <- unlist(unlist(df_rows, recursive=FALSE), recursive=FALSE)
  dds <- unlist(unlist(dmc, recursive=FALSE), recursive=FALSE)[names(df)]
  df <- do.call(rbind, df)
  df$dds <- dds
  df
}

