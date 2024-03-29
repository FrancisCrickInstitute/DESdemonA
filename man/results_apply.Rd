% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.r
\name{results_apply}
\alias{results_apply}
\title{Apply a function to a results list}
\usage{
results_apply(object, dataset, model, comparison, fn, depth = "comparison")
}
\arguments{
\item{object}{The list of results returned by get_result}

\item{dataset}{The name of the dataset you want to restrict the traversal to - optional, and if omitted all datasets will be traversed.}

\item{model}{The name of the model you want to restrict the traversal to - optional, and if omitted all models will be traversed.}

\item{comparison}{The name of the comparison you want to restrict the traversal to - optional, and if omitted all comparisons will be traversed.}

\item{fn}{The function that should be applied to each dds object encountered}

\item{depth}{The depth to which traversal should be carried out.}
}
\value{
A (nested) list of values that represent the return value of fn on the corresponding dds
}
\description{
Descend a result-tree
}
\details{
DESdemonA creates a Dataset > Model > Comparison nested list,
with everything stored at the third level of this list. results_apply
provides a way to traverse this list, applying a function to each dds object,
and reporting the results of that function aggregated to a particular level.

For instance, the counts are universal to every model within a dataset, and
also to every comparison done within each model in that dataset.  But the
actual fitted coeffecients are specific to each comparison.
}
\author{
Gavin Kelly
}
