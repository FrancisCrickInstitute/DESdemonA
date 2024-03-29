% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.r
\name{differential_heatmap}
\alias{differential_heatmap}
\title{Heatmaps of differential genes}
\usage{
differential_heatmap(ddsList, tidy_fn = NULL, param, caption)
}
\arguments{
\item{ddsList}{DESdemonA-generated list of \link[DESeq2:DESeqDataSet]{DESeq2::DESeqDataSet}s}

\item{tidy_fn}{A dplyr pipeline to transform the expression values
in a convenient manner. Use this for example to 'normalise' the
data so that all values are relative to a particular condition.}

\item{caption}{The function to use to create captions}
}
\value{

}
\description{
Heatmaps of expression for differential genes
}
\details{
For each set of results, plot the heatmap of counts, limited to
differential genes
}
\author{
Gavin Kelly
}
