% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.r
\name{export_biologic}
\alias{export_biologic}
\title{Export Biologic files}
\usage{
export_biologic(result_object, path)
}
\arguments{
\item{result_object}{Path to the rds of the DESdemonA result object}

\item{path}{Where to write the files for Biologic to read}
}
\value{
nothing
}
\description{
Generate text files require for Biologic
}
\details{
The visualisation part of 'Biologic' requires a set of text files
which represent the differential results and the models and
contrasts that were used to derive the results.  This function
generates them from a previous run of DESdemonA.
}
\author{
Gavin Kelly
}
