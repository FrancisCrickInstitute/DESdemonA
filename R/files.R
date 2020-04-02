##' Access extdata files whether they're in inst or not
##'
##' To build a package files should live in inst/extdata,
##' whereas to refer to them in the package, they should have the
##' inst in their path. devtools::system.file deals with this.
##' This function is a helper for the devtools-modified system.file but
##' with the added benefit that you can put whatever you want into "extdata" and
##' whatever you ever use will also get linked to from "inst/extdata" so that
##' extra files in extdata that aren't loaded by the package can easily be excluded.
##' 
##' @title Intelligent use of extdata files
##' @param ... The path elements, without any need for the 'extdata' part
##' @return The file-path, and a side-effect of creating a link in inst/extdata if necessary
##' @author Gavin Kelly
extdata <- function(...) {
  fname <- system.file("extdata", ..., package=getPackageName())
  if (!grepl(paste0(normalizePath(file.path("inst", "extdata", ...)), "$"), fname)) {
    file.symlink(fname, sub("/extdata/", "/inst/extdata/", fname))
  }
  fname
}

