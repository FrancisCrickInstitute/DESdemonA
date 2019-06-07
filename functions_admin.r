#' Supply pretty-printed git stats
#'
#'
#'

git_stats <- function() {
  vers <- try(suppressWarnings(system2("git", "log -1 --pretty=format:%h", stdout=TRUE, stderr=FALSE)), silent=TRUE)
  gitTag <- try(suppressWarnings(system2("git", "describe --tags --exact", stdout=TRUE, stderr=FALSE)), silent=TRUE)
  if (class(gitTag)=="try-error" || (!is.null(attr(gitTag, "status")) && attr(gitTag, "status")!=0)) {
    gitTag <- NA
  }
  if (length(vers)==0 || class(vers)=="try-error") {
    vers <- "uncontrolled"
    gitTag <- NA
  } else {
    if (any(grepl("^ M",system2("git", "status --porcelain", stdout=TRUE, stderr=FALSE)))) {
      vers <- sprintf("%s-M", vers)
      gitTag <- NA
    }
  }
  list(version=vers, tag=gitTag)
}

#' Open a device with a trackable name
#'
#' Can be used with 'pdf' 'write.table' etc
#'
#' @param dev contains the function which will be called (which itself will get its 'file' param filled in automatically)
#' @param file string will be used as a starting point for the filename.  If it doesn't contain an extension, will predict from 'dev'
#' @param ignore.chunk if set,  will over-rule the use of knitr override 
#' @return the filename actually used
vDevice <- function(dev=c, file="plot", ..., dir="results", use_chunk_label=FALSE) {
  if (!grepl("\\.", file)) {
    devstr <- as.character(substitute(dev)) # e.g. 'pdf', 'write.txt'
    devstr <- sub(".*\\.", "", devstr) # e.g. 'pdf', 'txt'
    file <- paste0(file, ".", devstr)
  }
  if (isTRUE(getOption('knitr.in.progress'))  &  use_chunk_label) {
    file <- paste0(knitr::fig_path(sub(".*\\.", "", file)) )
    dev(..., file=file)
    return(file)
  }
  g <- git_stats()
  gv <- g$version
  if (is.na(g$tag)) {
    dName <- file.path(dir, gv)
  } else {
    dName <- file.path(dir, g$tag)
  }
  if (!dir.exists(dName)) {
    dir.create(dName)
  }
  file=file.path(dName,sub("([^.]*)(.*)",paste0("\\1_", gv, "\\2"), file))
  dev(..., file= file)
  file
}


rmd_link_data <- function(fname) {
  paste0("\n\n<a href=\"", fname, "\"><i class=\"fa fa-file-download\"></i>Download data</a>")
}
