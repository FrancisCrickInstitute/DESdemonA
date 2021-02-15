##' Make captioning function that generates hyperlinks
##'
##' captioner returns a function that can be used in `fig.cap`.  Call
##' the resulting function with a string containing the caption text,
##' after any plot, to store a link to the pdf version of the
##' plot. Call the function without any arguments (e.g. in the
##' `fig.cap` argument of a chunk) to invoke the captioining mechanism
##' in the markdown.
##' @title Caption hyperlinking
##' @return
##' @author Gavin Kelly
captioner <- function() {
  local({
    captions <- character()
    function(caption) {
      if (missing(caption)) {
        ret <- captions
        captions <<- character()
        return(ret)
      } else {
        if (isTRUE(getOption('knitr.in.progress'))) {
          link_caption <- paste0("[", caption, "](", knitr::fig_path("pdf", number=length(captions)+1), ")")
        } else {
          link_caption <- caption
        }
        captions <<- c(captions, link_caption)
      }
    }
  })
}

##' Link GT tables to a csv file
##'
##' To be used in a `GT` pipeline, it will store the underlying table
##' data in a csv file under the given name, and insert a link in the
##' table's caption that points to the csv file.
##' @title Link GT tables to a csv file
##' @param data The GT object with a caption set
##' @param name The filename of the csv
##' @return Th GT object (invisibly)
##' @author Gavin Kelly
tab_link_caption <- function(data,name) {
  if (missing(name)) {
    heading <- gt:::dt_heading_get(data)
    name <- paste(heading$title, heading$subtitle)
  }
  caption <- gt:::dt_options_get_value(data = data, option = "table_caption")
  fname <- knitr::fig_path("csv", number=name)
  write.csv(gt:::dt_data_get(data), file=fname)
  gt:::dt_options_set_value(
    data,
    "table_caption",
    paste0("[", caption, "](", fname, ")")
  )
  invisible(data)
}

##' Deprecated
##'
##' To be replaced by tab_link_caption
##' @title deprecated
##' @param lbl the label
##' @return 
##' @author Gavin Kelly
klabel <- function(lbl) paste0(knitr::opts_current$get('label'), "-", lbl)

