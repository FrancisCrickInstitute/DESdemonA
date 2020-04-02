captioner <- function(filetype="pdf") {
  local({
    captions <- character()
    function(caption, dims=FALSE) {
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

klabel <- function(lbl) paste0(knitr::opts_current$get('label'), "-", lbl)

rmd_link_data <- function(fname) {
  paste0("\n\n<a href=\"", fname, "\"><i class=\"fa fa-file-download\"></i>Download data</a>")
}
