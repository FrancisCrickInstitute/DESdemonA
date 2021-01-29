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
}

klabel <- function(lbl) paste0(knitr::opts_current$get('label'), "-", lbl)

rmd_link_data <- function(fname) {
  paste0("\n\n<a href=\"", fname, "\"><i class=\"fa fa-file-download\"></i>Download data</a>")
}
