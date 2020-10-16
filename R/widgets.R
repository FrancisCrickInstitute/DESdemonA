#' Interactive PCA plot
#'
#' Generate an SVG scatterplot for rmarkdown reports where the user
#' can choose which dimensions to plot, and the colouring scheme.
#'
#' @import htmlwidgets
#'
#' @export
d3_iScatter <- function(data, xy=grepl("^PC[0-9]+$", names(data)), aes=!xy, sample_names=row.names(data), width = NULL, height = NULL, elementId = NULL) {

  x = list( # settings that will be forwarded
    df=as.data.frame(data),
    xy=xy,
    aes=aes,
    sample_names=sample_names
  )

  # create widget
  htmlwidgets::createWidget(
    name = 'd3_iScatter',
    x,
    width = width,
    height = height,
    package = getPackageName(),
    elementId = elementId
  )
}

#' Shiny bindings for d3_iScatter
#'
#' Output and render functions for using d3_iScatter within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a d3_iScatter
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name d3_iScatter-shiny
#'
#' @export
d3_iScatterOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'd3_iScatter', width, height, package = getPackageName())
}

#' @rdname d3_iScatter-shiny
#' @export
renderd3_iScatter <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, d3_iScatterOutput, env, quoted = TRUE)
}

d3_iScatter_html <- function(id, style, class, ...){
  htmltools::htmlTemplate(system.file("htmlwidgets/d3_iScatter_template.html",package=getPackageName()), id=id, style=style, class=class)
}




#' Interactive metric plot
#'
#' Generate an SVG scatterplot for rmarkdown reports where the user
#' can choose which dimensions to plot, and the colouring scheme.
#'
#' @import htmlwidgets
#'
#' @export
d3_iMetric <- function(data, width = NULL, height = NULL, elementId = NULL) {

  x = list( # settings that will be forwarded
    df=data
  )

  # create widget
  htmlwidgets::createWidget(
    name = 'd3_iMetric',
    x,
    width = width,
    height = height,
    package = getPackageName(),
    elementId = elementId
  )
}

#' Shiny bindings for d3_iMetric
#'
#' Output and render functions for using d3_iMetric within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a d3_iMetric
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name d3_iMetric-shiny
#'
#' @export
d3_iMetricOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'd3_iMetric', width, height, package = getPackageName())
}

#' @rdname d3_iMetric-shiny
#' @export
renderd3_iMetric <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, d3_iMetricOutput, env, quoted = TRUE)
}

d3_iMetric_html <- function(id, style, class, ...){
  htmltools::htmlTemplate(system.file("htmlwidgets/d3_iMetric_template.html",package=getPackageName()), id=id, style=style, class=class)
}
