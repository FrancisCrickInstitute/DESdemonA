## df2colorspace <- function(df) {
##   purrr::imap(df,
##        ~setNames(colorspace::qualitative_hcl(length(unique(.x)), row.names(colorspace::hcl_palettes("qualitative"))[match(.y, names(df))]),
##                  unique(.x)))
## }


df2colorspace <- function(df) {
  pal <- RColorBrewer::brewer.pal(12, "Set3")
  df <- dplyr::mutate_if(as.data.frame(df), is.character, as.factor)
  purrr::map2(df,
              cumsum(c(0,map(df, nlevels)))[1:length(df)],
              function(column, start_level) {
                if (is.factor(column)) {
                  setNames(pal[(seq(start_level, length=nlevels(column)) %% length(pal)) + 1],
                           levels(column))
                } else {
                  circlize::colorRamp2(range(column), c("white","red"))
                }
              }
              )
}
