##' rbind nested lists of data-frames
##'
##' Walk through a DESdemonA hierarchy of summaries, and join them
##' into one big summary table
##' @title rbind nested lists of data-frames
##' @param x the list of summaries
##' @param levels
##' @return The big table
##' @author Gavin Kelly
##' @export
rbind_summary <- function(x, levels=c("Dataset", "Design", "Comparison")) {
  if (class(x[[1]])=="list") {
    do.call(rbind, mapply(function(df, name) {
      extended <- cbind(name, df)
      names(extended)[1] <- levels[1]
      return(extended)
    }, lapply(x, rbind_summary, levels[-1]), names(x), SIMPLIFY=FALSE))
  } else {
    do.call(rbind, mapply(function(df, name) {
      extended <- cbind(name, df)
      names(extended)[1] <- levels[1]
      return(extended)
    }
   ,x, names(x), SIMPLIFY=FALSE))
  }
}


## classify_terms <- function(fml) {
##   has_random <- any(c("|", "||") %in% all.names(fml))
##   if (!has_random) {
##     return(list(fixed=all.vars(fml),groups=list() ))
##   }
##   rterms <- lme4:::barnames(lme4::findbars(fml))
##   rterms <- sub(":.*", "", rterms)
##   fvars <- setdiff(all.vars(lme4::subbars(fml)), rterms)
##   list(fixed=fvars, groups=rterms)
## }

classify_terms <- function(fml) {
  rhs <- as.character(fml)[length(fml)]
  spl <- strsplit(rhs, "|", fixed=TRUE)[[1]]
  rv <- all.vars(stats::reformulate(spl[1]))
  list(
    fixed=rv,
    groups=setdiff(all.vars(stats::reformulate(rhs)), rv)
  )
}

  
  
