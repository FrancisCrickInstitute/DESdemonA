nested_merge <- function(l1, l2) {
  if (length(l1)==0) {
    return(l2)
  }
  if (length(l2)==0) {
    return(l1)
  }
  if (class(l1[[1]])=="list") {
    is_common <- as.list(intersect(names(l1), names(l2)))
    names(is_common) <- unlist(is_common)
    merged <- lapply(is_common, function(l) {nestedMerge(l1[[l]], l2[[l]])})
    c(merged, l1[setdiff(names(l1), names(l2))], l2[setdiff(names(l2), names(l1))])
  } else {
    c(l1, l2)
  }
}

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
