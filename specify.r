list(
  sample_sets =  list(
    all = TRUE,
    no_cancer = fraction != "tumour"
  ),
  
  models = list(
    "Batch-corrected" = list(
      design = ~ fraction + pool,
      comparison = mult_comp(pairwise ~ fraction)
    ),
    ## "By Hand" = list(
    ##   design = ~ fraction + pool, 
    ##   comparison = list(                          # can be a list
    ##      "f1 v f2" = list("fraction", "f2","f1"), # of DESeq2::results(contrast=e.g.
    ##      "fraction" = ~ pool,                     # or LRT full=design, reduced=this
    ##      multcomp(pairwise ~ fraction | pool)     # and this will get expanded out to all pairwise
    ## ),
    ## "Complex" = list(
    ##   design = ~ fraction * pool,                               # No replicates, but just imagine if there were
    ##   comparison = list(
    ##    mult_comp(pairwise ~ fraction),                          # average across pools
    ##    mult_comp(pairwise ~ fraction | pool),                   # stratify over pools
    ##    mult_comp(trt.vs.ctrl ~fraction | pool, ref="tumour")    # All fraction_x vs tumour, per pool.  This would break on the 'no_cancer' dataset
    ##    mult_comp(pairwise ~ fraction * pool, interaction=TRUE)  # which genes have pool-dependent response to fraction.  All pairwise interactions.
    ##   )
    ## ),
    "Naive" = list(
      design = ~ fraction,
      comparison = mult_comp(pairwise ~ fraction)
  ),
  
  plot_scaler = function(y) {
    group_by(pool) %>%
      mutate(y =y / mean(y[fraction=="Control"]) ) %>%
      ungroup()
  }
)

#### Default if this file is missing, is
## list(
##   sample_sets = list(all=TRUE),
##   models=list(
##     "Naive" = list(
##       design = ~ covariate, # last column of design file
##       comparison = mult_comp(pairwise ~ covariate)
##     )
##   ),
##   plot_scale = function(y) {
##     y/mean(y)
##   }
## )

  
