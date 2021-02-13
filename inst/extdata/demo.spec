specification(
  sample_sets =  list( # Different combinations of samples
    all = sample_set(  # 'all' corresponds to all samples going into the dds object together
      subset=TRUE,     # if a sample evaluate's to TRUE then include it
      models=list(                                # We may have one or more design formulae
            "Batch-corrected" = model(            # specify each as
              design = ~ fraction + pool,         # a design formula
              comparisons = list(                 # and a list of comparisons (plurals are always lists)
                mult_comp(pairwise ~ fraction)    # all pairs of fractions are compared
              )
            )
      )
    ),
    no_cancer = sample_set(                   # Another grouping of samples
      subset = fraction != "tumour",           # this time, drop 'tumour'-fraction samples
      models = list(
        "Batch-corrected" = model(            # etc
          design = ~ fraction + pool,
          comparisons = list(
            mult_comp(pairwise ~ fraction)
          )
        ),
        "Naive" = model(
          design = ~ fraction,                # a different model, so adjusting for different terms
          comparisons = list(
            mult_comp(pairwise ~ fraction)
          )
        )
      )
    )
  ),
  models=list(),         ## will be added to every sample_set
  settings=settings(     ## analysis parameters
      alpha          = 0.05,    ## p-value cutoff
      lfcThreshold   = 0,       ## abs lfc threshold
      baseMeanMin    = 5,       ## discard transcripts with average normalised counts lower than this
      top_n_variable = 500,     ## For PCA
      showCategory   = 25,      ## For enrichment analyses
      seed           = 1,       ## random seed gets set at start of script, just in case.
      filterFun      = IHW::ihw ## NULL for standard DESeq2 results, otherwise  functions
      ),
  plot_scaler = function(y) { ## not currently called so no need to change
    group_by(pool) %>%
      mutate(.value = .value / mean(.value[fraction=="Control"]) ) %>%
      ungroup()
  }
)


#### Example 'comparison's for more complex design
## model(
##   design = ~ fraction + pool, 
##   comparisons = list(                          # can be a list
##     "f1 v f2" = list("fraction", "f2","f1"), # of DESeq2::results(contrast=e.g.
##     "fraction" = ~ pool,                     # or LRT full=design, reduced=this
##     mult_comp(pairwise ~ fraction | pool)     # and this will get expanded out to all pairwise
##   )
## )

## model(
##   design = ~ fraction * pool,                               # No replicates, but just imagine if there were
##   comparisons = list(
##     mult_comp(pairwise ~ fraction),                          # average across pools
##     mult_comp(pairwise ~ fraction | pool),                   # stratify over pools
##     mult_comp(trt.vs.ctrl ~fraction | pool, ref="tumour")    # All fraction_x vs tumour, per pool.  This would break on the 'no_cancer' dataset
##     mult_comp(pairwise ~ fraction * pool, interaction=TRUE)  # which genes have pool-dependent response to fraction.  All pairwise interactions.
##   )
## )

