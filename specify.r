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
    "Naive" = list(
      design = ~ fraction,
      comparison = mult_comp(pairwise ~ fraction)
    )
  ),
  
  plot_scaler = function(y) {
    group_by(pool) %>%
      mutate(y =y / mean(y[fraction=="Control"]) ) %>%
      ungroup()
  }
)

