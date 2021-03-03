specification(
    settings=settings(     ## analysis parameters
	alpha          = 0.01,    ## p-value cutoff
	lfcThreshold   = 0,       ## abs lfc threshold
	baseMeanMin    = 1,       ## discard transcripts with average normalised counts lower than this
	top_n_variable = 500,     ## For PCA
	showCategory   = 25,      ## For enrichment analyses
	seed           = 1,       ## random seed gets set at start of script, just in case.
	filterFun      = IHW::ihw ## NULL for standard DESeq2 results, otherwise  functions
    ),
  sample_sets =list( 
    all = sample_set(
      subset = TRUE,
      models = list(
        accurate = model(
          design = ~ stage *  genotype + batch,
          comparisons = list(
            mult_comp(revpairwise ~ stage | genotype  ),
            mult_comp(revpairwise ~ genotype  | stage)
          )
        ),
        powerful = model(
          design = ~stage *  genotype,
          comparisons = list(
            mult_comp(revpairwise ~ stage | genotype  ),
            mult_comp(revpairwise ~ genotype  | stage)
          )
        )
      )
    ),
    batch1=sample_set(
      subset = batch == "1",
      models = list( 
        powerful = model(
          design = ~stage *  genotype,
          comparisons = list(
            mult_comp(revpairwise ~ stage | genotype  ),
            mult_comp(revpairwise ~ genotype  | stage)
          )
        )
        )
    ),
    batch2=sample_set(
      subset = batch == "2",
      models = list( 
        powerful = model(
          design = ~stage *  genotype,
          comparisons = list(
            mult_comp(revpairwise ~ stage | genotype  ),
            mult_comp(revpairwise ~ genotype  | stage)
          )
          )
        )
      )
    )
  )

