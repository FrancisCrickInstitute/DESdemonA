#' ---
#' title: "{{project}}"
#' author: "{{author}}"
#' params:
#'   res_dir: "results"
#' output:
#'   bookdown::html_document2:
#'     toc: true
#'     code_folding: hide
#'     includes:
#'       in_header: styles.html
#' bibliography: R.bib
#' link-citations: yes
#' ---
#'

#' # Preface
#'
#' We load in all the necessary additional R packages, and set up some initial parameters.
#+ init,results='hide', warning=FALSE, error=FALSE, message=FALSE

library(tidyverse)       # Language
library(devtools)
library(openxlsx)        # IO
library(ggrepel)         # General Vis.
library(grid)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(DESeq2)          # Assay-specific
library(PoiClaClu)
library(rtracklayer)
library(tximport)
library(clusterProfiler)
library(ReactomePA)
library(GO.db)
library("IHW")
devtools::load_all()


param <- {{package}}::ParamList$new(
  title="title",
  script="analyse.r",
  seed = 1)
set.seed(param$get("seed"))

caption <- {{package}}::captioner()
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE,
                      dev=c("png","pdf"), out.width="80%",
                      results='asis')

#+ read
data(rsem_dds)
library(metadata(rsem_dds)$organism$org, character.only=TRUE)


ddsList <- list(
  all = rsem_dds
)

#Which features are zero across all samples in all datasets
all_zero <- Reduce(`&`, lapply(ddsList, function(x) apply(counts(x)==0,1, all)))
ddsList <- lapply(ddsList, `[`, !all_zero)


ddsList <- lapply(ddsList, estimateSizeFactors)


#'
#' # QC Visualisation {.tabset}
#'
#' These visualisations are carried out blind to the experimental
#' design. For the heatmaps, we select a subset of genes removing the
#' (assumed uninformative) genes with flat expression levels across
#' the whole sample set. We'd expect samples from the same
#' experimental group to cluster together, in the sense that they are
#' on the same branch of the tree. But this is a transcriptome-wide
#' picture, and even if the clustering is not perfect, there will
#' still be genes that are consistent with the expected expression
#' profiles, and they will be revealed in the differential analysis.
#'
#' The second heatmap focuses on the pairwise similarity of samples,
#' using a slightly different metric than the first, to emphasise that
#' there's not one way of calculating how similar two samples are,
#' transcriptome-wide.  The colours in the first heatmap are a gene's
#' expression in a sample, relative to its average expression.  The
#' second heatmap represents the similarity of each pair of samples.
#'
#' In the third plot, we attempt to give meta-genes (principal
#' components) a meaning in terms of the experimental factors. The
#' principal components are expression profiles that best summarise
#' the overall behaviour, the first being the best single 'gene'
#' summary.  And the plot shows which of these components relate to
#' which experimental factors - this should be taken fairly loosely,
#' and primarily guides which other plots we should generate to look
#' at how samples associate with each other.
#'
#' The remaining plots examine these components in the context of the
#' various experimental factors, either looking at the first two
#' components (which account for the largest sample:sample variation)
#' or the components selected to have the greatest association with
#' the factor.
#'
#+ qc-visualisation, fig.cap=caption()
param$set("top_n_variable", 500, "Only use {} genes for unsupervised clustering")


for (dataset in names(ddsList)) {
  {{package}}::qc_heatmap(
    ddsList[[dataset]], title=dataset,
    n=param$get("top_n_variable"),
    pc_x=1, pc_y=2,
    batch=~1,
    caption=caption
    )
}


#' 
#' # Find differential genes
#' 
#' We want to account for any differences in individual mouse
#' backgrounds when testing for differential expression between
#' groups.  We include terms for _both_ Group and Mouse in our model
#' to enable this.  We use DESeq2 [@pkg_DESeq2] to find differential
#' genes using the negative binomial distribution to model counts,
#' with IHW [@pkg_IHW] for multiple-testing correction with greater
#' power than Benjamini-Hochberg, and ashr [@pkg_ashr] for effect-size
#' shrinkage to ensure reported fold-changes are more robust.
#'
#+ differential, fig.cap=caption()

param$set("alpha", 0.05)
param$set("lfcThreshold", 0)


## Example set of designs and their comparisons.

mdlList <- list(
  "Standard" = list(
    design = ~ treatment,
    comparison = list(
      "DPN" = c("treatment","DPN","Control"),
      "OHT" = c("treatment","OHT","Control")
      )
  ),
  "Patient Adjusted" = list(
    design = ~ patient + treatment, 
    comparison = list(
      "Treated" = list(c("treatment_DPN_vs_Control", "treatment_OHT_vs_Control"), listValues=c(.5, -1)),
      mult_comp(revpairwise~treatment),
      "Patient" = ~treatment
    )
  )
)

## Each dataset gets the same set of models
ddsList <- map(ddsList, function(dds) {
  metadata(dds)$model <- mdlList
  dds
})

# or you could explicitly allocate a different set of
# models to each dataset separately

## For each dataset, fit all its models
dds_model_comp <- map(ddsList, {{package}}::fit_models, minReplicatesForReplace = Inf)

## Now put results in each 3rd level mcols(dds)$results
dds_model_comp <- map_depth(
  dds_model_comp, 3, {{package}}::get_result,
#  filterFun=NULL, # reinstate if we get IHR working
  alpha=param$get("alpha"))

xl_files <- {{package}}::write_results(dds_model_comp, param, dir=params$res_dir)
iwalk(xl_files,
     ~cat("\n\n<a class=\"download-excel btn btn-primary\" href=\"", sub(paste0("^", params$res_dir, "/?"), "", .x), "\"> Open Spreadsheet '", .y,"'</a>", sep="")
     )
#{{package}}::write_all_results(dds_model_comp, dir=params$res_dir)


#'
#' ## Summary Tables {.tabset}
#'
#' Here we summarise the results of the differential testing. As
#' mentioned above, there is only one design strategy, but we have a
#' choice as to which samples we use (determined by which tab we
#' select below), and which null hypothesis we test against (described
#' in the 'Comparison' column.)
#'
#' In the 'Significant' column we tally the number of significant
#' genes (for pairwise comparisons, separated into up or down, where
#' A-B being labelled up means expression is higher in A; for omnibus
#' comparisons, a broad categorisation of the most extreme groups, as
#' described above.) We also tally the total number of genes
#' exhibiting that behaviour (ie not necessarily statistically
#' signficant) - these might not always add up to the same value, as
#' there is an independent filter of low-signal genes whose effect
#' varies from comparison to comparison.
#' 
#+ summary, fig.cap=caption()

summaries <- map_depth(dds_model_comp, 3, {{package}}::summarise_results)


per_dataset <- map(summaries, {{package}}::rbind_summary,
                  levels=c("Design","Comparison")) %>%
  map(function(x) {
    levels(x$Group) <- sub(".*<(.+<.+<0)$", "\\1", levels(x$Group))
    levels(x$Group) <- sub("^(0<.+?<.+?)<.*$", "\\1", levels(x$Group))
    levels(x$Group) <- sub("(.*?<)(.*0.*)(<.*)", "\\10\\3", levels(x$Group))
    droplevels(subset(x, Group!="Zero Count" & Group!="Low Count"))
  }
  )

    

for (dataset in names(per_dataset)) {
  cat("### ", dataset, " \n", sep="")
  print(knitr::kable(per_dataset[[dataset]],
                     options=list(dom="t"),
                     caption=paste0("Size of differential gene-lists for ", dataset),
                     label=klabel(paste0("summary", dataset))))
}


#' # Scatterplot of differential genes {.tabset}
#'
#' Here we plot the fold-change across all genes on the vertical axis,
#' against their overall expression on the horizontal. Statistically
#' significant genes are highlighted in colour, and we attempt to
#' auto-label as many gene symbols as might be legible.
#'
#' The fold-changes are regularised ('shrunk') to reflect their
#' reliability: genes with large variability between replicates are
#' moved towards the horizontal axis, to de-emphasise them and to
#' provide a more reliable prediction of behaviour in future
#' replications.
#'
#' Again, select the dataset you wish to examine via the tabs - there are
#' sub-tabs for different experimental designs, if any (e.g. removing/ignoring
#' batch effects).
#'
#+ differential-MA, fig.cap=caption()

for (dataset in names(dds_model_comp)) {
  cat("## ", dataset, " {.tabset} \n", sep="") 
  for (mdl in names(dds_model_comp[[dataset]])) {
    cat("### ", mdl, "\n", sep="") 
    differential_MA(dds_model_comp[[dataset]][[mdl]], caption=caption)
  }
}




#'
#' # Differential Heatmaps {.tabset}
#'
#' Here we present heatmaps of the differential genelists.  Note
#' that these will, by definition, divide the samples along lines
#' of the experimental groups, so should be interpreted differently
#' from the QC heatmaps which were blind to the experimental design.

#+ differential-heatmap, fig.cap=caption()
for (dataset in names(dds_model_comp)) {
  cat("## ", dataset, " {.tabset} \n", sep="") 
  for (mdl in names(dds_model_comp[[dataset]])) {
    cat("### ", mdl, "\n", sep="")
    {{package}}::differential_heatmap(dds_model_comp[[dataset]][[mdl]],
                         . %>% rownames_to_column() %>%
                           mutate(.value=.value - mean((.value[Diagnosis=="CONTROL"]))) %>%
                           dplyr::select(Diagnosis, .value, rowname) %>%
                           column_to_rownames(),
                         caption=caption
                         )
  }
}


#' # Terms Of Use
#'
#' The Crick has a [publication
#' policy](https://intranet.crick.ac.uk/our-crick/library-information-services/pages/guidelines-publication)
#' and we expect to be included on publications, regardless of funding
#' arrangements. Any use of these results in publication must be
#' discussed with BABS regarding authorship. If not authorship then
#' the BABS analyst must receive a named acknowledgement. Please also
#' cite the following sources which have enabled the analysis to be
#' carried out.
#'
#' # Bibliography
#'
#' 
knitr::knit_exit()

#'
#' # Enrichment Analysis
#'
#' Here we look at which [Reactome](https://reactome.org/) and [GO molecular functions](http://geneontology.org/)
#' are enriched in the various genelists we have created.
#'
#+ enrich-init, fig.cap=caption()
param$set("showCategory", 25, "Only show top {} enriched categories in plots")

#' ## Reactome Enrichment {.tabset} 
#'
#+ enrich-reactome, fig.cap=caption(), eval=FALSE
enrich_plots <- map_depth(dds_model_comp, 2, {{package}}::enrichment,
  fun="enrichPathway", showCategory = param$get("showCategory"), max_width=30)

enrich_plots <- map(enrich_plots, function(x) x[!sapply(x, length)==0])
for (dataset in names(enrich_plots)) {
    cat("### ", dataset, " {.tabset} \n", sep="") 
    for (mdl in names(enrich_plots[[dataset]])) {
      cat("#### ", mdl, " \n", sep="")
      print(enrich_plots[[dataset]][[mdl]]$plot)
      lbl <- paste("Reactome for", mdl, dataset)
      caption(lbl)
      print(
        knitr::kable(enrich_plots[[dataset]][[mdl]]$table,
                     options=list(dom="t"),
                     caption=lbl,
                     label=klabel(paste("REACTOME", dataset, mdl))
                     )
      )
    }
}

#' ## GO MF Enrichment {.tabset} 
#'
#+ enrich-GO, fig.cap=caption(), eval=FALSE
enrich_plots <- map_depth(dds_model_comp, 2, {{package}}::enrichment,
  fun="enrichGO", showCategory = param$get("showCategory"), max_width=30)
enrich_plots <- map(enrich_plots, function(x) x[!sapply(x, length)==0])

for (dataset in names(enrich_plots)) {
    cat("### ", dataset, " {.tabset} \n", sep="") 
    for (mdl in names(enrich_plots[[dataset]])) {
      cat("#### ", mdl, " \n", sep="")
      print(enrich_plots[[dataset]][[mdl]]$plot)
      lbl <- paste("GO MF for", mdl, dataset)
      caption(lbl)
      print(
        knitr::kable(enrich_plots[[dataset]][[mdl]]$table,
                     options=list(dom="t"),
                     caption=lbl,
                     label=klabel(paste("GO MF", dataset, mdl))
                     )
        )
    }
}

            
  



#' # Terms Of Use
#'
#' The Crick has a [publication
#' policy](https://intranet.crick.ac.uk/our-crick/library-information-services/pages/guidelines-publication)
#' and we expect to be included on publications, regardless of funding
#' arrangements. Any use of these results in publication must be
#' discussed with BABS regarding authorship. If not authorship then
#' the BABS analyst must receive a named acknowledgement. Please also
#' cite the following sources which have enabled the analysis to be
#' carried out.
#'
#' # Bibliography
#'
#' 
