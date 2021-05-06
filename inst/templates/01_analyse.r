#' ---
#' title: "{{{project}}}"
#' author: "{{{author}}}"
#' params:
#'   res_dir: "results"
#'   spec_file: ""
#'   spec_suffix: ""
#' output:
#'   bookdown::html_document2:
#'     toc: true
#'     code_folding: hide
#'     includes:
#'       in_header: styles.html
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
library(gt)
library(DESeq2)          # Assay-specific
library(PoiClaClu)
library(rtracklayer)
library(tximport)
library(clusterProfiler)
library(ReactomePA)
library(GO.db)
library(IHW)
library(DESdemonA)



fig_caption <- DESdemonA::captioner()

knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE,
                      dev=c("png","pdf"), out.width="90%",
                      results='asis', fig.cap=expression(fig_caption())
                      )
if (!isTRUE(getOption('knitr.in.progress'))) {
  params <- list(
    spec_file=dir(pattern="*.spec")[1],
    res_dir="results",
    spec_suffix="")
}

#+ read

data(list=paste0("dds", params$spec_suffix))

library(metadata(dds)$organism$org, character.only=TRUE)
metadata(dds)$template_git <- packageDescription("DESdemonA")$git_last_commit

specs   <- DESdemonA::load_specs(file=params$spec_file, context=dds)

param <- DESdemonA::ParamList$new(defaults=specs$settings)
# Calling the setter without a value will pick up the default (from the spec), _and_ alert it in the markdown, and return the default
set.seed(param$set("seed"))
param$set("title", "{{{project}}}")
param$set("script", file.path(getwd(),"01_analyse.r"))

ddsList <- DESdemonA::build_dds_list(dds, specs)


#Which features are zero across all samples in all datasets
all_zero <- Reduce(f=`&`,
                   x=lapply(ddsList, function(x) apply(counts(x)==0,1, all)),
                   init=TRUE)
ddsList <- lapply(ddsList, `[`, !all_zero)


ddsList <- lapply(ddsList, estimateSizeFactors)
param$set("baseMeanMin")
if (param$get("baseMeanMin")>0) {
  ddsList <- lapply(ddsList,
                   function(x) x[rowMeans(counts(x, normalized=TRUE)) >= param$get("baseMeanMin"),]
                   )
}


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
#+ qc-visualisation


ddsList <- lapply(ddsList,
                 DESdemonA::add_dim_reduct
                 )

param$set("top_n_variable")
for (dataset in names(ddsList)) {
  DESdemonA::qc_heatmap(
    ddsList[[dataset]], title=dataset,
    n=param$get("top_n_variable"),
    pc_x=1, pc_y=2,
    batch = ~1,
    caption=fig_caption
    )
}


#' 
#' # Find differential genes
#' 
#'  We use DESeq2 [@pkg_DESeq2] to find differential genes using the
#' negative binomial distribution to model counts, with IHW [@pkg_IHW]
#' for multiple-testing correction with greater power than
#' Benjamini-Hochberg, and ashr [@pkg_ashr] for effect-size shrinkage
#' to ensure reported fold-changes are more robust.
#'
#+ differential

param$set("alpha")
param$set("lfcThreshold")
param$set("filterFun")

## For each dataset, fit all its models
dds_model_comp <- map(ddsList, DESdemonA::fit_models, minReplicatesForReplace = Inf)

## Now put results in each 3rd level mcols(dds)$results
dds_model_comp <- map_depth(
  dds_model_comp, 3, DESdemonA::get_result,
  filterFun=param$get("filterFun"),
  alpha=param$get("alpha"),
  lfcThreshold=param$get("lfcThreshold")
)

dds_env <- new.env()
dds_name <- paste0(basename(tools::file_path_sans_ext(params$spec_file)),"_dds")
assign(dds_name, dds_model_comp, envir=dds_env)
save(list=dds_name,
     file=file.path("data", paste0(dds_name, ".rda")),
     envir=dds_env,
     eval.promises=FALSE
     )
rm(list=dds_name, envir=dds_env)
rm(dds_env)

xl_files <- DESdemonA::write_results(dds_model_comp, param, dir=params$res_dir)
iwalk(xl_files,
     ~cat("\n\n<a class=\"download-excel btn btn-primary\" href=\"", sub(paste0("^", params$res_dir, "/?"), "", .x), "\"> Open Spreadsheet '", .y,"'</a>", sep="")
     )
#DESdemonA::write_all_results(dds_model_comp, dir=params$res_dir)




#'
#' ## Summary Tables {.tabset}
#'
#' Here we summarise the results of the differential testing. Each
#' sample-set gets its own table (flick between tabs to choose which),
#' within which the different models, and their null hypotheses, are
#' enumerated.
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
#+ summary

summaries <- map_depth(dds_model_comp, 3, DESdemonA::summarise_results)


per_dataset <- map(summaries, DESdemonA::rbind_summary,
                  levels=c("Design","Comparison")) %>%
  map(function(x) {
    levels(x$Group) <- sub(".*<(.+<.+<0)$", "\\1", levels(x$Group))
    levels(x$Group) <- sub("^(0<.+?<.+?)<.*$", "\\1", levels(x$Group))
    levels(x$Group) <- sub("(.*?<)(.*0.*)(<.*)", "\\10\\3", levels(x$Group))
    droplevels(subset(x, Group!="Zero Count" & Group!="Low Count"))
  }
  )

pval_frame <- map_depth(
  dds_model_comp,
  2,
  ~map_dfr(.x, function(comp) {as.data.frame(mcols(comp)$results)}, .id="comparison")
)

for (dataset in names(per_dataset)) {
  options(htmltools.preserve.raw=TRUE)
  cat("### ", dataset, " \n", sep="")
  per_dataset[[dataset]] %>%
  gt(caption=paste0("Size of gene-lists for ", dataset)) %>%
    tab_header(title="Genelist summary",
               subtitle=dataset) %>%
    DESdemonA::tab_link_caption() %>%
    print()

  for (model in names(pval_frame[[dataset]])) {
    pl <- ggplot(pval_frame[[dataset]][[model]], aes(x=pvalue, fill=baseMean<1)) +
      geom_density(alpha=0.5) +
      facet_wrap(~comparison) +
      theme_bw()
    print(pl)
    fig_caption(caption=paste0("p-Diagnostic for", dataset, model))
  }
  
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
#' Again, select the dataset you wish to examine via the tabs - there
#' are sub-tabs for different experimental designs, if any
#' (e.g. removing/ignoring batch effects).
#'
#+ differential-MA

for (dataset in names(dds_model_comp)) {
  cat("## ", dataset, " {.tabset} \n", sep="") 
  for (mdl in names(dds_model_comp[[dataset]])) {
    cat("### ", mdl, "\n", sep="") 
    differential_MA(dds_model_comp[[dataset]][[mdl]], caption=fig_caption)
  }
}




#'
#' # Differential Heatmaps {.tabset}
#'
#' Here we present heatmaps of the differential genelists.  Note
#' that these will, by definition, divide the samples along lines
#' of the experimental groups, so should be interpreted differently
#' from the QC heatmaps which were blind to the experimental design.

#+ differential-heatmap
for (dataset in names(dds_model_comp)) {
  cat("## ", dataset, " {.tabset} \n", sep="") 
  for (mdl in names(dds_model_comp[[dataset]])) {
    cat("### ", mdl, "\n", sep="")
    DESdemonA::differential_heatmap(dds_model_comp[[dataset]][[mdl]],
                         . %>% mutate(.value=.value - mean(.value)),
                         caption=fig_caption
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
#+ enrich-init
param$set("showCategory")

#' ## Reactome Enrichment {.tabset} 
#'
#+ enrich-reactome,  eval=FALSE
enrich_plots <- map_depth(dds_model_comp, 2, DESdemonA::enrichment,
  fun="enrichPathway", showCategory = param$get("showCategory"), max_width=30)

enrich_plots <- map(enrich_plots, function(x) x[!sapply(x, length)==0])
for (dataset in names(enrich_plots)) {
    cat("### ", dataset, " {.tabset} \n", sep="") 
    for (mdl in names(enrich_plots[[dataset]])) {
      cat("#### ", mdl, " \n", sep="")
      print(enrich_plots[[dataset]][[mdl]]$plot)
      lbl <- paste("Reactome for", mdl, dataset)
      fig_caption(lbl)
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
#+ enrich-GO,  eval=FALSE
enrich_plots <- map_depth(dds_model_comp, 2, DESdemonA::enrichment,
  fun="enrichGO", showCategory = param$get("showCategory"), max_width=30)
enrich_plots <- map(enrich_plots, function(x) x[!sapply(x, length)==0])

for (dataset in names(enrich_plots)) {
    cat("### ", dataset, " {.tabset} \n", sep="") 
    for (mdl in names(enrich_plots[[dataset]])) {
      cat("#### ", mdl, " \n", sep="")
      print(enrich_plots[[dataset]][[mdl]]$plot)
      lbl <- paste("GO MF for", mdl, dataset)
      fig_caption(lbl)
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
