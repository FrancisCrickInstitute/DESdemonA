#' ---
#' title: "{{{project}}}"
#' author: "{{{author}}}"
#' params:
#'   res_dir: "results"
#'   spec_file: "{{{specfile}}}"
#'   count_source: "{{{count_source}}}"
#'   param_call: !r as.list(sys.call(1))$params
#' output:
#'   bookdown::html_document2:
#'     toc: true
#'     code_folding: hide
#'     includes:
#'       in_header: !expr system.file("templates/styles.html", package="DESdemonA")
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
library(RColorBrewer)
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
library(equatiomatic)
library(DESdemonA)



fig_caption <- DESdemonA::captioner()

knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE,
                      dev=c("ragg_png","pdf"), out.width="90%",
                      fig.width=14, fig.height=10,
                      results='asis', fig.cap=expression(fig_caption())
                      )
if (!isTRUE(getOption('knitr.in.progress'))) {
  params <- list(
    spec_file=dir(pattern="*.spec")[1],
    res_dir="results",
    count_source=gsub("counts_(.*).rda", "\\1", dir("data", pattern="counts_.*.rda")[1]))
}

#+ read

if (is.character(params$count_source)) {
  data(list=paste0("counts_", params$count_source))
} else {
  dds <- params$count_source
}

if (!is.null(metadata(dds)$organism$org)) {
  library(metadata(dds)$organism$org, character.only=TRUE)
}
metadata(dds)$template_git <- packageDescription("DESdemonA")$git_last_commit

specs   <- DESdemonA::load_specs(file=params$spec_file, context=dds)

param <- DESdemonA::ParamList$new(defaults=specs$settings)
# Calling the setter without a value will pick up the default (from the spec), _and_ alert it in the markdown, and return the default
set.seed(param$set("seed"))
param$set("title", "{{{project}}}")
param$set("script", file.path(getwd(),"01_analyse.r"))
param$set("spec", sub("\\.spec$", "", params$spec_file), "Using analysis plan '{}'.")
if (is.character(params$count_source)) {
  param$set("count_source", params$count_source, "Using alignment settings '{}'")
} else {
  param$set("count_source", params$param_call$count_source, "Using counts from '{}'")
}  
  
ddsList <- DESdemonA::build_dds_list(dds, specs)


#Which features are zero across all samples in all datasets
all_zero <- Reduce(f=`&`,
                  x=map_des(ddsList, function(x) apply(counts(x)==0,1, all)),
                  init=TRUE)
ddsList <- map_des(ddsList, function(dds) dds[!all_zero,])



ddsList <- map_des(ddsList, estimateSizeFactors)

param$set("baseMeanMin")
if (param$get("baseMeanMin")>0) {
  ddsList <- map_des(ddsList,
                   function(x) x[rowMeans(counts(x, normalized=TRUE)) >= param$get("baseMeanMin"),]
                   )
}

#' # Input Summary
#'
#' The sample annotations are as follows:
#'
#+ inputs
DESdemonA::table1(dds, ddsList) %>%
  print()

#' We may examine the samples in different combinations, and leave out
#' certain samples.  In the above table, the columns under the 'In
#' Subset' group indicate these combinations are listed, and the
#' samples' inclusions are indicated.  Each of those combinations may be analysed in
#' potentially several ways, as formulated here:
#'

for (dataset in names(specs$sample_sets)) {
  mdls <- lapply(specs$sample_sets[[dataset]]$models, function(x) x$design)
  mdls <- mdls[sapply(specs$sample_sets[[dataset]]$model, function(x) "comparisons" %in% names(x))]
  for (mdl in names(mdls)) {
    cat("\n\n####", dataset, mdl, "{-}\n\n")
    print(extract_eq(lm(
      update(mdls[[mdl]], Expression ~ .),
      cbind(as.data.frame(colData(ddsList[[dataset]])), Expression=1:ncol(ddsList[[dataset]]))
    )
    ))
  }
}



#' 
#' We also have the following defaults set.  Whenever they are used
#' or changed in the, analysis that will be highlighted in the text
#' alongside a 'wrench' icon:
#' 
#+ settings
data.frame(
  Option=names(specs$settings),
  Value=sapply(specs$settings, deparse)
) %>%
  gt(caption="Analysis settings") %>%
  DESdemonA::tab_link_caption() %>%
  print()

 


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


ddsList <- map_des(ddsList, function(x) DESdemonA::add_dim_reduct(x))

for (assay in c("raw", "norm", "vst")) {
  DESdemonA::write_assay(ddsList, assay=assay, path=params$res_dir)
}
                       
param$set("top_n_variable")
param$set("clustering_distance_rows")
param$set("clustering_distance_columns")

map_des(ddsList, 
            ~DESdemonA::qc_heatmap(.x,
              title=dataset_name(.x),
              caption=fig_caption,
              param=param$publish()
            )
            )



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
dds_model_comp <- map_des(
  ddsList,
  function(x) DESdemonA::fit_models(x,
                             minReplicatesForReplace = Inf)
)


## Now put results in each 3rd level mcols(dds)$results
dds_model_comp <- map_des(
  dds_model_comp,
  DESdemonA::get_result,
  filterFun=param$get("filterFun"),
  alpha=param$get("alpha"),
  lfcThreshold=param$get("lfcThreshold")
)



dds_name <- paste0(basename(tools::file_path_sans_ext(params$spec_file)),"_x_", params$count_source)
save(dds_model_comp,
     file=file.path("data", paste0(dds_name, ".rda")),
     eval.promises=FALSE
     )
saveRDS(dds_model_comp, file=file.path("data", paste0(dds_name, ".rds")))


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


summaries <- dds_model_comp %>%
  map_des(DESdemonA::summarise_results) %>%
  map_des(function(x) {bind_rows(x, .id="Comparison")}, depth="model") %>%
  map_des(function(x) {bind_rows(x, .id="Model")}, depth="dataset") %>%
  map(function(x) {
    levels(x$Group) <- sub(".*<(.+<.+<0)$", "\\1", levels(x$Group))
    levels(x$Group) <- sub("^(0<.+?<.+?)<.*$", "\\1", levels(x$Group))
    levels(x$Group) <- sub("(.*?<)(.*0.*)(<.*)", "\\10\\3", levels(x$Group))
    droplevels(subset(x, Group!="Zero Count" & Group!="Low Count"))
  }
  )

pval_frame <- dds_model_comp %>%
  map_des(~as.data.frame(mcols(.x)$results)) %>%
  map_des(function(x) {bind_rows(x, .id="Comparison")}, depth="model")
                      
                      

for (dataset in names(summaries)) {
  options(htmltools.preserve.raw=TRUE)
  cat("### ", dataset, " \n", sep="")
  summaries[[dataset]] %>%
    bookdown_label(dataset) %>%
    gt(caption=paste0("Size of gene-lists for ", dataset),
       groupname_col="Model") %>%
    tab_header(title="Genelist summary",
               subtitle=dataset) %>%
    DESdemonA::tab_link_caption() %>%
    print() %>%
    bookdown_label()

  for (model in names(pval_frame[[dataset]])) {
    pl <- ggplot(pval_frame[[dataset]][[model]], aes(x=pvalue, fill=baseMean<1)) +
      geom_density(alpha=0.5) +
      facet_wrap(~Comparison) +
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
    DESdemonA::differential_heatmap(
      ddsList=dds_model_comp[[dataset]][[mdl]],
      param=param$publish(),
      caption=fig_caption
    )
  }
}



#'
#' # Over-representation Analysis
#'
#' Here we look at which [Reactome](https://reactome.org/) and [GO molecular functions](http://geneontology.org/)
#' are over-represented in the various genelists we have created.
#'
param$set("showCategory")

#' ## Reactome Over-representation {.tabset} 
#'
#+ over-reactome

over_rep_plots <-
  map_des(dds_model_comp,
          ~DESdemonA::over_representation(.x, fun="enrichPathway", showCategory = param$get("showCategory"), max_width=30),
          depth="model")

for (dataset in names(over_rep_plots)) {
  has_over_rep <- !sapply(over_rep_plots[[dataset]], is.null)
  if (!any(has_over_rep)) {
    next
  }
  cat("### ", dataset, " {.tabset} \n", sep="")
  for (mdl in names(over_rep_plots[[dataset]][has_over_rep])) {
    cat("#### ", mdl, " \n", sep="")
    print(over_rep_plots[[dataset]][[mdl]]$plot)
    lbl <- paste("Reactome for", mdl, dataset)
    fig_caption(lbl)
    over_rep_plots[[dataset]][[mdl]]$table %>%
      bookdown_label(dataset) %>%
      gt(caption=lbl) %>%
      tab_header(title="Reactome",
                 subtitle=paste(dataset, mdl)) %>%
      DESdemonA::tab_link_caption() %>%
      print() %>%
      bookdown_label()
  }
}

#' ## GO MF Over-representation {.tabset} 
#'
#+ over-GO

over_rep_plots <-
  map_des(dds_model_comp,
          ~DESdemonA::over_representation(.x, fun="enrichGO", showCategory = param$get("showCategory"), max_width=30),
          depth="model")

for (dataset in names(over_rep_plots)) {
  has_over_rep <- !sapply(over_rep_plots[[dataset]], is.null)
  if (!any(has_over_rep)) {
    next
  }
  cat("### ", dataset, " {.tabset} \n", sep="")
  for (mdl in names(over_rep_plots[[dataset]][has_over_rep])) {
    cat("#### ", mdl, " \n\n", sep="")
    print(over_rep_plots[[dataset]][[mdl]]$plot)
    lbl <- paste("GO MF for", dataset, mdl)
    fig_caption(lbl)
    over_rep_plots[[dataset]][[mdl]]$table %>%
      bookdown_label(dataset) %>%
      gt(caption=lbl) %>%
      tab_header(title="GO Molecular Function",
                 subtitle=paste(dataset, mdl)) %>%
      DESdemonA::tab_link_caption() %>%
      print() %>%
      bookdown_label()
  }
}



#'
#' # Enrichment Analysis
#'
#' Here we look at which [Reactome](https://reactome.org/) and [GO
#' molecular functions](http://geneontology.org/) are enriched in the
#' various genelists we have created. Enrichment analysis looks at the
#' fold-changes of the genes of interest (pathway/go term/...) and sees if they
#' are different from the fold-changes of the remaining genes.
#'

#' ## Reactome Enrichment {.tabset} 
#'
#+ enrich-reactome

enrich <- map_des(dds_model_comp,
                        ~DESdemonA::enrichment(.x, fun="gsePathway")
                 ) %>%
  DESdemonA:::trim_map()
  

for (dataset in names(enrich)) {
  cat("### ", dataset, "\n\n")
  for (model in names(enrich[[dataset]])) {
    cat("#### ", model, "\n\n")
    for (comparison in names(enrich[[dataset]][[model]])) {
      obj <- enrich[[dataset]][[model]][[comparison]]
      as.data.frame(obj) %>%
        bookdown_label(paste(dataset, model, comparison, sep="-")) %>%
        gt(caption=comparison) %>%
        cols_hide(columns=c(core_enrichment)) %>%
        fmt_scientific(
          columns = c(pvalue, p.adjust, qvalues),
          decimals=2) %>%
        fmt_number(
          columns = c(enrichmentScore, NES),
          decimals=2) %>%
        tab_header(title="Reactome Enrichment",
                   subtitle=paste(dataset, model, comparison)) %>%
        DESdemonA::tab_link_caption() %>%
        print() %>%
        bookdown_label()
    }
  }
}

#' ## GO Enrichment {.tabset} 
#'
#+ enrich-go

enrich <- map_des(dds_model_comp,
                        ~DESdemonA::enrichment(.x, fun="gseGO")
                        ) %>%
  DESdemonA:::trim_map()
  

for (dataset in names(enrich)) {
  cat("### ", dataset, "\n\n")
  for (model in names(enrich[[dataset]])) {
    cat("#### ", model, "\n\n")
    for (comparison in names(enrich[[dataset]][[model]])) {
      obj <- enrich[[dataset]][[model]][[comparison]]
      as.data.frame(obj) %>%
        bookdown_label(paste(dataset, model, comparison, sep="-")) %>%
        gt(caption=comparison) %>%
        cols_hide(columns=c(core_enrichment)) %>%
        fmt_scientific(
          columns = c(pvalue, p.adjust, qvalues),
          decimals=2) %>%
        fmt_number(
          columns = c(enrichmentScore, NES),
          decimals=2) %>%
        tab_header(title="GO Enrichment",
                   subtitle=paste(dataset, model, comparison)) %>%
        DESdemonA::tab_link_caption() %>%
        print() %>%
        bookdown_label()
    }
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
