#' ---
#' title: "{{title}}"
#' author: "{{author}}"
#' output:
#'   bookdown::html_document2:
#'     toc: true
#'     code_folding: hide
#'     css: styles.css
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


param <- babsRNASeq::ParamList$new(
  title="title",
  script="analyse.r",
  seed = 1)
set.seed(param$get("seed"))

caption <- babsRNASeq::captioner()
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
#+ qc-visualisation, fig.cap=caption()
param$set("top_n_variable", 500, "Only use {} genes for unsupervised clustering and PCA")


for (dataset in names(ddsList)) {
  babsRNASeq::qc_heatmap(
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
#' to enable this.  We use [DESeq2 @pkg_DESeq2] to find differential
#' genes using the negative binomial distribution to model counts,
#' with [IHW @pkg_IHW] for multiple-testing correction with greater
#' power than Benjamini-Hochberg, and [ashr @pkg_ashr] for effect-size
#' shrinkage to ensure reported fold-changes are more robust.
#'
#' #+ differential, fig.cap=caption()

param$set("alpha", 0.05)
param$set("lfcThreshold", 0)
cat("\n\n## Summary Tables {.tabset}\n\n")


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
      "DPN" = c("treatment","DPN","Control"),
      "OHT" = c("treatment","OHT","Control"),
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
dds_model_comp <- map(ddsList, babsRNASeq::fit_models)
## Now put results in each 3rd level mcols(dds)$results
dds_model_comp <- map_depth(
  dds_model_comp, 3, babsRNASeq::get_result,
  alpha=param$get("alpha"))


## Summarise results
summaries <- map_depth(dds_model_comp, 3, babsRNASeq::summarise_results)


per_dataset <- map(summaries, babsRNASeq::rbind_summary,
    levels=c("Design","Comparison"))
for (dataset in names(per_dataset)) {
  cat("### ", dataset, " \n", sep="")
  print(knitr::kable(per_dataset[[dataset]],
                     options=list(dom="t"),
                     caption=paste0("Size of differential gene-lists for ", dataset),
                     label=klabel(paste0("summary", dataset))))
}




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
#+ enrich-reactome, fig.cap=caption()
enrich_plots <- map_depth(dds_model_comp, 2, babsRNASeq::enrichment,
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
#+ enrich-GO, fig.cap=caption()
enrich_plots <- map_depth(dds_model_comp, 2, babsRNASeq::enrichment,
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

            
  

#' # Scatterplot of differential genes {.tabset}
#'
#+ differential-MA, fig.cap=caption()

for (dataset in names(dds_model_comp)) {
  cat("## ", dataset, " {.tabset} \n", sep="") 
  for (mdl in names(dds_model_comp[[dataset]])) {
    cat("### ", mdl, "\n", sep="") 
    differential_MA(dds_model_comp[[dataset]][[mdl]], caption=caption)
  }
}



# This section very sketchy, not yet abstracted?

#'
#' # Differential Heatmaps
#'
#' Here we present heatmaps of the differential genelists.  Note
#' that these will, by definition, divide the samples along lines
#' of the exprimental groups, so should be interpreted differently
#' from the QC heatmaps which were blind to the experimental design.

#+ differential-heatmap, fig.cap=caption()

if (FALSE) {
differential_heatmap(dds_model_comp$all$Pooled,
                     . %>% group_by(Mouse) %>%
                       mutate(.value=.value - (.value[Group=="1"])) %>%
                       dplyr::select(Group, Treatment, Lympho, Mouse, .value) %>%
                       group_by(Group) %>%
                       arrange(Mouse),
                     title="All Samples")
}


## *** Output filtered results

xl_files <- babsRNASeq::write_results(dds_model_comp, param)
iwalk(xl_files,
     ~cat("\n\n[Download Spreadsheet '", .y,"'](", .x, ")\n\n", sep="")
     )
#babsRNASeq::write_all_results(dds_model_comp)
