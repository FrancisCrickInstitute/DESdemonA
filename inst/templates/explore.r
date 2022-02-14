#' ---
#' title: "desdemona"
#' author: ""
#' params:
#'   res_dir: "results"
#'   spec_file: ""
#'   spec_suffix: ""
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
library(grid)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(gt)
library(DESeq2)          # Assay-specific
library(PoiClaClu)
library(rtracklayer)
library(tximport)
library(GO.db)
library(IHW)
library(DESdemonA)


for (spec_file in dir(pattern="*.spec")) {
params <- list(
  spec_file=spec_file,
  res_dir="results",
  spec_suffix="")


#+ read
sample_sheet <- read.csv("experiment_table.csv", header=TRUE, check.names=FALSE, stringsAsFactors=TRUE)
dds <- makeExampleDESeqDataSet(
  n = 1000, m = nrow(sample_sheet),
  betaSD = 0, interceptMean = 4, interceptSD = 2, dispMeanRel = function(x) 4/x + 0.1,
  sizeFactors = rep(1, nrow(sample_sheet))
)
mcols(dds)$entrez <-
  mcols(dds)$symbol <-
  row.names(dds)

design(dds) <- ~1
colData(dds) <- DataFrame(sample_sheet)
specs   <- DESdemonA::load_specs(file=params$spec_file, context=dds)

param <- DESdemonA::ParamList$new(defaults=specs$settings)
# Calling the setter without a value will pick up the default (from the spec), _and_ alert it in the markdown, and return the default
set.seed(param$set("seed"))
param$set("title", "desdemona")
param$set("script", file.path(getwd(),"01_analyse.r"))
param$set("spec", sub("\\.spec$", "", params$spec_file), "Using analysis plan '{}'.")
param$set("spec_suffix", params$spec_suffix, "Using alignment settings '{}'")
  
ddsList <- DESdemonA::build_dds_list(dds, specs)


#Which features are zero across all samples in all datasets
all_zero <- Reduce(f=`&`,
                  x=map_des(ddsList, function(x) apply(counts(x)==0,1, all)),
                  init=TRUE)
ddsList <- map_des(ddsList, function(dds) dds[!all_zero,])
ddsList <- map_des(ddsList, estimateSizeFactors)


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

explore <- map_des(
  dds_model_comp,
  function(dds) list(
    metadata=metadata(dds),
    df=as.data.frame(colData(dds))
  )
)
saveRDS(explore, file=sub(".spec$", "_explore.rds", spec_file))
}
