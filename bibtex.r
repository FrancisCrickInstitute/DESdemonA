library(tidyverse)       # Language
library(devtools)
library(openxlsx)        # IO
library(ggrepel)         # General Vis.
library(grid)
library(RColorBrewer)
library(ComplexHeatmap)
library(DESeq2)          # Assay-specific
library(PoiClaClu)
library(rtracklayer)
library(tximport)
library(clusterProfiler)
library(ReactomePA)
library(GO.db)
library(IHW)
library(ashr)
my_bib <- function(cite) {
  bib <- sub("\\{,",
            paste0("{pkg_", cite, ","),
      toBibtex(citation(cite)[1]))
  bib <- bib[!grepl("note = \\{R package version", bib)]
  bib
}
write(my_bib("base"), file="R.bib")
for (i in setdiff(.packages(), "babsRNASeq")) {
  write(my_bib(i), file="R.bib", append=TRUE)
}
