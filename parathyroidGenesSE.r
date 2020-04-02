library(tidyverse)
library(tximport)
library(DESeq2)
library(parathyroidSE)
devtools::load_all()

data(parathyroidGenesSE)
rsem_dds <- DESeqDataSet(parathyroidGenesSE, design=~1)
rsem_dds <- collapseReplicates(rsem_dds, interaction(dds$patient, dds$treatment, dds$time))



align_settings <- c("genus: Homo", "species: Sapiens")
organism <- list(genus = grep("genus: .*", align_settings, value=TRUE),
                species = grep("species: .*", align_settings, value=TRUE))
organism <- map(organism, ~ gsub(".*: ", "", .))
organism$Gs <- paste0(toupper(substr(organism$genus, 1, 1)), tolower(substr(organism$species, 1, 1)))
organism$org <-  unique(grep(paste0("org\\.",organism$Gs),row.names(installed.packages()),  value=TRUE))
if (length(organism$org)!=1) {
  stop(paste("Can't find annotation for", organism$genus, organism$species))
}
metadata(rsem_dds)$organism <- organism

library(organism$org, character.only=TRUE)
mcols(rsem_dds)$symbol <- mapIds(
  eval(parse(text = organism$org)),
  keys=row.names(rsem_dds),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first")[row.names(rsem_dds)]

mcols(rsem_dds)$entrez <- mapIds(
  eval(parse(text = organism$org)),
  keys=row.names(rsem_dds),
  column="ENTREZID",
  keytype="ENSEMBL",
  multiVals="first")[row.names(rsem_dds)]


metadata(rsem_dds)$labels <- c("patient", "treatment","time")


usethis::use_data(rsem_dds, overwrite=TRUE)


