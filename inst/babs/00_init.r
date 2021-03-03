library(readxl)
library(tidyverse)
library(tximport)
library(DESeq2)
library(DESdemonA)

sample_data <- readxl::read_excel(system.file("extdata/{{metadata_file}}.xlsx", package="babsrnaseq"), sheet=1)

sample_sheet <- sample_data %>%
  dplyr::mutate(
           rsem = file.path(system.file("extdata/{{rsem}}", package="babsRNASeq"),
                            paste0(LIMSID, "_R1.genes.results")) # that annoying nfcore '_R1' suffix may disappear
         )

sample_sheet[sapply(sample_sheet, is.character)] <- lapply(sample_sheet[sapply(sample_sheet, is.character)], 
                                       as.factor)

txi <- tximport(as.character(sample_sheet$rsem), type="rsem")
txi$length[txi$length==0] <- 1
rsem_dds <- DESeqDataSetFromTximport(txi, sample_sheet, ~ 1) # ok to have fixed trivial design - it'll get swapped out
ind <- rowSums(counts(rsem_dds)!=0) > 0
rsem_dds <- rsem_dds[ind,]

txi[c("abundance", "counts", "length")] <- map(txi[c("abundance", "counts", "length")], ~ .[ind,])


organism <- list(genus   = "{{genus}}",
                 species = "{{species}}")
organism$Gs <- paste0(toupper(substr(organism$genus, 1, 1)),
                     tolower(substr(organism$species, 1, 1)))
organism$org <-  unique(grep(paste0("org\\.",organism$Gs),row.names(installed.packages()),  value=TRUE))
if (length(organism$org)!=1) {
  stop(paste("Can't find annotation for", organism$genus, organism$species))
}
metadata(rsem_dds)$organism <- organism
metadata(rsem_dds)$template_git <- packageDescription("babsrnaseq")$git_last_commit

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

usethis::use_data(rsem_dds, overwrite=TRUE)

