library(tidyverse)
library(tximport)
library(DESeq2)
library(openxlsx)
devtools::load_all()

## samples <- read.xlsx(system.file("extdata/{{ASF_XL_FILE}}", package="babsRNASeq"),
##                     sheetIndex=1, check.names=FALSE) %>%
##   select_if( ~ length(unique(.x)) !=1) %>%
##   dplyr::select(-`Sample TS/BA/Gel File Name`) %>%
##   rename("Sample limsid"="sample", "Sample Name"="name") %>%
##   separate(name, into=c("{{EXPERIMENTAL_CONDITIONS}}"), sep="_") %>%
##   dplyr::select(-treatment, -num, -si)

sample_sheet <- read.xlsx(system.file("extdata/{{ASF_XL_FILE}}", package="babsRNASeq"),
                         startRow=2) %>%
  dplyr::select(starts_with("Sample.")) %>%
  rename_all(~ sub("^Sample\\.", "", .)) %>%
  dplyr::select(-Genotype) %>%
  dplyr::mutate(
    Lympho=case_when(
      Treatment=="baseline" ~ "baseline",
      Treatment=="KOL" ~ "KO",
      TRUE ~ "WT"),
    Treatment=case_when(
      Treatment=="wtL/Cx5" ~ "Cx5",
      Treatment=="wtL/IgG" ~ "IgG",
      TRUE ~ "ctrl"),
    rsem=file.path(system.file("extdata/rsem", package="babsRNASeq"), paste0(limsid, ".genes.results")),
    Mouse=Replicate.Group
  ) %>%
  dplyr::select(-Replicate.Group)


txi <- tximport(sample_sheet$rsem, type="rsem")
txi$length[txi$length==0] <- 1
rsem_dds <- DESeqDataSetFromTximport(txi, sample_sheet, ~ Group + Mouse)
ind <- rowSums(counts(rsem_dds)!=0) > 0
rsem_dds <- rsem_dds[ind,]
txi[c("abundance", "counts", "length")] <- map(txi[c("abundance", "counts", "length")], ~ .[ind,])


align_settings <- readLines("align.yml")
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
mcols(rsem_dds)symbol <- mapIds(
  eval(parse(text = organism$org)) 
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


metadata(rsem_dds)$labels <- c("Treatment", "Lympho", "Group", "Mouse")


usethis::use_data(rsem_dds, overwrite=TRUE)


