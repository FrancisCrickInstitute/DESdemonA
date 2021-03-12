defaults <- list(
  nfcore      = "results",                                                  ## directory that contains nfcore output
  metadata    = system.file("extdata/metadata.xlsx", package="babsrnaseq"), ## Path to  spreadsheet containing covariates
  file_col    = "filename",                                                 ## Column in that spreadsheet that identifies the basename of count file for each sample
  counts      = "results/star_rsem",                                        ## directory that contains the quantified gene counts
  org_package = NULL                                                        ## NULL will derive org.*.eg.db from nfcore info, otherwise specifcy string identifying annotation package of relevant species.
)


library(readxl)
library(tidyverse)
library(tximport)
library(DESeq2)
library(DESdemonA)
library(R.utils)


args <- R.utils::cmdArgs(args=defaults)

sample_sheet <- readxl::read_excel(args$metadata), sheet=1)

if (is.null(sample_sheet[[args$file_col]])) {
  sample_sheet[[args$file_col]]  <- paste0(sample_sheet$LIMSID, "_R1") # that annoying nfcore '_R1' suffix may disappear
}  
  
sample_sheet[sapply(sample_sheet, is.character)] <- lapply(sample_sheet[sapply(sample_sheet, is.character)], 
                                       as.factor)

txi <- tximport(
  file.path(
    args$counts,
    paste0(as.character(sample_sheet[[args$file_col]]), ".genes.results")
   ),
  type="rsem")
txi$length[txi$length==0] <- 1

init_dds <- DESeqDataSetFromTximport(txi, sample_sheet, ~ 1) # ok to have fixed trivial design - it'll get swapped out
ind <- rowSums(counts(init_dds)!=0) > 0
init_dds <- init_dds[ind,]

txi[c("abundance", "counts", "length")] <- map(txi[c("abundance", "counts", "length")], ~ .[ind,])


if (is.null(args$org_package)) {
  run_info <- readLines(file.path(args$nfcore, "pipeline_info", "pipeline_report.txt"))
  gtf_file <- sub(".*/", "", grep(" - gtf:", run_info, value=TRUE))
  organism <- list(binomial=strsplit(gtf_file, ".", fixed=TRUE)[[1]][1])
  organism <- list(
    genus   = sub("_.*", "", organism$binomial),
    species = sub(".*_", "", organism$binomial)
  )
  organism$Gs <- paste0(toupper(substr(organism$genus, 1, 1)),
                       tolower(substr(organism$species, 1, 1)))
  organism$org <-  unique(grep(paste0("org\\.",organism$Gs),row.names(installed.packages()),  value=TRUE))
  if (length(organism$org)!=1) {
    stop(paste("Can't find annotation for", organism$genus, organism$species))
  }
} else {
  organism <- list(org=args$org_package)
}
metadata(init_dds)$organism <- organism
metadata(init_dds)$template_git <- packageDescription("babsrnaseq")$git_last_commit

library(organism$org, character.only=TRUE)

mcols(init_dds)$symbol <- mapIds(
  eval(parse(text = organism$org)), 
  keys=row.names(init_dds),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first")[row.names(init_dds)]

mcols(init_dds)$entrez <- mapIds(
  eval(parse(text = organism$org)),
  keys=row.names(init_dds),
  column="ENTREZID",
  keytype="ENSEMBL",
  multiVals="first")[row.names(init_dds)]

usethis::use_data(init_dds, overwrite=TRUE)

