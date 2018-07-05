## * Code

## ** Setup

package_list <- c("xlsx", "dplyr", "tidyr","readr",
             "purrr", "pheatmap", "devtools","PoiClaClu","rtracklayer",
             "tximport", "ggplot2", "ggrepel", "grid","RColorBrewer",
             "xlsx", "jsonlite", "DESeq2")
all(vapply(package_list, require, logical(1), character=TRUE, quietly=TRUE)) ||
  stop("Not all packages present")

if (!require(crick.kellyg, quietly=TRUE)) {
  if (Sys.getenv("my_r_package")!="") {
    load_all(Sys.getenv("my_r_package"))
  } else {
    install_github("macroscian/R-Package")
  }
}

params <- params_init(
  title="RNASeq Analysis",
  script = "analyse.r",
  seed=1,
  alpha=0.05,
  lfcThreshold=0
)


set.seed(seed=params()$seed)
dirs <- derivedDirs(publish=NA,
                   subResults=params()$fVersion)


align_settings <- readLines("align.yml")
organism <- list(genus = grep("genus: .*", align_settings, value=TRUE),
                species = grep("species: .*", align_settings, value=TRUE))
organism <- map(organism, ~ gsub(".*: ", "", .))
organism$Gs <- paste0(toupper(substr(organism$genus, 1, 1)), tolower(substr(organism$species, 1, 1)))
organism$org <-  grep(paste0("org\\.",organism$Gs),row.names(installed.packages()),  value=TRUE)
if (length(organism$org)!=1) {
  stop(paste("Can't find annotation for", organism$genus, organism$species))
}
library(organism$org, character.only=TRUE)

               
## ** Read in Data
design <- read.csv("design.csv", as.is=TRUE) %>%
  mutate(id=sample, sample=gsub(".*(PER.*)_S.*", "\\1", file1)) %>%
  dplyr::select(sample, id)

samples <- read.table("data/sample_sheet.txt", header=TRUE, sep="\t", as.is=TRUE) %>%
  left_join(design) %>%
  separate(name, into=c("Irrelevant", "Condition","Replicate"), sep="\\.") %>%
  dplyr::select(-Irrelevant)
txi <- tximport(sprintf("results_fp0/star/%s.genes.results", samples$id), type="none", txIn = FALSE, geneIdCol = "gene_id",  abundanceCol = "TPM", countsCol = "expected_count", lengthCol = "effective_length", importer=read_tsv)
txi$length[txi$length==0] <- 1
dat <- DESeqDataSetFromTximport(txi, samples, ~ Time)
ind <- rowSums(counts(dat)!=0) > 0
dat <- dat[ind,]
txi[c("abundance", "counts", "length")] <- map(txi[c("abundance", "counts", "length")], ~ .[ind,])


symbol <- mapIds(get(organism$org),
                keys=row.names(dat),
                column="SYMBOL",
                keytype="ENSEMBL",
                multiVals="first")

mcols(dat)$symbol <- symbol[row.names(dat)]

ddsList <- list(all=DESeq(dat))


## ** QC Visualisation
## tmp <- tmp[,tmp$sample!="RM08"]
## tmp$VD <- droplevels(tmp$VD)
## tmp$D3 <- droplevels(tmp$D3)
## ddsList$omit_11_8=DESeq(tmp)

exprList <- lapply(ddsList, . %>% vst %>% assay)

top <- lapply(exprList,
             . %>% apply(1, sd) %>% order(decreasing=TRUE) %>% head(n=500)
             )

for (i in names(exprList)) {
  plotDat <- exprList[[i]][top[[i]],]
  plotDat <- plotDat-rowMeans(plotDat)
  colDat <- as.data.frame(colData(ddsList[[i]])) %>%
    dplyr::select(Time, Replicate, Concentration)
  rownames(colDat) <- colnames(plotDat) <- ddsList[[i]]$sample
  ph <- pheatmap(plotDat,
           annotation_col=colDat,
           show_rownames=FALSE, show_colnames=TRUE,
           silent=TRUE);dev.off()
  pl <- ph$gtable
  clst <- cutree(ph$tree_row, 4)
  inds <- map(1:4, ~ which(clst==.))
  inds <- inds[sapply(inds, length)!=1]
  
  pl_sub <- map(inds, ~ pheatmap(plotDat[.,],
           annotation_col=colDat,
           show_rownames=FALSE, show_colnames=TRUE,
           silent=TRUE)$gtable);dev.off()

  poisd <- PoissonDistance(t(counts(ddsList[[i]], normalized=TRUE)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <-
    with(colData(ddsList[[i]]), paste(Time, Replicate,  sep=" - "))
  colnames(samplePoisDistMatrix) <- row.names(colDat)
  
  sq <- pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           annotation_col=colDat,
           silent=TRUE)$gtable;dev.off()
  vDevice(pdf, paste0("heatmap_", i))
  grid::grid.newpage(); grid::grid.draw(pl) 
  grid::grid.newpage();grid::grid.draw(sq)
  lapply(pl_sub, function(x) {grid::grid.newpage(); grid::grid.draw(x)})
  dev.off()#

  pc <- prcomp(t(exprList[[i]][top[[i]],]), scale=FALSE)
  percentVar <- pc$sdev^2/sum(pc$sdev^2)
  loadings <- data.frame(symbol=mcols(ddsList[[1]])$symbol[top[[i]]],pc$rotation)
  vDevice(txt, paste0("loadings_", i), loadings)

  
  colDat <- as.data.frame(colData(ddsList[[i]]))
  pc.df <- data.frame(PC1=pc$x[,1], PC2=pc$x[2,], colDat)
  vDevice(pdf, paste0("pca_", i))
  gg <- ggplot(pc.df, aes(x=PC1, y=PC2, col=Time))  + geom_point() +
    geom_text_repel(aes(label=sample)) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    coord_fixed()
  print(gg)
  print(gg+aes(y=Concentration))
  dev.off()
}





tpm <- txi$abundance
colnames(tpm) <- paste(ddsList$all$Time, ddsList$all$Replicate, sep="_")
vDevice(txt, "tpm", tpm)
nrc <- counts(ddsList[[1]], normalized=TRUE)
colnames(nrc) <- colnames(tpm)
vDevice(txt, "Normalised_counts", nrc)


## ** Find differential genes

summaries <- list()
res <- list()

## *** First analysis approach
analysis <- params("analysis" = "Across_time")
mdlList <- list(dat=list(AP=~AP, VD=~VD, D3=~D3),
               omit_11=list(AP=~AP, D3=~D3))
desList <- map(ddsList, ~DESeq(., reduced=~1, test="LRT"))

resultList <- map(desList, results)
coefList <- map(desList, coef)

ord <- function(m) {
  m[,1] <- 0
  id <- sub("[A-Z].*?_(.*)_vs.*", "\\1", colnames(m))
  id[1] <- sub(".*_vs_", "", colnames(m)[2])
  apply(m, 1, function(x) paste(id[order(x)], collapse="<"))
}

resultList <- imap(resultList, function(r,i) {
  r$symbol <- mcols(desList[[i]])$symbol
  r$group <- ord(coefList[[i]])
  r$tpm5 <- rowSums(tpm>5)
  r[!(names(r) %in% c("log2FoldChange", "lfcSE"))] 
})

tab <- map(resultList, ~table(Group=.$group, Significant=.$padj<params()$alpha))
tab <- imap(tab,~ cbind(as.data.frame(.x), dataset=.y)) %>%
  bind_rows %>%
  spread(Significant, Freq) %>%
  rename(`FALSE`="not", `TRUE`="Significant") %>%
  mutate(Total=not+Significant) %>%
  dplyr::select(-not) %>%
  arrange(desc(Significant/Total))



res[[analysis]] <- resultList
summaries[[analysis]] <- tab


vDevice(pdf,"barchart_behaviours")
tab2 <- map(resultList, ~table(Group=.$group, Significant=.$padj<params()$alpha)) %>%
  imap(~ cbind(as.data.frame(.x), dataset=.y)) %>%
  bind_rows
tab2$Group <- factor(tab2$Group, levels=tab$Group)
ggplot(tab2 , aes(x=Group, y=Freq, fill=Significant)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()




## ** Generate results spreadsheet


## *** Output filtered results
vDevice(saveRDS, "sessionInfo.rds", sessionInfo())
si <- session_info()
for (analysis in names(res)) {
  xl <- vDevice(file=paste0(analysis, ".xlsx"))
  dframe <- t(as.data.frame(params()))
  write.xlsx(as.data.frame(dframe), file=xl, sheetName="Parameters", col.names=FALSE)
  #dframe <- do.call(rbind, lapply(ddsList, compose(as.data.frame, colData)))
  dframe <- do.call(rbind, imap(ddsList, ~cbind(Dataset=.y, as.data.frame(colData(.x)))))
  ## Design
  write.xlsx(dframe, file=xl, sheetName="Design", col.names=TRUE,row.names=FALSE,  append=TRUE)
  dframe <- data.frame(summaries[[analysis]], check.names=FALSE)
  write.xlsx(dframe, file=xl, sheetName="Summary", col.names=TRUE,row.names=TRUE,  append=TRUE)
  ## Differential gene-lists
  for (i in names(res[[analysis]])) {
    dframe <- subset(res[[analysis]][[i]], padj<params()$alpha)#[c("stat", "pvalue", "padj",   "symbol", "group")]
    dframe <- dframe[order(dframe$stat),]
    write.xlsx2(dframe, file=xl, sheetName=i, append=TRUE)
  }
  write.xlsx2(as.data.frame(si$packages), file=xl, sheetName="R Packages", row.names=FALSE, append=TRUE)
  write.xlsx2(data.frame(setting = names(si$platform), value = unlist(si$platform), stringsAsFactors = FALSE), file=xl, sheetName="R Settings", row.names=FALSE, append=TRUE)
}

for (i in names(res)) {
  for (j in names(res[[i]])) { 
    vDevice(txt, file=paste(i,j, sep="-"), x=as.data.frame(res[[i]][[j]]) %>% dplyr::select(stat, symbol, group, tpm5))
  }
}
