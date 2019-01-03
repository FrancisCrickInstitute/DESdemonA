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

## *** Useful Functions
## Only use certain samples
applySubset <- function(dds, subs) {
  tmp <- dds[,subs]
  colData(tmp) <- droplevels(colData(tmp))
  tmp
}

## Set the design model
applyDesign <- function(dds, mdl) {
  design(dds) <- mdl
  DESeq(dds, test="Wald")
}
## Do likelihood ratio test, and classify in order of effect size
## Doesn't work for interactions, obviously
applyLRT <- function(dds, mdl) {
  design(dds) <- mdl$full
  dds <- DESeq(dds, test="LRT", full=mdl$full, reduced=mdl$reduced)
  cols <- resultsNames(dds)
  inBoth <- intersect(attr(terms(mdl$full), "term.labels"),
                     attr(terms(mdl$reduced), "term.labels"))
  testTerm <- setdiff(attr(terms(mdl$full), "term.labels"),
                     attr(terms(mdl$reduced), "term.labels"))
  if (length(inBoth)) {
    cols <- cols[!grepl(inBoth, cols)]
  }
  toOrd <- mcols(dds)[cols]
  toOrd[,1] <- 0
  ords <-  apply(toOrd, 1, order)
  str <- matrix(levels(dds[[testTerm]])[ords], nrow=nrow(ords), ncol=ncol(ords))
  mcols(dds)$class <- apply(str, 2, paste, collapse="<")
  dds
}
## apply contrast, and transfer across interesting mcols from the dds
applyResults <- function(dds, mcols=c("symbol"),...) {
  r <- results(dds, ...)
  r[mcols] <- mcols(dds)[mcols]
  if ("LRTPvalue" %in% names(mcols(dds))) {
    r$class <- mcols(dds)$class
    r$class[is.na(r$padj) | is.na(r$pvalue) | r$baseMean==0] <- NA
    r$log2FoldChange <- NULL
    r$lfcSE <- NULL
  }  else {
    r$class <- ifelse(r$log2FoldChange >0, "Up", "Down")
  }
  ind <- which(r$padj<metadata(r)$alpha)
  r$class[ind] <- paste0(r$class[ind], "*")
  r$class[is.na(r$pvalue)] <- "Outlier"
  r$class[is.na(r$padj)] <- "Low Count"
  r$class[r$baseMean==0] <- "Zero Count"
  r
}

set.seed(seed=params()$seed)


align_settings <- readLines("align.yml")
organism <- list(genus = grep("genus: .*", align_settings, value=TRUE),
                species = grep("species: .*", align_settings, value=TRUE))
organism <- map(organism, ~ gsub(".*: ", "", .))
organism$Gs <- paste0(toupper(substr(organism$genus, 1, 1)), tolower(substr(organism$species, 1, 1)))
organism$org <-  unique(grep(paste0("org\\.",organism$Gs),row.names(installed.packages()),  value=TRUE))
if (length(organism$org)!=1) {
  stop(paste("Can't find annotation for", organism$genus, organism$species))
}
library(organism$org, character.only=TRUE)

               
## ** Read in Data
design <- read.csv("data/design.csv", as.is=TRUE) 

samples <- read.xlsx("data/{{ASF_XL_FILE}}",sheetIndex=1, check.names=FALSE) %>%
  select_if( ~ length(unique(.x)) !=1) %>%
  dplyr::select(-`Sample TS/BA/Gel File Name`) %>%
  rename("Sample limsid"="sample", "Sample Name"="name") %>%
  separate(name, into=c("{{EXPERIMENTAL_CONDITIONS}}"), sep="_") %>%
  dplyr::select(-treatment, -num, -si)

txi <- tximport(sprintf("results/star/%s.genes.results", samples$id), type="rsem")
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

exprList <- lapply(ddsList, . %>% vst %>% assay)

top <- lapply(exprList,
             . %>% apply(1, sd) %>% order(decreasing=TRUE) %>% head(n=500)
             )

for (i in names(exprList)) {
  plotDat <- exprList[[i]][top[[i]],]
  plotDat <- plotDat-rowMeans(plotDat)
  vars <- all.vars(design(ddsList[[i]]))
  colDat <- as.data.frame(colData(ddsList[[i]])) %>%
    dplyr::select(!!vars) %>%
    as.data.frame
  rownames(colDat) <- colnames(plotDat) <- ddsList[[i]]$sample
  ph <- pheatmap(plotDat,
           annotation_col=colDat,
           show_rownames=FALSE, show_colnames=TRUE,
           silent=TRUE)
  
  poisd <- PoissonDistance(t(counts(ddsList[[i]], normalized=TRUE)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <-unite(colDat, col="exp_group", !!vars, remove=FALSE)[[1]]
  colnames(samplePoisDistMatrix) <- row.names(colDat)
  
  sq <- pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           annotation_col=colDat,
           silent=TRUE)$gtable
  vDevice(pdf, paste0("heatmap_", i))
  grid::grid.newpage(); grid::grid.draw(ph) 
  grid::grid.newpage();grid::grid.draw(sq)
  dev.off()#

  pc <- prcomp(t(exprList[[i]][top[[i]],]), scale=FALSE)
  percentVar <- pc$sdev^2/sum(pc$sdev^2)

  
  pc.df <- data.frame(PC1=pc$x[,1], PC2=pc$x[2,], colDat, sample=row.names(colDat))
  vDevice(pdf, paste0("pca_", i))
  gg <- ggplot(pc.df, aes(x=PC1, y=PC2, col=Treatment))  + geom_point() +
    geom_text_repel(aes(label=sample)) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
  print(gg)
  dev.off()
}








## ** Find differential genes

summaries <- list()
res <- list()


## *** LRT template
mdlList <- list("Pooled" = list(full = ~ {{CONDITION}},
                                reduced = ~ 1),
                "Normed" = list(full = ~ {{CONDITION}} + {{BATCH}},
                                reduced = ~ {{BATCH}})
ddsByDesign <- list(dds=ddsList, mdl=mdlList) %>%
  cross %>%
  map(lift(applyLRT))
names(ddsByDesign) <-list(names(ddsList), names(mdlList)) %>%
  cross %>%
  map(paste, collapse="_")
resultList <- lapply(ddsByDesign, applyResults, alpha=params("alpha"))

## *** Wald template
mdlList <- list("Pooled" = ~ {{CONDITION}}, "Normed" = ~ {{CONDITION}} + {{BATCH}})
ddsByDesign <- list(dds=ddsList, mdl=mdlList) %>%
  cross %>%
  map(lift(applyDesign))
names(ddsByDesign) <-list(names(ddsList), names(mdlList)) %>%
  cross %>%
  map(paste, collapse="\t")
contrList <- list("KO v WT"=c("CONDITION","KO","WT")
                  )
resultList <- list(dds=ddsByDesign, contrast=contrList) %>%
  cross %>%
  map(lift(applyResults, alpha=params("alpha")))
names(resultList) <- list(dds=ddsByDesign, contrast=contrList) %>%
  map(names) %>%
  cross %>%
  map(paste, collapse="\t")

## *** Summarise Results

tab <- map(resultList, ~table(Group=sub("\\*$","",.$class), Significant=grepl("\\*$",.$class)))
tab <- imap(tab,~ cbind(as.data.frame(.x), dataset=.y)) %>%
  bind_rows %>%
  spread(Significant, Freq) %>%
  rename(`FALSE`="not", `TRUE`="Significant") %>%
  mutate(Total=not+Significant) %>%
  dplyr::select(-not) %>%
  arrange(desc(Significant/Total))
tab <- tab %>%
  separate(dataset, into=c("Samples","Design","Comparison"), sep="\\t")

res$standard <- resultList
summaries$standard <- tab





## ** Generate results spreadsheet

################################################################
#### Enrichment Analysis
################################################################
vDevice(pdf, file="Enrichment Plots")
genes <- lapply(res, function(x) lapply(x, function(y) na.omit(entrez[row.names(y)[grepl("\\*", y$class)]])))
genes <- flatten(genes)
## Reactome
fit <- compareCluster(genes, fun="enrichPathway", organism="mouse", universe=na.omit(entrez))
pl <- dotplot(fit, showCategory=25)
pl + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6),
           axis.text.y = element_text(size=8))
## GO
fit<- compareCluster(genes, OrgDb = org.Mm.eg.db, universe=na.omit(entrez))
df <- fortify(fit)
levels(df$Description) <- gsub("RNA polymerase II", "RNAPii", levels(df$Description))
levels(df$Description) <- gsub(".*,","", levels(df$Description))
ggplot(df, aes(x = Cluster, y = Description, size = GeneRatio,  color = p.adjust)) +
  geom_point() + scale_color_gradient(low = "red",  high = "blue") +
  ylab("") + ggtitle("GO Enrichment") + theme_dose(6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6),
        axis.text.y = element_text(angle=45, size=8))
dev.off()

################################################################
#### Results Heatmaps
################################################################


wt_heatmap_dat <- assay(vst(dds_all))
wt_heatmap_cdat <- colData(dds_all)
row.names(wt_heatmap_cdat) <- wt_heatmap_cdat$name
colnames(wt_heatmap_dat) <- wt_heatmap_cdat$name


merge_hclust <- function(hclist) {
    d <- as.dendrogram(hclist[[1]])
    for (i in 2:length(hclist)) {
        d <- merge(d, as.dendrogram(hclist[[i]]))
    }
    as.hclust(d)
}


vDevice(pdf, file="Heatmap of results")
for (r in names(res)) {
  for (i in names(res[[r]])) {
    wt_heatmap_rdat <- data.frame(class=res[[r]][[i]]$class, row.names=row.names(wt_heatmap_dat))
    
    gene_ind <- grep("\\*", res[[r]][[i]]$class)
#    max_vst_per_gene <- apply(wt_heatmap_dat[gene_ind,], 1, max)
#    gene_ind <- gene_ind[max_vst_per_gene>quantile(max_vst_per_gene, .8)]
    group_sizes <- table(res[[r]][[i]]$class[gene_ind]) %>% sort(decreasing=TRUE)
    group_size <- group_sizes[res[[r]][[i]]$class[gene_ind]]
    gene_ind <- gene_ind[order(group_size, decreasing=TRUE)]
    gene_class <- factor(res[[r]][[i]]$class[gene_ind], levels=names(group_sizes))
    if (length(group_sizes)>6) {
      levels(gene_class)[cumsum(group_sizes)/sum(group_sizes) > 0.8] <- "Other*"
    }
    #  biggest_class <- sub("\\*","",res[[r]][[i]]$class[gene_ind][1])
    biggest_class <- sub("\\*", "",  levels(gene_class)[1])
    if (grepl("<", biggest_class)) {
      sample_order <- strsplit(biggest_class, "<")[[1]]
      wt_heatmap_cdat$tissue <- factor(as.character(wt_heatmap_cdat$tissue), levels=sample_order)
      sample_ind <- order(match(wt_heatmap_cdat$tissue, sample_order), wt_heatmap_cdat$experiment)
    } else {
      sample_ind <- TRUE
    }

    plList <- tapply(gene_ind, gene_class, function(x)
      pheatmap(wt_heatmap_dat[x,sample_ind,drop=FALSE],  cluster_cols=FALSE, silent=TRUE)$tree_row)
    
    row_hclust <- merge_hclust(plList)
    pl <- pheatmap(wt_heatmap_dat[gene_ind,sample_ind],
                   annotation_col=as.data.frame(wt_heatmap_cdat)[sample_ind, c( "tissue"),drop=FALSE],
                   #                                  annotation_row=wt_heatmap_rdat[gene_ind,,drop=FALSE],
                   scale="row",
                   cluster_rows=row_hclust, # cluster_cols=!(grepl("<", biggest_class)),
                   show_rownames=FALSE, show_colnames=FALSE,
                   silent=TRUE,
                   main=i)$gtable
    ind <- pl$layout$name!="row_tree"
    pl$layout <- pl$layout[ind, , drop = FALSE]
    pl$grobs <- pl$grobs[ind]
    pl <- gtable::gtable_trim(pl)
    grid::grid.draw(pl)
    grid::grid.newpage()
  }
}
dev.off()

## *** Output filtered results
si <- session_info()
for (analysis in names(res)) {
  summaries <- do.call(cbind, lapply(res[[analysis]], function(x) table(x$class)))
  colnames(summaries) <- names(res[[analysis]])
  ind <- rownames(summaries) %in% c("Low Count", "Zero Count", "Outlier")
  summaries <- rbind(summaries[ind,,drop=FALSE], summaries[!ind,,drop=FALSE])
  dframe <- t(as.data.frame(params()))
  fname <- vDevice(write.xlsx, glue("differential_genelists_{analysis}.xls"),
                  x=as.data.frame(dframe), sheetName="Parameters", col.names=FALSE)
  ## Design
  write.xlsx(colData(dds_all), file=fname, sheetName="Design", col.names=TRUE,row.names=FALSE,  append=TRUE)
  write.xlsx(summaries, file=fname, sheetName="Class Sizes", col.names=TRUE,row.names=TRUE,  append=TRUE)
  ## Differential gene-lists
  for (i in names(res[[analysis]])) {
    dframe <- as.data.frame(res[[analysis]][[i]]) %>%
      rownames_to_column("id") %>%
      filter(padj<params()$alpha) %>%
      arrange(desc(stat))
    write.xlsx2(dframe, file=fname, sheetName=i, append=TRUE, row.names=FALSE)
  }
  write.xlsx2(as.data.frame(si$packages), file=xl, sheetName="R Packages", row.names=FALSE, append=TRUE)
  write.xlsx2(data.frame(setting = names(si$platform), value = unlist(si$platform), stringsAsFactors = FALSE), file=xl, sheetName="R Settings", row.names=FALSE, append=TRUE)
}
vDevice(saveRDS, "sessionInfo.rds", sessionInfo())

for (i in names(res)) {
  for (j in names(res[[i]])) { 
    vDevice(txt, file=paste(i,j, sep="-"), x=as.data.frame(res[[i]][[j]]) %>% dplyr::select(stat, symbol, class))
  }
}
