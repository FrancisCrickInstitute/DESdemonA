##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Geneset enrichment analysis
##' @param ddsList 
##' @param fun 
##' @param showCategory 
##' @param max_width 
##' @return 
##' @author Gavin Kelly
##' @export
enrichment <- function(dds, fun="gseGO") {
    ord <- order(mcols(dds)$results$shrunkLFC, decreasing=TRUE)
    ord <- ord[!is.na(mcols(dds)$entrez)]
    genelist <- setNames(mcols(dds)$results$shrunkLFC, mcols(dds)$entrez)[ord]
    if (fun=="gseGO") {
      res <- gseGO(
        genelist,
        ont = "MF",
        OrgDb=get(metadata(dds)$organism$org))
    }
    if (fun=="gsePathway") {
      reactome_org <- DESdemonA:::org2reactome(metadata(dds)$organism$org)
      res <- gsePathway(
        genelist,
        organism = reactome_org)
    }
    if (nrow(res)>0) {
      res
    } else {
      NULL
    }
}




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Geneset over-representation analysis
##' @param ddsList 
##' @param fun 
##' @param showCategory 
##' @param max_width 
##' @return 
##' @author Gavin Kelly
##' @export
over_representation <- function(ddsList, fun, showCategory, max_width=30) {
  genes <- lapply(ddsList, function(dds) {
    res <- mcols(dds)$results
    na.omit(res$entrez[grepl("\\*", res$class)])
  })
  genes <- genes[sapply(genes, length)!=0]
  if (length(genes)<1) {
    return(NULL)
  }
  if (fun=="enrichGO") {
    reactome <- try(eval(substitute(compareCluster(genes, fun=fun, OrgDb=metadata(ddsList[[1]])$organism$org, universe=na.omit(metadata(ddsList[[1]])$entrez)), list(fun=fun))), silent=TRUE)
  } else {
    reactome_org <- DESdemonA:::org2reactome(metadata(ddsList[[1]])$organism$org)
    reactome <- try(eval(substitute(compareCluster(genes, fun=fun, organism=reactome_org, universe=na.omit(metadata(ddsList[[1]])$entrez)), list(fun=fun))), silent=TRUE)
  }
  if (inherits(reactome,"try-error")) {
    return(NULL)
  }
  enrich_table <- as.data.frame(reactome)[c("Cluster", "ID", "Description","GeneRatio","BgRatio")]
  reactome@compareClusterResult$Description <- ifelse(
    nchar(reactome@compareClusterResult$Description)>max_width,
    reactome@compareClusterResult$ID,
    reactome@compareClusterResult$Description)
  pl <- dotplot(reactome, showCategory=showCategory) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6),
          axis.text.y = element_text(size=8)
          )
  if (nrow(enrich_table)>0) {
    list(plot=pl, table=enrich_table)
  } else {
    NULL
  }
}


org2reactome <- function(org) {
  orgs <- c(anopheles = "org.Ag.eg.db", arabidopsis = "org.At.tair.db", 
           bovine = "org.Bt.eg.db", canine = "org.Cf.eg.db", celegans = "org.Ce.eg.db", 
           chicken = "org.Gg.eg.db", chimp = "org.Pt.eg.db", coelicolor = "org.Sco.eg.db", 
           ecolik12 = "org.EcK12.eg.db", ecsakai = "org.EcSakai.eg.db", 
           fly = "org.Dm.eg.db", gondii = "org.Tgondii.eg.db", human = "org.Hs.eg.db", 
           malaria = "org.Pf.plasmo.db", mouse = "org.Mm.eg.db", 
           pig = "org.Ss.eg.db", rat = "org.Rn.eg.db", rhesus = "org.Mmu.eg.db", 
           xenopus = "org.Xl.eg.db", yeast = "org.Sc.sgd.db", zebrafish = "org.Dr.eg.db"
           )
  names(orgs[orgs==org])
}
  
