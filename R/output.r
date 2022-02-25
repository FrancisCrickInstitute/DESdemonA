##' Xlsx reporting of results
##'
##' Store all the differential gene-lists and supporting materials
##' in a multi-worksheet spreadsheet
##' @title XLSX report of results
##' @param ddsList A list of [DESeq2::DESeqDataSet-class()] objects generated by DESdemonA containing results in the mcols slot
##' @param param The parameter object used to generate the results
##' @param dir Directory to store the results in
##' @return A list of file paths to the excel files
##' @author Gavin Kelly
#' @export
write_results <- function(ddsList, param, dir=".", assays=NULL) {
  si <- session_info()
  crick_colours <-list(
    primary=list(red="#e3001a",yellow="#ffe608",blue="#4066aa",green="#7ab51d",purple="#bb90bd"),
    secondary=list(pink="#fadbd9", yellow="#fff6a7", green="#adcf82", orange="#ffe7ab", blue="#bee2e6"),
    spare=list(blue="#95ceed"))
  hs1 <- createStyle(fgFill = crick_colours$secondary$orange, textDecoration = "italic",
                    border = "Bottom")
  hs2 <- createStyle(fgFill = crick_colours$secondary$blue, textDecoration = "italic",
                    border = "Bottom")
  summaries <- map_depth(ddsList, 3, DESdemonA::summarise_results)
  out <- lapply(ddsList, function(x) "")
  for (dataset in names(ddsList)) {
    wb <- openxlsx::createWorkbook(title=param$get("title"),
                                  creator="Gavin Kelly")
    tmp <- param$describe()
    dframe <- data.frame(id=names(tmp), description=unlist(tmp))
    sn <- "Parameters"
    addWorksheet(wb, sn)
    writeData(wb, sn, dframe,rowNames=FALSE, colNames=TRUE)
    ## Design
    samples_used <- as.data.frame(colData(ddsList[[dataset]][[1]][[1]]))
    sn <- "Design"
    addWorksheet(wb, sn)
    writeData(wb, sn, samples_used, headerStyle=hs2)
    sn <- "Class Sizes"
    addWorksheet(wb, sn)
    dframe <- DESdemonA::rbind_summary(
      summaries[[dataset]],
      levels=c("Design","Comparison")
    )
    writeData(wb, sn, dframe, headerStyle=hs2)
    ## Differential gene-lists
    comparison_name_lookup <- list()
    addWorksheet(wb, "Comparison Key")
    for (design_ind in 1:length(ddsList[[dataset]])) {
      for (contrast_name in names(ddsList[[dataset]][[design_ind]])) {
        dframe <- as.data.frame(mcols(ddsList[[dataset]][[design_ind]][[contrast_name]])$results)
        for (assay_name in assays) {
          if (!assay_name %in% assayNames(ddsList[[dataset]][[design_ind]][[contrast_name]])) {
            warning(assay_name, " not an assay, so not added to output")
            next
          }
          this_assay <- assay(ddsList[[dataset]][[design_ind]][[contrast_name]], assay_name)
          if (length(assays)>1) {
            names(this_assay) <- paste(this_assay, names(this_assay), sep="_")
          }
          dframe <- cbind(dframe, this_assay)
        }
        dframe <- dframe %>%
          tibble::rownames_to_column("id") %>%
#          dplyr::filter(padj<param$get("alpha")) %>%
          dplyr::arrange(desc(abs(shrunkLFC))) %>%
          dplyr::select(-padj)
        if (length(ddsList[[dataset]])==1) {
          sn <- contrast_name
        } else {
          sn <- paste0(contrast_name, ", ", names(ddsList[[dataset]])[design_ind])
        }
        if (nchar(sn)>31) {
          alpha_key <- DESdemonA:::to_letter(length(comparison_name_lookup)+1)
          comparison_name_lookup[[alpha_key]] <- sn
          sn <- alpha_key
        }
        addWorksheet(wb, sn, tabColour=crick_colours$secondary[[design_ind]])
        writeData(wb, sn, dframe, headerStyle=hs1, withFilter=TRUE)
        groupRows(wb, sn, rows=which(!grepl("\\*$", dframe$class))+1, hidden=TRUE)
        filtCol <- match("class", names(dframe))
        if (!is.na(filtCol)) {
          filt_string <- sprintf(
            '><filterColumn colId="%s"><customFilters><customFilter val="*~*"/></customFilters></filterColumn></autoFilter>',
            filtCol-1
          )
          sheet_n <- match(sn, names(wb))
          wb$worksheets[[sheet_n]]$autoFilter <- sub("/>$", filt_string, wb$worksheets[[sheet_n]]$autoFilter)
        }
      }
    }
    if (length(comparison_name_lookup)==0) {
      removeWorksheet(wb, "Comparison Key")
    } else {
      writeData(wb, "Comparison Key", data.frame(Key=names(comparison_name_lookup),
                                                 Comparison=unlist(comparison_name_lookup)
                                                 )
                )
    }
    ## sn <- "GO terms"
    ## addWorksheet(wb, sn)
    ## writeData(wb, sn, go_df, headerStyle=hs1)
    sn <- "R Packages"
    addWorksheet(wb, sn)
    writeData(wb, sn, as.data.frame(si$packages), headerStyle=hs2)
    sn <- "R Details"
    addWorksheet(wb, sn)
    writeData(wb, sn, data.frame(setting = names(si$platform),
                                 value = unlist(si$platform),
                                 stringsAsFactors = FALSE),
              headerStyle=hs2)
    out[[dataset]] <- file.path(dir, paste0("differential_", param$get("spec"), "_", dataset, ".xlsx"))
    (saveWorkbook(wb, out[[dataset]], overwrite=TRUE))
  }
  out
}

##' Store results as text files
##'
##' Save unfiltered versions of the results in text files
##' @title Store results as text files
##' @param ddsList A list of [DESeq2::DESeqDataSet-class()] objects generated by DESdemonA containing results in the [S4Vectors::mcols()] slot
##' @param dir Directory to store the results in
##' @return 
##' @author Gavin Kelly
#' @export
write_all_results <- function(ddsList, dir=".") {
  for (i in names(ddsList)) {
    for (j in names(ddsList[[i]])) { 
      for (k in names(ddsList[[i]][[j]])) { 
        readr::write_excel_csv(
          path=file.path(dir, sprintf("allgenes_%s_%s_%s.csv", i, j, k)),
          x=as.data.frame(mcols(ddsList[[i]][[j]][[k]])$results) %>% dplyr::select(log2FoldChange, stat, symbol, class))
      }
    }
  }
}



to_letter <- function(i, so_far="") {
  if (i<27)
    paste0(LETTERS[i], so_far)
  else
    to_letter((i %/% 26), LETTERS[((i-1) %% 26) + 1])
}


##' Produce Table 1 of sample metadata
##'
##' Make a GT object representing all the original metadata and, for
##' each sample subset, the sample's presence and any extra columns introduced by a
##' transform.
##' 
##' @title Table One
##' @param dds The reference dds object
##' @param ddsList The DESdemonA list of datasets derived from the reference object
##' @return a GT object
##' @author Gavin Kelly
#' @export
table1 <- function(dds, ddsList) {
  samp <- as.data.frame(colData(dds))
  # Extra columns introduced by datasets
  extra_meta <- lapply(
    names(ddsList),
    function(x) merge(
      setNames(data.frame(ifelse(colnames(dds) %in% colnames(ddsList[[x]]), "✓", ""), row.names=row.names(colData(dds))), x),
      as.data.frame(colData(ddsList[[x]])[, setdiff(names(colData(ddsList[[x]])), names(colData(dds)))]),
      by=0, all.x=TRUE
      )[,-1,drop=FALSE]
  )
  
  extra_meta <- extra_meta[sapply(extra_meta, function(x) ncol(x)!=0 & nrow(x)!=0)]
  # Generate tab_spanner parameters for each datasets extra columns
  ts_extra <- lapply(names(extra_meta),
                    function(x) list(name=x,
                              cols=1:ncol(extra_meta[[x]]))
                    )
  # Cummulatively offset the columns to be spanned
  ts_extra <- Reduce(function(ts, x) {
    x$cols <- x$cols+max(ts$cols)
    x},
    ts_extra,
    accumulate=TRUE)
  
  df <- cbind(samp, do.call(cbind, unname(extra_meta)))
  clabel <- setNames(as.list(names(df)), 1:ncol(df))
  names(df) <- names(clabel)
             
  gt(data=df,
     caption="Sample annotation") %>%
    cols_label(.list=clabel) %>%
    tab_spanner(label="Metadata",
                columns=seq_along(samp)) %>%
    Reduce(f=function(gti, x) tab_spanner(gti, label=x$name,columns=x$cols + ncol(samp)), x=ts_extra, init=.) %>%
    DESdemonA::tab_link_caption()
}

##' Generate text files require for Biologic
##'
##' The visualisation part of 'Biologic' requires a set of text files
##' which represent the differential results and the models and
##' contrasts that were used to derive the results.  This function
##' generates them from a previous run of DESdemonA.
##' 
##' @title Export Biologic files
##' @param result_object Path to the rds of the DESdemonA result object
##' @param path Where to write the files for Biologic to read
##' @return nothing
##' @author Gavin Kelly
#' @export
export_biologic <- function(result_object, path) {
  obj <- readRDS(result_object)
  comparison_table <- list()
  definition_table <- list()
  models_list <- list()
  coldata_list <- list()
  rsem_raw <- NULL
  rsem_norm <- NULL
  for (subset_name in names(obj)) {
    this_subset <- obj[[subset_name]]
    first_dds <- this_subset[[1]][[1]]
    coldata_list <- c(coldata_list, list(as.data.frame(colData(first_dds))))
    this_raw_counts <- counts(first_dds, norm=FALSE)
    this_norm_counts <- counts(first_dds, norm=TRUE)
    novel_samples <- setdiff(colnames(this_raw_counts), colnames(rsem_raw))
    if (length(novel_samples)!=0) {
      rsem_raw <- cbind(rsem_raw, this_raw_counts[, novel_samples])
      rsem_norm <- cbind(rsem_norm, this_norm_counts[, novel_samples])
    }
    for (model_name in names(this_subset)) {
      this_model <- this_subset[[model_name]]
      for (comparison_name in names(obj[[subset_name]][[model_name]])) {
        dds <- this_model[[comparison_name]]
        ## Results file
        results_frame <- as.data.frame(mcols(dds)$results)
        write.table(results_frame,
#                    file=file.path(path, paste0(name_sanitizer(paste(subset_name, model_name, comparison_name, sep="_")), ".txt")),
                    file=file.path(path, paste0(name_sanitizer(paste(comparison_name, sep="_")), ".txt")),
                    sep="\t", col.names=NA
                    )
        ## Model file
        is_lrt <- rlang::is_formula(metadata(dds)$models$comparisons[[1]])
        comparison_row <- data.frame(
          comparison=name_sanitizer(metadata(dds)$dmc$comparison),
          test=ifelse(is_lrt, "LRT", "Wald"),
          type=ifelse(is_lrt, "LRT", "DGE"),
          model=capture.output(dput(metadata(dds)$model$design)),
          reducedModel=ifelse(is_lrt, capture.output(dput(metadata(dds)$model$comparisons[[1]])), "")
        )
        comparison_table <- c(comparison_table, list(comparison_row))
        ## Definition file
        models_list <- c(models_list,list(metadata(dds)$model$lm$model[,-1,drop=FALSE]))
      }
    }
  }

  #  colnames(rsem_raw) <- coldata_frame(
  write.table(rsem_raw,
              file=file.path(path, "rsem_raw.txt"),
              sep="\t", row.names=TRUE, col.names=NA
              )

  write.table(rsem_norm,
              file=file.path(path, "rsem_norm.txt"),
              sep="\t", row.names=TRUE, col.names=NA
              )
  
  comparison_frame <- do.call(rbind, comparison_table)
  write.table(comparison_frame,
              file=file.path(path, "design.model.file.txt"),
              sep="\t", row.names=FALSE
              )
  
  sample_id_list <- lapply(coldata_list, function(df) data.frame(sampleID=df[[1]], sample.id=df$sample_label))
  definition_frame <- do.call(cbind, sample_id_list[!duplicated(sample_id_list)])
  coldata_list <- lapply(coldata_list, function(df) {
    df[!grepl("^\\.PCA\\.PC", names(df))]
  })
  coldata_frame <- do.call(cbind, coldata_list[!duplicated(coldata_list)])
  n_unique <- sapply(coldata_frame, function(x) length(unique(x)))
  interesting <- (n_unique!=1 & n_unique!=nrow(coldata_frame)) | sapply(coldata_frame, is.numeric)
  definition_frame <- cbind(data.frame(sample.groups= do.call(paste, coldata_frame[,interesting,drop=FALSE])),
                           definition_frame)
  comparisons_for_def <- as.data.frame(matrix("", nrow=nrow(definition_frame), ncol=nrow(comparison_frame),
                                             dimnames=list(row.names(definition_frame), name_sanitizer(comparison_frame$comparison))
                                             ))
  definition_frame <- cbind(definition_frame, comparisons_for_def)
  write.table(definition_frame,
              file=file.path(path, "design.dge.lrt.definition.file.txt"),
              sep="\t", row.names=FALSE
              )

}

name_sanitizer <- function(str) {  gsub("[^a-zA-Z0-9\\._]", "_", str) }

