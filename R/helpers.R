#' Merge Nested Data Frames
#'
#' This function merges nested data frames in a list into a single data frame. It also adds new columns to the merged data frame
#' indicating the path of each original data frame in the list structure.
#'
#' @param nested_list A list containing nested data frames.
#' @param parent_names A character vector specifying the names of the parent elements in the list structure. Default is an empty vector.
#' @param columnsToDrop A character vector specifying the names of the columns to drop from the data frames. Default is an empty vector.
#' @param forceColumns A character vector specifying the names of the columns to force in the data frames. If a data frame does not
#' have a column specified in forceColumns, the column will be added with NA values. Default is an empty vector.
#' @return A data frame resulting from merging the nested data frames in nested_list.
#' @examples
#' \dontrun{
#' merge_nested_df(nested_list, parent_names = c(), columnsToDrop = c(), forceColumns = c())
#' }
#' @export
#' @importFrom stats setNames
merge_nested_df <- function(nested_list, parent_names = c(), columnsToDrop = c(), forceColumns = c()) {
  # This will hold merged data frames
  merged <- do.call("rbind", lapply(names(nested_list), function(name) {
    current_path <- c(parent_names, name)
    element <- nested_list[[name]]
    if (typeof(element) == "S4") {
      element <- as.data.frame(element)
    }
    if (is.data.frame(element)) {
      # ignore data frames with no rows
      if (nrow(element) == 0) {
        return(NULL)
      }
      # if forceColumns is not empty, and the data frame dose not have the columns, add them
      if (length(forceColumns) > 0) {
        for (col in forceColumns) {
          if (!col %in% names(element)) {
            element[[col]] <- NA
          }
        }
      }
      # drop columns specified in columnsToDrop
      element <- element[, !names(element) %in% columnsToDrop]
      # If it's a data frame, add the path as new columns
      for (i in seq_along(current_path)) {
        element[[paste0("Level", i)]] <- current_path[i]
      }
      element$row <- rownames(element)
      rownames(element) <- NULL
      return(element)
    } else if (is.list(element)) {
      # If it's a list, recurse
      return(merge_nested_df(element, current_path, columnsToDrop))
    } else {
      stop("Encountered a non-data.frame and non-list element in the structure.")
    }
  }))
  if (!is.null(merged)) {
    merged <- dplyr::relocate(merged, row)
  }
  return(merged)
}

#' @title Create a DataTable with specified options
#'
#' @description This function creates a DataTable with the specified options using the DT and htmltools packages.
#'
#' @param data A data frame or data table to be displayed.
#' @param pageLength An integer specifying the number of rows to display per page. Default is 10.
#' @param dom A string specifying the DataTables DOM positioning. Default is "Bfrtip".
#' @param ... Additional arguments to be passed to the DT::datatable function.
#'
#' @return A tag list containing the DataTable.
#' @export
#' @importFrom DT datatable
#' @importFrom htmltools tagList
getTable <- function(data, pageLength = 10, dom = "Bfrtip", ...) {
  htmltools::tagList(DT::datatable(data, options = list(pageLength = pageLength, scrollX = TRUE, dom = dom), style = "bootstrap4", ...))
}

#' @title Get a DESeq2 Table with Specified Columns and Rounded Values
#'
#' @description This function filters and selects specific columns from a DESeq2 result data frame, rounds certain
#' columns to a specified number of decimal places, and then returns a DataTable of the result.
#'
#' @param data A data frame or data table containing DESeq2 result data.
#' @param shrunken A logical value indicating whether to use shrunken log2FoldChange and lfcSE values. Default is TRUE.
#'
#' @return A DataTable containing the filtered and selected DESeq2 result data with rounded values.
#' @export
#' @importFrom dplyr filter select mutate
#' @importFrom rlang .data
getDESeq2Table <- function(data, shrunken = TRUE) {
  showColumns <- c("SYMBOL", "log2FoldChange", "padj", "lfcSE", "sig_exp")
  if (shrunken) {
    # replace log2FoldChange with shrunken_log2FoldChange
    showColumns <- c("SYMBOL", "shrunken_log2FoldChange", "padj", "shrunken_lfcSE", "shrunken_sig_exp")
  }
  res <- data |>
    dplyr::select(dplyr::all_of(showColumns)) |>
    dplyr::rename_with(~ gsub("^shrunken_", "", .x)) |>
    dplyr::filter(.data$sig_exp != "NO") |>
    dplyr::mutate(
      log2FoldChange = round(.data$log2FoldChange, 2),
      padj = formatC(.data$padj, format = "e", digits = 3),
      lfcSE = round(.data$lfcSE, 3)
    )
  return(getTable(res))
}

#' @title Get a Limma Table with Specified Columns and Rounded Values
#'
#' @description This function filters and selects specific columns from a Limma result data frame,
#' rounds certain columns to a specified number of decimal places, and then returns a DataTable of the result.
#'
#' @param data A data frame or data table containing Limma result data.
#'
#' @return A DataTable containing the filtered and selected Limma result data with rounded values.
#' @export
#' @importFrom dplyr filter select mutate
#' @importFrom rlang .data
getLimmaTable <- function(data) {
  res <- data |>
    dplyr::filter(.data$sig_exp != "NO") |>
    dplyr::select(c("SYMBOL", "logFC", "adj.P.Val", "t", "sig_exp")) |>
    dplyr::mutate(
      logFC = round(.data$logFC, 2),
      adj.P.Val = round(.data$adj.P.Val, 5),
      t = round(.data$t, 2)
    )
  return(getTable(res))
}

#' @title Get an EdgeR Table with Specified Columns and Rounded Values
#'
#' @description This function filters and selects specific columns from an EdgeR result data frame,
#' rounds certain columns to a specified number of decimal places, and then returns a DataTable of the result.
#'
#' @param data A data frame or data table containing EdgeR result data.
#'
#' @return A DataTable containing the filtered and selected EdgeR result data with rounded values.
#' @export
#' @importFrom dplyr filter select mutate
#' @importFrom rlang .data
getEdgeRTable <- function(data) {
  res <- data |>
    dplyr::filter(.data$sig_exp != "NO") |>
    dplyr::select(c("SYMBOL", "logFC", "FDR", "F", "sig_exp")) |>
    dplyr::mutate(
      logFC = round(.data$logFC, 2),
      FDR = round(.data$FDR, 5),
      F = round(.data$F, 2)
    )
  return(getTable(res))
}

#' Save Homer List
#'
#' This function saves a list of genes that are upregulated or downregulated.
#' The function creates a directory, filters the data for upregulated and downregulated genes,
#' and writes these lists to separate text files.
#'
#' @param data A data frame containing gene information.
#' @param geneCol The name of the column in the data frame that contains the gene symbols. Default is "SYMBOL".
#' @param sigCol The name of the column in the data frame that contains the significance of gene expression. Default is "shrunken_sig_exp".
#' @param saveDir The directory where the output files will be saved.
#' @param name The base name for the output files. The function will append "_UP.txt" for upregulated genes and "_DOWN.txt" for downregulated genes.
#'
#' @return No return value. This function is used for its side effects of creating files.
#' @export
SaveHomerList <- function(data, geneCol = "SYMBOL", sigCol = "shrunken_sig_exp", saveDir, name) {
  dir.create(saveDir, showWarnings = FALSE)
  subData <- dplyr::select(data, !!sym(geneCol), !!sym(sigCol)) |>
    dplyr::filter(!!sym(sigCol) %in% c("UP", "DOWN"))
  upGenes <- dplyr::filter(subData, !!sym(sigCol) == "UP") |>
    dplyr::pull(!!sym(geneCol))
  downGenes <- dplyr::filter(subData, !!sym(sigCol) == "DOWN") |>
    dplyr::pull(!!sym(geneCol))
  write.table(upGenes, file.path(saveDir, paste0(name, "_UP.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(downGenes, file.path(saveDir, paste0(name, "_DOWN.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' Create Sample Sheet
#'
#' This function takes a directory of raw data, a path for the sample sheet, a name for the sample sheet, and a strandness.
#' It gets the list of fastq/fq files in the directory, groups the files by the sample name (directory name),
#' creates a data frame with the files, and writes the sample sheet to a file.
#'
#' @param rawDataDir A character string for the directory of raw data.
#' @param sampleSheetPath A character string for the path for the sample sheet.
#' @param sampleSheetName A character string for the name for the sample sheet.
#' @param strandedness A character string for the strandedness Default is "auto".
#' @param nameFrom A function to extract the sample name from the file name. Default is NULL.
#'
#' @return NULL. The function writes the sample sheet to a file and does not return anything.
#' @export
#'
#' @importFrom dplyr mutate rename
#' @importFrom tibble rownames_to_column
#' @importFrom utils write.table
CreateSampleSheet <- function(rawDataDir, sampleSheetPath, sampleSheetName, strandedness = "auto", nameFrom = NULL) {
  # Get the list of files in the directory
  files <- list.files(path = rawDataDir, full.names = TRUE, recursive = TRUE)
  # Get the list of files that are fastq/fq files
  fastqFiles <- files[grep(".fastq.gz|.fq.gz", files)]
  if (is.null(nameFrom)) {
    # Group the files by the sample name (directory name)
    fastqFiles <- split(fastqFiles, basename(dirname(fastqFiles)))
  } else {
    # Group the files by the name before the second underscore
    # fastqFiles <- split(fastqFiles, sub("^(.*?)_(.*?)_(.*)$", "\\1_\\2", basename(fastqFiles)))
    fastqFiles <- split(fastqFiles, nameFrom(fastqFiles))
  }
  # Create a data frame with the fastqFiles
  sampleSheet <- fastqFiles |>
    as.data.frame() |>
    t() |>
    as.data.frame() |>
    dplyr::rename(fastq_1 = "V1", fastq_2 = "V2") |>
    tibble::rownames_to_column("sample") |>
    dplyr::mutate(strandedness = strandedness)
  # Write the sample sheet to a file
  write.table(sampleSheet, file.path(sampleSheetPath, paste0(sampleSheetName, ".csv")), row.names = FALSE, quote = FALSE, sep = ",")
  return(sampleSheet)
}

#' Get Labels
#' This function filters the label data for a specific split type, sorts the data based on fold change and adjusted p-value,
#' and returns the top 5 genes based on fold change and p-value.
#' @param labelData A data frame containing the label data.
#' @param splitType A character string specifying the split type to filter the data.
#' @param sortOrder A character string specifying the sort order for the data. Default is "desc".
#' @param fc A character string specifying the column name in the data for the fold change.
#' @param padj A character string specifying the column name in the data for the adjusted p-value.
#' @param geneLabelCol A character string specifying the column name in the data for the gene labels.
#' 
#' @return A data frame containing the top 5 genes based on fold change and p-value.
#' @export
#' @importFrom dplyr filter arrange
#' @importFrom rlang sym
getLabels <- function(labelData, splitType = NULL, sortOrder, fc, padj, geneLabelCol) {
  if (!is.null(splitType)) {
    labelData <- dplyr::filter(labelData, split == splitType)
  }
  if (sortOrder == "desc") {
    labels <- labelData |>
      dplyr::arrange(dplyr::desc(!!sym(fc)), dplyr::desc(!!sym(padj)))
  } else {
    labels <- labelData |>
      dplyr::arrange(!!sym(fc), !!sym(padj))
  }
  labels <- labels |>
    dplyr::filter(!!sym(fc) * ifelse(sortOrder == "desc", 1, -1) > 0) |>
    dplyr::filter(!is.na(!!sym(geneLabelCol))) |>
    head(5)
  return(labels)
}