#' Get GEO Metadata
#'
#' This function retrieves the metadata for a given GEO ID.
#' It uses the GEOquery and Biobase packages to fetch and process the data.
#' If no count data or metadata is found for the given GEO ID, it will return a message or stop the function respectively.
#' The metadata is then cleaned up by renaming all columns, removing ":ch1" and replacing spaces with underscores.
#'
#' @param geoID A character string specifying the GEO ID for which to retrieve metadata.
#' @return A data frame containing the metadata for the given GEO ID.
#' @examples
#' \dontrun{
#' GetGEOMetadata("GSE12345")
#' }
#' @export
#' @importFrom GEOquery getGEO
#' @importFrom Biobase pData phenoData exprs assayData
#' @importFrom dplyr rename_all
GetGEOMetadata <- function(geoID) {
  # check geoID is a string starting with "GSE"
  if (!is.character(geoID) || !grepl("^GSE", geoID)) {
    stop("geoID should be a character string starting with 'GSE'.")
  }
  gse <- GEOquery::getGEO(GEO = geoID, GSEMatrix = TRUE, AnnotGPL = TRUE)
  metadata <- Biobase::pData(phenoData(gse[[1]]))
  if (nrow(Biobase::exprs(gse[[1]])) == 0) {
    message("No count data found for the given GEO ID.")
  }
  if (nrow(metadata) == 0) {
    stop("No metadata found for the given GEO ID.")
  }
  metadata <- dplyr::rename_all(metadata, ~ sub(":ch1", "", .)) |>
    dplyr::rename_all(~ gsub(" ", "_", .))
  return(metadata)
}

#' Create Summarized Experiment
#'
#' This function creates a SummarizedExperiment object from given counts and metadata.
#' It uses the SummarizedExperiment and DataFrame functions from the SummarizedExperiment package.
#'
#' @param counts A matrix or data frame of counts where rows are features (e.g., genes) and columns are samples.
#' @param metadata A data frame containing the metadata for the samples. The row names should match the column names of the counts.
#' @return A SummarizedExperiment object containing the counts and metadata.
#' @examples
#' \dontrun{
#' CreateSummarizedExperiment(counts, metadata)
#' }
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
CreateSummarizedExperiment <- function(counts, metadata) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("counts should be a matrix or a data frame.")
  }
  if (is.data.frame(counts)) {
    message("Converting counts to a matrix.")
    counts <- as.matrix(counts)
  }
  if (!is.data.frame(metadata)) {
    stop("metadata should be a data frame.")
  }
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(metadata)
  )
  return(se)
}

#' Get Salmon File List
#'
#' This function retrieves a list of Salmon quant.sf files from a given directory.
#' It uses the list.files function to search for files with the "quant.sf" pattern.
#' The names of the files in the list are set to the base name of the parent directory of each file.
#' If no Salmon quant.sf files are found in the given directory, the function will stop with a message.
#'
#' @param dirPath A character string specifying the directory path to search for Salmon quant.sf files.
#' @return A named list of Salmon quant.sf files. The names of the list elements are the base names of the parent directories of the files.
#' @examples
#' \dontrun{
#' GetSalmonFileList("/path/to/directory")
#' }
#' @export
GetSalmonFileList <- function(dirPath) {
  if (!dir.exists(dirPath)) {
    stop("The given directory path does not exist.")
  }
  fileList <- list.files(dirPath, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
  names(fileList) <- basename(dirname(fileList))
  if (length(fileList) == 0) {
    stop("No Salmon quant.sf files found in the given directory.")
  }
  return(fileList)
}

#' Get DESeqDataSet from tximport
#'
#' This function creates a DESeqDataSet object from a list of Salmon quant.sf files using the tximport and DESeq2 packages.
#' The tximport function is used to import transcript-level quantifications, and the DESeqDataSetFromTximport function is used to create a DESeqDataSet object.
#'
#' @param salmonFiles A named list of Salmon quant.sf files. The names of the list elements should be the sample names.
#' @param type A character string specifying the type of quantification. This should be "salmon" for Salmon quant.sf files.
#' @param tx2gene A two-column data frame linking transcript IDs to gene IDs. The first column should be named "TXNAME" and the second column should be named "GENEID".
#' @param metadata A data frame containing the metadata for the samples. The row names should match the names of the salmonFiles list.
#' @param design A formula specifying the experimental design for the differential expression analysis.
#' @param txOut A logical value indicating whether to output the transcript-level counts. Default is FALSE.
#' @return A DESeqDataSet object containing the imported transcript-level quantifications and the sample metadata.
#' @examples
#' \dontrun{
#' GetddsTxi(salmonFiles, "salmon", tx2gene, metadata, ~condition)
#' }
#' @export
#' @importFrom tximport tximport
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom plyr is.formula
GetddsTxi <- function(salmonFiles, type = "salmon", tx2gene, metadata, design, txOut = FALSE) {
  if (!is.data.frame(metadata)) {
    stop("metadata should be a data frame.")
  }
  if (!is.character(salmonFiles) || !all(names(salmonFiles) %in% rownames(metadata))) {
    stop("salmonFiles should be a named character vector with names matching the sample names in the metadata.")
  }
  if (!plyr::is.formula(design)) {
    stop("design should be a formula.")
  }
  txi <- tximport::tximport(
    files = salmonFiles,
    type = type,
    tx2gene = tx2gene,
    txOut = txOut
  )
  ddsTxi <- DESeq2::DESeqDataSetFromTximport(
    txi,
    colData = metadata,
    design = design
  )
  return(ddsTxi)
}

#' Filter DESeqDataSet
#'
#' This function filters a DESeqDataSet object to keep only the rows (features) that have a minimum number of counts in a minimum number of samples.
#' It uses the counts function from the DESeq2 package to get the count matrix from the DESeqDataSet object.
#'
#' @param ddsTxi A DESeqDataSet object.
#' @param minCounts An integer specifying the minimum number of counts a feature must have in a sample to be kept. Default is 10.
#' @param groupSize An integer specifying the minimum number of samples a feature must have the minimum number of counts in to be kept. Default is 1.
#' @return A DESeqDataSet object with the lowly expressed features filtered out.
#' @examples
#' \dontrun{
#' Filterdds(ddsTxi, minCounts = 10, groupSize = 1)
#' }
#' @export
#' @importFrom DESeq2 counts
Filterdds <- function(ddsTxi, minCounts = 10, groupSize = 1) {
  if (!inherits(ddsTxi, "DESeqDataSet")) {
    stop("ddsTxi should be a DESeqDataSet object.")
  }
  if (!is.numeric(groupSize) || groupSize < 1) {
    stop("groupSize should be a positive integer.")
  }
  if (groupSize > ncol(DESeq2::counts(ddsTxi))) {
    stop("groupSize should be less than the number of samples.")
  }
  if (!is.numeric(minCounts) || minCounts < 1) {
    stop("minCounts should be a positive integer.")
  }
  keep <- rowSums(DESeq2::counts(ddsTxi) >= minCounts) >= groupSize
  message(paste0(sum(!keep), " features filtered out."))
  ddsTxi <- ddsTxi[keep, ]
  message(paste0(nrow(ddsTxi), " features remaining."))
  return(ddsTxi)
}
