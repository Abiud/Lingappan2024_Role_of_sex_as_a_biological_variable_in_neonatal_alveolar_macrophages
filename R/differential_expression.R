#' Get Differential Expression Analysis with DESeq2
#'
#' This function performs differential expression analysis on a given DESeqDataSet using DESeq2. It also applies
#' log2 fold change shrinkage using the apeglm method.
#'
#' @param dds A DESeqDataSet object containing the count data and the design formula.
#' @param comparisonName A character string specifying the name of the comparison to be used in the results function of DESeq2.
#' @param alpha A numeric value specifying the significance level for the adjusted p-values.
#' @param foldChangeThreshold A numeric value specifying the log2 fold change threshold for determining significant genes.
#' @param verbose A logical value indicating whether the MA plots should be displayed. Default is TRUE.
#' @return A data frame containing the results of the differential expression analysis. It includes columns for the statistics,
#' p-values, adjusted p-values, and significance expression for both the original and shrunken log2 fold changes.
#' @examples
#' \dontrun{
#' GetDESeq2DE(dds, "condition_treatment_vs_control", 0.05, 1)
#' }
#' @export
#' @importFrom DESeq2 results lfcShrink
#' @importFrom dplyr select rename_with left_join
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom BiocGenerics plotMA
#' @importFrom graphics par
GetDESeq2DE <- function(dds, comparisonName, alpha, foldChangeThreshold, verbose = TRUE) {
  if (missing(alpha)) {
    stop("The alpha argument is missing.")
  }
  if (missing(foldChangeThreshold)) {
    stop("The foldChangeThreshold argument is missing.")
  }
  if (!inherits(dds, "DESeqDataSet")) {
    stop("dds should be a DESeqDataSet object.")
  }
  res <- DESeq2::results(dds, name = comparisonName, alpha = alpha)
  # Get shrunken log2 fold changes
  resShrunk <- DESeq2::lfcShrink(dds, coef = comparisonName, type = "apeglm")

  if (verbose) {
    graphics::par(mfrow = c(1, 2))
    BiocGenerics::plotMA(res, alpha = alpha, main = paste0("MA plot for ", comparisonName))
    BiocGenerics::plotMA(resShrunk, alpha = alpha, main = paste0("MA plot for shrunken ", comparisonName))
  }

  res <- as.data.frame(res) |>
    tibble::rownames_to_column("geneId")

  resShrunk <- as.data.frame(resShrunk) |>
    dplyr::select(.data$log2FoldChange, .data$lfcSE) |>
    dplyr::rename_with(~ paste0("shrunken_", .x)) |>
    tibble::rownames_to_column("geneId")

  res <- dplyr::left_join(res, resShrunk, by = "geneId") |>
    GetGeneFormat()

  res$padj[is.na(res$padj)] <- .Machine$double.xmin

  res <- addSignificance(res, alpha = alpha, fc = foldChangeThreshold, padjCol = "padj", fcCol = "log2FoldChange") |>
    addSignificance(alpha = alpha, fc = foldChangeThreshold, padjCol = "padj", fcCol = "shrunken_log2FoldChange", sigCol = "shrunken_sig_exp")
  return(res)
}

#' Add Significance Expression to Results
#'
#' This function adds a new column to the results of differential expression analysis indicating the significance expression of each gene.
#'
#' @param res A data frame containing the results of differential expression analysis. It should have columns for the adjusted p-values and log2 fold changes.
#' @param alpha A numeric value specifying the significance level for the adjusted p-values.
#' @param fc A numeric value specifying the log2 fold change threshold for determining significant genes.
#' @param padjCol A character string specifying the column name in res for the adjusted p-values.
#' @param fcCol A character string specifying the column name in res for the log2 fold changes.
#' @param sigCol A character string specifying the name of the new column to be added for the significance expression. Default is "sig_exp".
#' @return A data frame containing the results of the differential expression analysis with the new column for the significance expression.
#' @examples
#' \dontrun{
#' addSignificance(res, 0.05, 1, "padj", "log2FoldChange", "sig_exp")
#' }
#' @export
#' @importFrom dplyr mutate case_when
#' @importFrom rlang !! sym .data :=
addSignificance <- function(res, alpha, fc, padjCol, fcCol, sigCol = "sig_exp") {
  if (!inherits(res, "data.frame")) {
    stop("res should be a data.frame object.")
  }
  if (missing(alpha)) {
    stop("The alpha argument is missing.")
  }
  if (missing(fc)) {
    stop("The fc argument is missing.")
  }
  # Check if required columns are present in the dataframe
  if (!(padjCol %in% names(res))) {
    stop(paste("The padjCol argument:", padjCol, "is not a column name in the dataframe."))
  }
  if (!(fcCol %in% names(res))) {
    stop(paste("The fcCol argument:", fcCol, "is not a column name in the dataframe."))
  }
  res %>%
    dplyr::mutate(
      !!sigCol := dplyr::case_when(
        .data[[padjCol]] < alpha & .data[[fcCol]] > fc ~ "UP",
        .data[[padjCol]] < alpha & .data[[fcCol]] < -fc ~ "DOWN",
        TRUE ~ "NO"
      )
    )
}

#' Get Gene Format
#'
#' This function converts the gene IDs in the results of differential expression analysis to gene symbols using the clusterProfiler
#' package. It supports conversion for human and mouse genes.
#'
#' @param res A data frame containing the results of differential expression analysis. It should have a column for the gene IDs.
#' @return A data frame containing the results of the differential expression analysis with the gene IDs converted to gene symbols.
#' @examples
#' \dontrun{
#' GetGeneFormat(res)
#' }
#' @export
#' @importFrom clusterProfiler bitr
#' @importFrom dplyr rename left_join
#' @importFrom tibble column_to_rownames
GetGeneFormat <- function(res) {
  if (!inherits(res, "data.frame")) {
    stop("res should be a data.frame object.")
  }
  # check if rownames are human gene symbols or mouse gene symbols
  symbols <- NULL
  if (all(grepl("^ENSG", res$geneId))) {
    symbols <- clusterProfiler::bitr(res$geneId, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
  } else if (all(grepl("^ENSMUSG", res$geneId))) {
    symbols <- clusterProfiler::bitr(res$geneId, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
  } else {
    message("The rownames of the dataframe are not Ensembl gene IDs.")
  }
  res <- dplyr::rename(res, ENSEMBL = "geneId")
  if (!is.null(symbols)) {
    res <- dplyr::left_join(res, symbols, by = "ENSEMBL", multiple = "first")
  }
  res <- tibble::column_to_rownames(res, "ENSEMBL")
  return(res)
}

#' Get DGE Object
#'
#' This function creates a DGEList object from count data and group information, filters out lowly expressed genes, and calculates normalization factors.
#'
#' @param dds A DESeqDataSet object containing the count data and metadata information.
#' @param group A character string specifying the column name in colData for the group information.
#' @return A DGEList object containing the count data, group information, and normalization factors.
#' @examples
#' \dontrun{
#' GetDGEObject(counts, group)
#' }
#' @export
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
GetDGEObject <- function(dds, group = "Experiment") {
  if (missing(dds)) {
    stop("argument \"dds\" is missing, with no default")
  }
  if (!inherits(dds, "DESeqDataSet")) {
    stop("dds should be a DESeqDataSet object.")
  }
  if (!(group %in% colnames(colData(dds)))) {
    stop("group should be a column in colData.")
  }
  if (!is.factor(colData(dds)[[group]])) {
    stop("group should be a factor. With first level as control and second level as treatment.")
  }
  # Create a DGEList object, which is a container for storing count data and associated information for differential expression analysis
  dge <- edgeR::DGEList(counts = counts(dds), samples = colData(dds))
  dge$samples$group <- factor(colData(dds)[[group]])
  # Filter out lowly expressed genes
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, ]
  # Calculate normalization factors to scale the raw library sizes
  dge <- edgeR::calcNormFactors(dge)
  return(dge)
}

#' Get Differential Expression Analysis with edgeR
#'
#' This function performs differential expression analysis on a given DGEList using edgeR. It also adds gene symbols
#' to the results and a new column indicating the significance expression of each gene.
#'
#' @param dge A DGEList object containing the count data, group information, and normalization factors.
#' @param alpha A numeric value specifying the significance level for the adjusted p-values.
#' @param foldChangeThreshold A numeric value specifying the log2 fold change threshold for determining significant genes.
#' @param designMatrix A formula specifying the experimental design for the differential expression analysis.
#' @param comparisonName A character string specifying the name of the comparison to be used in the results function of edgeR.
#' @param verbose A logical value indicating whether the contrast matrix should be displayed. Default is TRUE.
#' @return A data frame containing the results of the differential expression analysis. It includes columns for the
#' statistics, p-values, adjusted p-values.
#'
#' @export
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom stats model.matrix
#' @importFrom knitr kable
GetEdgeRDE <- function(dge, alpha, foldChangeThreshold, designMatrix, comparisonName, verbose = TRUE) {
  if (!inherits(dge, "DGEList")) {
    stop("dge should be a DGEList object.")
  }
  if (missing(alpha)) {
    stop("The alpha argument is missing.")
  }
  if (missing(foldChangeThreshold)) {
    stop("The foldChangeThreshold argument is missing.")
  }
  if (missing(designMatrix)) {
    stop("The designMatrix argument is missing.")
  }
  if (missing(comparisonName)) {
    stop("The comparisonName argument is missing.")
  }
  # Create a design matrix, which specifies the experimental design for the differential expression analysis
  design <- stats::model.matrix(designMatrix, data = dge$samples)
  # rename (Intercept) to Intercept
  colnames(design)[1] <- "Intercept"
  # Estimate the dispersion parameter for each gene
  dge <- edgeR::estimateDisp(dge, design)
  # Fit the statistical model
  fit <- edgeR::glmQLFit(dge, design)
  # Perform the statistical test for differential expression
  qlf <- edgeR::glmQLFTest(fit, contrast = limma::makeContrasts(contrasts = comparisonName, levels = design))
  if (verbose) {
    message("Contrast matrix:")
    print(knitr::kable(design))
    message("Contrast: ", comparisonName)
  }
  # Extract the top-ranked genes
  res <- edgeR::topTags(qlf, n = Inf)$table
  # Add gene symbols to the results
  res <- tibble::rownames_to_column(res, "geneId") |>
    GetGeneFormat() |>
    addSignificance(alpha = alpha, fc = foldChangeThreshold, padjCol = "FDR", fcCol = "logFC")
  return(res)
}

#' Get Differential Expression Analysis with limma
#'
#' This function performs differential expression analysis on a given DGEList using limma. It also adds gene symbols
#' to the results and a new column indicating the significance expression of each gene.
#'
#' @param dge A DGEList object containing the count data, group information, and normalization factors.
#' @param alpha A numeric value specifying the significance level for the adjusted p-values.
#' @param foldChangeThreshold A numeric value specifying the log2 fold change threshold for determining significant genes.
#' @param designMatrix A formula specifying the experimental design for the differential expression analysis.
#' @param comparisonName A character string specifying the name of the comparison to be used in the results function of limma.
#' @param verbose A logical value indicating whether the contrast matrix should be displayed. Default is TRUE.
#' @return A data frame containing the results of the differential expression analysis. It includes columns for the
#' statistics, p-values, adjusted p-values.
#' @examples
#' \dontrun{
#' GetLimmaDE(dge, 0.05, 1)
#' }
#' @export
#' @importFrom limma voom lmFit eBayes makeContrasts topTable
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom stats model.matrix
GetLimmaDE <- function(dge, alpha, foldChangeThreshold, designMatrix, comparisonName, verbose = TRUE) {
  if (!inherits(dge, "DGEList")) {
    stop("dge should be a DGEList object.")
  }
  if (missing(alpha)) {
    stop("The alpha argument is missing.")
  }
  if (missing(foldChangeThreshold)) {
    stop("The fc argument is missing.")
  }
  if (missing(designMatrix)) {
    stop("The designMatrix argument is missing.")
  }
  if (missing(comparisonName)) {
    stop("The comparisonName argument is missing.")
  }
  # Transform the count data to log2-counts per million (CPM), estimate the mean-variance relationship and compute the observational-level weights
  v <- limma::voom(dge, plot = verbose, design = stats::model.matrix(designMatrix, data = dge$samples))
  # rename (Intercept) to Intercept
  # colnames(v$design)[1] <- "Intercept"
  # Fit the linear model
  fit <- limma::lmFit(v)
  # Compute empirical Bayes statistics
  fit <- limma::eBayes(fit)
  contrast.matrix <- limma::makeContrasts(contrasts = comparisonName, levels = v$design)
  if (verbose) {
    message("Contrast matrix:")
    print(knitr::kable(v$design))
    message("Contrast: ", comparisonName)
    print(knitr::kable(contrast.matrix))
  }
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  res <- limma::topTable(fit2, adjust = "fdr", n = Inf)
  # Add gene symbols to the results
  res <- tibble::rownames_to_column(res, "geneId") |>
    GetGeneFormat() |>
    addSignificance(alpha = alpha, fc = foldChangeThreshold, padjCol = "adj.P.Val", fcCol = "logFC")
  return(res)
}

#' Get Differential Expression Results
#'
#' This function gets the differential expression results from DESeq2, edgeR, and limma packages.
#' It then merges these results into a single summary data frame.
#'
#' @param dds A DESeqDataSet object.
#' @param dge A DGEList object.
#' @param comparisonNameDESeq2 A character string specifying the name of the comparison to be used in the results function of DESeq2.
#' @param comparisonNameDGE A character string specifying the name of the comparison to be used in the results function of edgeR and limma.
#' @param designDGE A formula specifying the experimental design for the differential expression analysis in edgeR and limma.
#' @param alpha A numeric value specifying the significance level for the differential expression tests.
#' @param foldChangeThreshold A numeric value specifying the fold change threshold for the differential expression tests.
#' @param verbose A logical value indicating whether the contrast matrix should be displayed. Default is TRUE.
#'
#' @return A list containing the differential expression results from DESeq2, edgeR, and limma, and a summary data frame.
#' @export
#'
#' @examples
#' \dontrun{
#' dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
#' dge <- DGEList(counts = counts, group = group)
#' res <- GetDE(dds, dge, "condition", 0.05, 1)
#' }
GetDE <- function(dds, dge, comparisonNameDESeq2, comparisonNameDGE, designDGE, alpha, foldChangeThreshold, verbose = TRUE) {
  resDESeq2 <- GetDESeq2DE(dds, comparisonNameDESeq2, alpha, foldChangeThreshold, verbose = verbose)
  resEdgeR <- GetEdgeRDE(dge, alpha, foldChangeThreshold, designMatrix = designDGE, comparisonName = comparisonNameDGE, verbose = verbose)
  resLimma <- GetLimmaDE(dge, alpha, foldChangeThreshold, designMatrix = designDGE, comparisonName = comparisonNameDGE, verbose = verbose)
  forceColumns <- c(colnames(resDESeq2), colnames(resEdgeR), colnames(resLimma))
  summary <- merge_nested_df(list(DESeq2 = resDESeq2, edgeR = resEdgeR, limma = resLimma), forceColumns = forceColumns) |>
    dplyr::rename(Package = .data$Level1, ENSEMBL = .data$row) |>
    dplyr::relocate(.data$Package, .after = "ENSEMBL")
  return(list(DESeq2 = resDESeq2, edgeR = resEdgeR, limma = resLimma, summary = summary))
}
