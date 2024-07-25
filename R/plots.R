#' Sample Distance Heatmap
#'
#' This function creates a heatmap of the sample distances for a given DESeqDataSet object.
#' It uses the vst function from the DESeq2 package to transform the count data, and the pheatmap function from the ComplexHeatmap package to create the heatmap.
#' The sample distances are calculated using the dist function.
#'
#' @param dds A DESeqDataSet object.
#' @param labelColumns A character vector specifying the column names in the colData of the DESeqDataSet to use for the sample labels in the heatmap. Default is NULL.
#' @return A heatmap of the sample distances.
#' @examples
#' \dontrun{
#' sampleDistanceHeatmap(dds, labelColumns = c("condition"))
#' }
#' @export
#' @importFrom DESeq2 vst
#' @importFrom ComplexHeatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats dist
#' @importFrom SummarizedExperiment assay colData
sampleDistanceHeatmap <- function(dds, labelColumns) {
  # Check if the 'dds' object is of the class 'DESeqDataSet'
  if (!inherits(dds, "DESeqDataSet")) {
    stop("dds should be a DESeqDataSet object.")
  }

  # Apply variance stabilizing transformation (VST) to the 'dds' object
  vsd <- DESeq2::vst(dds, blind = FALSE)
  # Compute the Euclidean distance between samples
  sampleDists <- dist(t(SummarizedExperiment::assay(vsd)))
  # Convert the distance object to a matrix
  sampleDistMatrix <- as.matrix(sampleDists)
  # Get the column names of the 'vsd' object
  labels <- colnames(vsd)
  
  # If 'labelColumns' is provided
  if (!missing(labelColumns)) {
    # Loop over each element in 'labelColumns'
    for (i in labelColumns) {
      # If 'i' is a column name in the 'vsd' object
      if (i %in% colnames(SummarizedExperiment::colData(vsd))) {
        # Concatenate the current labels with the values in column 'i'
        labels <- paste(labels, colData(vsd)[[i]], sep = "-")
      } else {
        # If 'i' is not a column name in the 'vsd' object, stop the execution and display an error message
        stop(paste("Invalid label provided.", i, "is not a column in the colData of the DESeqDataSet."))
      }
    }
  }
  
  # Remove the column names of the distance matrix
  colnames(sampleDistMatrix) <- NULL
  # Set the row names of the distance matrix to 'labels'
  rownames(sampleDistMatrix) <- labels
  # Generate a color palette
  colors <- grDevices::colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  # Create a heatmap of the distance matrix
  ComplexHeatmap::pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
  )
}

#' PCA Plot
#'
#' This function creates a PCA plot for a given DESeqDataSet object.
#' It uses the vst and plotPCA functions from the DESeq2 package to transform the count data and create the PCA plot.
#'
#' @param dds A DESeqDataSet object.
#' @param intgroup A character vector specifying the column names in the colData of the DESeqDataSet to use for the
#' sample grouping in the PCA plot. Default is c("Condition", "Sex").
#' @param title A character string specifying the title of the plot. Default is an empty string.
#' @param batch A boolean flag indicating whether to remove batch effects. Default is FALSE.
#' @param formula A formula specifying the design matrix for the batch effect. Default is NULL.
#' @param batchCol A character string specifying the column name in the colData of the DESeqDataSet to use for the batch effect. Default is NULL.
#' @param remove A character vector specifying the sample names to remove from the plot. Default is NULL.
#' @return A PCA plot.
#' @examples
#' \dontrun{
#' pcaPlot(dds, intgroup = c("Condition", "Sex"))
#' }
#' @export
#' @importFrom DESeq2 vst plotPCA
#' @importFrom ggplot2 ggplot aes geom_point geom_text theme_bw
#' @importFrom rlang sym .data
#' @importFrom stats model.matrix
#' @importFrom limma removeBatchEffect
#' @importFrom SummarizedExperiment assay colData
pcaPlot <- function(dds, intgroup = c("Condition", "Sex"), title = "", batch = FALSE, formula = NULL, batchCol = NULL, remove = NULL) {
  # Check if the 'dds' object is of the class 'DESeqDataSet'
  if (!inherits(dds, "DESeqDataSet")) {
    stop("dds should be a DESeqDataSet object.")
  }
  
  # Apply variance stabilizing transformation (VST) to the 'dds' object
  vsd <- DESeq2::vst(dds, blind = FALSE)

  if (!is.null(remove)) {
    # remove samples passed in the vector remove
    vsd <- vsd[, !(colnames(vsd) %in% remove)]
  }
  
  if (batch) {
    mat <- assay(vsd)
    mm <- model.matrix(formula, colData(vsd))
    mat <- limma::removeBatchEffect(mat, batch = vsd[[batchCol]], design = mm)
    SummarizedExperiment::assay(vsd) <- mat
  }
  
  # Perform Principal Component Analysis (PCA) on the 'vsd' object
  # 'intgroup' specifies the variables to plot
  # 'returnData = TRUE' returns the data used for plotting
  # 'ntop = 500' uses the top 500 genes by variance
  pca <- DESeq2::plotPCA(vsd, intgroup = intgroup, returnData = TRUE, ntop = 500)

  # Create a scatter plot of the first two principal components
  # The points are colored by the first variable in 'intgroup'
  # The labels of the points are the row names of the 'pca' object
  plot <- ggplot2::ggplot(pca, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = !!sym(intgroup[1]))) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = rownames(pca)), nudge_x = 0.6, nudge_y = 0.6) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title)
  
  # If 'intgroup' has more than one variable
  if (length(intgroup) > 1) {
    # Add a different shape to the points for the second variable in 'intgroup'
    plot <- plot + ggplot2::aes(shape = !!sym(intgroup[2]))
  }
  # Return the plot
  plot
}

#' Volcano Plot
#'
#' This function creates a volcano plot for a given data set.
#' It uses the ggplot2 and ggrepel packages to create the plot and add labels to the points.
#'
#' @param data A data frame containing the data to plot. It should have columns for the fold change, adjusted p-value, and significance expression.
#' @param fc A character string specifying the column name in the data for the fold change.
#' @param padj A character string specifying the column name in the data for the adjusted p-value.
#' @param sig_exp A character string specifying the column name in the data for the significance expression.
#' @param yintercept A numeric value specifying the y-intercept for the horizontal line in the plot. Default is -log10(0.05).
#' @param xintercept A numeric vector specifying the x-intercepts for the vertical lines in the plot. Default is c(-0.6, 0.6).
#' @param title A character string specifying the title of the plot. Default is "".
#' @param xlab A character string specifying the x-axis label. Default is "log2(Fold Change)".
#' @param ylab A character string specifying the y-axis label. Default is "-log10(P-adj)".
#' @param plotly A logical value indicating whether to convert the ggplot object to a plotly object. Default is FALSE.
#' @param colors A list specifying the colors to use for the points in the plot. Default is list(UP = "#00806B", DOWN = "#E1AD25", NO = "grey").
#' @param maxPvalue A numeric value specifying the maximum p-value to include in the plot. If specified, points with a p-value greater
#' than thisvalue will be removed from the plot.
#' @param label A character vector specifying the gene symbols to label in the plot. If not specified, the top 5 and bottom 5 genes
#' based on fold change and p-value will be labeled.
#' @param highlightLabels A logical value indicating whether to highlight the labeled points in red. Default is FALSE.
#' @param geneLabelCol A character string specifying the column name in the data for the gene names. Default is "SYMBOL".
#' @param geneIDCol A character string specifying the column name in the data for the gene IDs. Default is "rownames".
#' @param maxLogFC A numeric value specifying the maximum absolute log2 fold change to include in the plot. Default is 10.
#' @param splitAt A numeric value specifying the threshold for splitting the plot into two parts based on the adjusted p-value.
#' @return A ggplot or plotly object of the volcano plot.
#' @examples
#' \dontrun{
#' volcanoPlot(data, "log2FoldChange", "padj", "sig_exp")
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual geom_vline geom_hline coord_cartesian labs theme_bw theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom plotly ggplotly
#' @importFrom dplyr filter arrange pull
#' @importFrom rlang sym
volcanoPlot <- function(
    data, fc, padj, sig_exp, yintercept = -log10(0.05), xintercept = c(-0.6, 0.6), title = "", xlab = "log2(Fold Change)",
    ylab = "-log10(P-adj)", plotly = FALSE, colors = list(UP = "#00806B", DOWN = "#E1AD25", NO = "grey"),
    maxPvalue = 20, maxLogFC = 10, label, highlightLabels = FALSE, geneLabelCol = "SYMBOL", geneIDCol = "rownames", splitAt = NULL) {
  # If the geneIDCol is "rownames", convert the rownames to a column
  if (geneIDCol == "rownames") {
    data <- tibble::rownames_to_column(data, var = geneIDCol)
  }
  
  # Create a temporary symbol column in the data
  data <- dplyr::mutate(data, tmpSymbol = ifelse(is.na(!!sym(geneLabelCol)), !!sym(geneIDCol), !!sym(geneLabelCol)))
  # Add genes with absolute fold change greater than or equal to maxLogFC to the removedGenes
  removedGenes <- dplyr::filter(data, abs(!!sym(fc)) >= maxLogFC) |> dplyr::pull("tmpSymbol")
  # Add genes with -log10(padj) greater than maxPvalue to the removedGenes
  removedGenes <- c(removedGenes, dplyr::filter(data, -log10(!!sym(padj)) > maxPvalue) |> dplyr::pull("tmpSymbol"))
  
  # Filter the data to keep only genes with -log10(padj) less than or equal to maxPvalue
  data <- dplyr::filter(data, -log10(!!sym(padj)) <= maxPvalue)
  # Further filter the data to keep only genes with absolute fold change less than maxLogFC and non-NA padj
  data <- dplyr::filter(data, abs(!!sym(fc)) < maxLogFC & !is.na(padj))
  
  # Filter the data to keep only genes with significant expression (UP or DOWN)
  labelData <- dplyr::filter(data, !!sym(sig_exp) %in% c("UP", "DOWN"))
  # Add a new column 'split' to labelData based on the splitAt value
  if (!is.null(splitAt)) {
    labelData$split <- ifelse(-log10(labelData[[padj]]) <= splitAt, "Below", "Above")
  }
  # If labels are provided, filter the labelData to keep only those genes
  if (!missing(label)) {
    labelData <- dplyr::filter(labelData, !!sym(geneLabelCol) %in% label)
  } else {
    if (!is.null(splitAt)) {
      # Use the helper function to get labels
      upLabels_below <- getLabels(labelData, "Below", "desc", fc, padj, geneLabelCol)
      downLabels_below <- getLabels(labelData, "Below", "asc", fc, padj, geneLabelCol)
      upLabels_above <- getLabels(labelData, "Above", "desc", fc, padj, geneLabelCol)
      downLabels_above <- getLabels(labelData, "Above", "asc", fc, padj, geneLabelCol)
      # Combine all labels
      labelData <- rbind(upLabels_below, downLabels_below, upLabels_above, downLabels_above)
    } else {
      upLabels <- getLabels(labelData, NULL, "desc", fc, padj, geneLabelCol)
      downLabels <- getLabels(labelData, NULL, "asc", fc, padj, geneLabelCol)
      labelData <- rbind(upLabels, downLabels)
    }
  }

  if (!is.null(splitAt)) {
    data$split <- ifelse(-log10(data[[padj]]) <= splitAt, "Below", "Above")
    data_below <- dplyr::filter(data, split == "Below")
    data_above <- dplyr::filter(data, split == "Above")

    p_below <- ggplot2::ggplot(data_below, ggplot2::aes(!!sym(fc), -log10(!!sym(padj)), , colour = !!sym(sig_exp))) +
      ggplot2::geom_point(alpha = 0.4, size = 1.75) +
      ggplot2::scale_colour_manual(values = colors) +
      ggplot2::geom_vline(xintercept = xintercept, linetype = "longdash") +
      ggplot2::geom_hline(yintercept = yintercept, linetype = "longdash") +
      ggplot2::coord_cartesian(xlim = c(min(data[[fc]]), max(data[[fc]]))) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::theme_bw() +
      ggrepel::geom_text_repel(data = dplyr::filter(labelData, split == "Below"), ggplot2::aes(label = !!sym(geneLabelCol)), color = "black", box.padding = 0.5, max.overlaps = Inf) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8))

    p_above <- ggplot2::ggplot(data_above, ggplot2::aes(!!sym(fc), -log10(!!sym(padj)), , colour = !!sym(sig_exp))) +
      ggplot2::geom_point(alpha = 0.4, size = 1.75) +
      ggplot2::scale_colour_manual(values = colors) +
      ggplot2::geom_vline(xintercept = xintercept, linetype = "longdash") +
      ggplot2::coord_cartesian(xlim = c(min(data[[fc]]), max(data[[fc]]))) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = title) +
      ggrepel::geom_text_repel(data = dplyr::filter(labelData, split == "Above"), ggplot2::aes(label = !!sym(geneLabelCol)), color = "black", box.padding = 0.5, max.overlaps = Inf) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8)) +
      ggplot2::labs(subtitle = if (length(removedGenes) > 0 & length(removedGenes) <= 10) {
        paste("Genes with |log2(FC)| >=", maxLogFC, " or -log10(", padj, ") > ", maxPvalue, " :", paste(removedGenes, collapse = ", "))
      } else if (length(removedGenes) > 10) {
        paste("Genes with |log2(FC)| >=", maxLogFC, " or -log10(", padj, ") > ", maxPvalue, " : ", length(removedGenes), " genes removed")
      }) +
      ggplot2::scale_y_log10() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank()
      ) +
      ggplot2::theme(legend.position = "none")

    p <- p_above / p_below
    p <- p + patchwork::plot_layout(heights = c(1, 3))
  } else {
    # Create a volcano plot with ggplot2
    p <- ggplot2::ggplot(data, ggplot2::aes(!!sym(fc), -log10(!!sym(padj)), , colour = !!sym(sig_exp))) +
      ggplot2::geom_point(alpha = 0.4, size = 1.75) +
      {
        if (highlightLabels) ggplot2::geom_point(data = labelData, colour = "red", size = 2.5)
      } + # nolint: brace_linter.
      ggplot2::scale_colour_manual(values = colors) +
      ggplot2::geom_vline(xintercept = xintercept, linetype = "longdash") +
      ggplot2::geom_hline(yintercept = yintercept, linetype = "longdash") +
      ggplot2::coord_cartesian(xlim = c(min(data[[fc]]), max(data[[fc]]))) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_bw() +
      ggrepel::geom_text_repel(data = labelData, ggplot2::aes(label = !!sym(geneLabelCol)), color = "black", box.padding = 0.5, max.overlaps = Inf) +
      ggplot2::labs(subtitle = if (length(removedGenes) > 0 & length(removedGenes) <= 10) {
        paste("Genes with |log2(FC)| >=", maxLogFC, " or -log10(", padj, ") > ", maxPvalue, " :", paste(removedGenes, collapse = ", "))
      } else if (length(removedGenes) > 10) {
        paste("Genes with |log2(FC)| >=", maxLogFC, " or -log10(", padj, ") > ", maxPvalue, " : ", length(removedGenes), " genes removed")
      }) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8))
  }

  # If plotly is TRUE, convert the ggplot to a plotly plot
  if (plotly) {
    p <- plotly::ggplotly(p)
  }
  # Return the plot
  return(p)
}

#' Plot GSE Comparison Lollipop
#'
#' This function takes a data frame and creates a lollipop plot for GSE comparison.
#' It allows for optional filtering of the data, and customization of the plot title and colors.
#'
#' @param data A data frame containing GSE data.
#' @param nameCol A character string specifying the column name to be used for the lollipop labels.
#' @param title A character string specifying the title of the plot. Default is an empty string.
#' @param filter A list specifying the filter conditions for each group. Default is NULL.
#' @param splitCol A character string specifying the column name to be used for splitting the data into groups.
#' @param colors A character vector specifying the colors to be used for the groups. Default is c("#94609F", "#96CECE").
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read.csv("gse_data.csv")
#' plot <- plotGSEComparisonLollipop(
#'   data, "Description", "GSE Comparison",
#'   list(Group1 = c("Pathway1", "Pathway2"), Group2 = c("Pathway3", "Pathway4")), "Group"
#' )
#' print(plot)
#' }
#' @importFrom dplyr filter arrange mutate bind_rows
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_size labs theme_bw ggtitle scale_fill_gradient scale_y_discrete
#' @importFrom stringr str_trunc
#' @importFrom rlang sym .data
plotGSEComparisonLollipop <- function(data, nameCol, title = "", filter, splitCol, colors = c("#94609F", "#96CECE")) {
  # If a filter is provided
  if (!missing(filter)) {
    # Apply the filter to the data
    # For each element in the filter, select the rows where 'nameCol' is in the filter element and 'splitCol' equals the filter name
    tmpDf <- lapply(names(filter), function(i) {
      dplyr::filter(data, !!sym(nameCol) %in% filter[[i]] & !!sym(splitCol) == i)
    })
    # Combine the filtered data frames into one data frame
    data <- dplyr::bind_rows(tmpDf)
  }

  # Arrange the data by 'NES' and convert 'nameCol' to a factor with levels in the order they appear in the data
  data <- dplyr::arrange(data, .data$NES) |>
    dplyr::mutate(!!sym(nameCol) := factor(!!sym(nameCol), levels = unique(!!sym(nameCol))))

  # Create a lollipop chart with ggplot2
  data |>
    ggplot2::ggplot(aes(x = .data$NES, y = .data$Description, fill = -log10(.data$p.adjust))) +
    # Add segments from x = 0 to x = 'NES' for each 'Description'
    ggplot2::geom_segment(ggplot2::aes(x = 0, y = .data$Description, xend = .data$NES, yend = .data$Description), linewidth = 1) +
    # Add points at (x = 'NES', y = 'Description') with size proportional to 'setSize'
    ggplot2::geom_point(ggplot2::aes(size = .data$setSize), pch = 21) +
    # Scale the size of the points from 4 to 12
    ggplot2::scale_size(range = c(4, 12)) +
    # Set the labels for the x and y axes and the color legend
    ggplot2::labs(x = "NES", y = "Pathway", color = "-log10(P adjusted)") +
    # Use a white background
    ggplot2::theme_bw() +
    # Set the title of the plot
    ggplot2::ggtitle(title) +
    # Use a gradient color scale from blue to red for the fill color
    ggplot2::scale_fill_gradient(low = "#1c6fff", high = "#fd3b3b") +
    # Truncate the y labels to 50 characters
    ggplot2::scale_y_discrete(label = function(x) stringr::str_trunc(x, 50))
}

#' Expression Plot
#'
#' This function takes a DESeqDataSet, a vector of gene symbols, and a boolean flag for bar plot.
#' It checks if the rownames of the DESeqDataSet are Ensembl gene IDs and converts the gene symbols to Ensembl IDs if necessary.
#' It then calculates the normalized counts for each gene, averages the counts for genes sharing the same symbol,
#' and creates a bar plot or a jitter plot of the normalized counts.
#'
#' @param dds A DESeqDataSet containing the count data.
#' @param normalized A boolean flag indicating whether to use normalized counts. Default is TRUE.
#' @param genes A character vector of gene symbols.
#' @param bar A boolean flag indicating whether to create a bar plot. Default is FALSE.
#' @param color A character string specifying the column name in the colData of the DESeqDataSet to use for the color of
#' the points in the plot. Default is "Sex"
#' @param shape A character string specifying the column name in the colData of the DESeqDataSet to use for the shape of the points in the plot. Default is NULL.
#' @param colors A character vector specifying the colors to be used in the plot. Default is c("#94609F", "#96CECE").
#' @param title A character string specifying the title of the plot. Default is "Gene expression".
#' @param nrow An integer specifying the number of rows in the facet grid. Default is NULL.
#' @param ncol An integer specifying the number of columns in the facet grid. Default is NULL.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
#' genes <- c("BRCA1", "BRCA2", "TP53")
#' plot <- expressionPlot(dds, genes, bar = TRUE)
#' print(plot)
#' }
#' @importFrom BiocGenerics counts
#' @importFrom clusterProfiler bitr
#' @importFrom dplyr filter left_join group_by summarise_if ungroup
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual labs facet_wrap theme_bw geom_jitter scale_color_manual
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble rownames_to_column
#' @importFrom rlang .data
expressionPlot <- function(dds, normalized = TRUE, genes, type = "scatter", color = "Sex", shape = NULL, colors = c("#94609F", "#96CECE"), title = "Gene expression", nrow = NULL, ncol = NULL) {
  # Check if all row names of the dataframe 'dds' start with "ENSG" (human Ensembl gene IDs) or "ENSMUSG" (mouse Ensembl gene IDs)
  if (all(grepl("^ENSG", rownames(dds)))) {
    # If true, convert gene symbols to Ensembl IDs using the human database
    genes <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db") |>
      # Filter out rows with NA Ensembl IDs
      dplyr::filter(!is.na(.data$ENSEMBL))
  } else if (all(grepl("^ENSMUSG", rownames(dds)))) {
    # If all row names start with "ENSMUSG", convert gene symbols to Ensembl IDs using the mouse database
    genes <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Mm.eg.db") |>
      # Filter out rows with NA Ensembl IDs
      dplyr::filter(!is.na(.data$ENSEMBL))
  } else {
    # If row names do not match either of the above patterns, print a message
    message("The rownames of the dataframe are not Ensembl gene IDs.")
    countData <- BiocGenerics::counts(dds, normalized = normalized) |>
      as.data.frame() |>
      tibble::rownames_to_column("SYMBOL") |>
      dplyr::filter(.data$SYMBOL %in% genes)
  }

  if (all(grepl("^ENSG", rownames(dds))) | all(grepl("^ENSMUSG", rownames(dds)))) {
    # Get normalized counts from 'dds', convert to dataframe, and add Ensembl IDs as a column
    countData <- BiocGenerics::counts(dds, normalized = normalized) |>
      as.data.frame() |>
      tibble::rownames_to_column("ENSEMBL") |>
      # Filter rows where Ensembl ID is in 'genes'
      dplyr::filter(.data$ENSEMBL %in% genes$ENSEMBL) |>
      # Join 'genes' dataframe by Ensembl ID
      dplyr::left_join(genes, by = "ENSEMBL") |>
      # Group by gene symbol
      dplyr::group_by(.data$SYMBOL) |>
      # For each group, calculate the mean of numeric columns, ignoring NA values
      dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) |>
      # Ungroup the dataframe
      dplyr::ungroup()
  }

  # Convert dataframe from wide to long format, with 'SYMBOL' as the identifier variable
  # Join with the column data of 'dds' by 'Sample'
  countData <- tidyr::pivot_longer(countData, -.data$SYMBOL, names_to = "Sample", values_to = "Normalized_Counts") |>
    dplyr::left_join(SummarizedExperiment::colData(dds), by = "Sample", copy = TRUE)

  # Create a ggplot object with 'SYMBOL' on the x-axis and 'Normalized_Counts' on the y-axis
  # If 'bar' is true, create a bar plot
  if (type == "bar") {
    p <- ggplot2::ggplot(countData, ggplot2::aes(x = .data$SYMBOL, y = .data$Normalized_Counts, label = .data$Sample, fill = !!sym(color))) +
      # Add bars with fill color determined by 'color'
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual("legend", values = c("Female" = "#94609F", "Male" = "#96CECE")) +
      ggplot2::labs(x = paste0(if (normalized) "Normalized " else "", "Counts"), y = "Sample") +
      # Facet the plot by 'SYMBOL', with each facet in a separate row
      ggplot2::facet_wrap(~ .data$SYMBOL, scales = "free", nrow = nrow, ncol = ncol) +
      ggplot2::theme_bw()
  } else if (type == "boxplot") {
    if (!is.null(shape)) {
      p <- ggplot2::ggplot(countData, ggplot2::aes(x = .data$SYMBOL, y = .data$Normalized_Counts, fill = !!sym(color), shape = !!sym(shape)))
    } else {
      p <- ggplot2::ggplot(countData, ggplot2::aes(x = .data$SYMBOL, y = .data$Normalized_Counts, fill = !!sym(color)))
    }
    p <- p +
      # Add bars with fill color determined by 'color'
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge()) +
      ggplot2::scale_fill_manual("legend", values = colors) +
      ggplot2::labs(y = paste0(if (normalized) "Normalized " else "", "Counts"), x = "Gene") +
      # Facet the plot by 'SYMBOL', with each facet in a separate row
      ggplot2::facet_wrap(~ .data$SYMBOL, scales = "free", nrow = nrow, ncol = ncol) +
      ggplot2::theme_bw()
  } else {
    p <- ggplot2::ggplot(countData, ggplot2::aes(x = .data$SYMBOL, y = .data$Normalized_Counts, label = .data$Sample))
    # If 'shape' is not NULL, add it as an aesthetic
    if (!is.null(shape)) {
      p <- p +
        # Add points
        ggplot2::geom_jitter(ggplot2::aes(color = !!sym(color), shape = !!sym(shape), size = 2), width = 0.3, alpha = 0.8)
    } else {
      p <- p +
        # Add points without 'shape'
        ggplot2::geom_jitter(ggplot2::aes(color = !!sym(color), size = 2), width = 0.3, alpha = 0.8)
    }
    p <- p +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(x = "Gene", y = paste0(if (normalized) "Normalized " else "", "Counts")) +
      # Facet the plot by 'SYMBOL', with each facet in a separate column
      ggplot2::facet_wrap(~ .data$SYMBOL, scales = "free", nrow = nrow, ncol = ncol) +
      ggplot2::theme_bw() +
      # Remove x axis labels and ticks
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  }
  # Add a title to the plot and use a white background
  p <- p +
    ggplot2::labs(title = title)
  # Return the plot
  p
}

#' P-value Distribution
#'
#' This function takes a numeric vector of p-values and creates a histogram of the p-values.
#' The histogram has 50 breaks, a sky blue color, a white border, and labels for the main title and the x-axis.
#'
#' @param pvalues A numeric vector of p-values.
#'
#' @return A histogram of the p-values.
#' @export
#'
#' @examples
#' \dontrun{
#' pvalues <- runif(1000, min = 0, max = 1)
#' pvalueDistribution(pvalues)
#' }
#' @importFrom graphics hist
pvalueDistribution <- function(pvalues) {
  hist(pvalues, breaks = 50, col = "skyblue", border = "white", main = "Histogram of p-values", xlab = "p-value")
}

#' Get Venn Diagrams of DE Methods
#'
#' This function takes a list of DE results from different methods (DESeq2, edgeR, and limma) and a boolean flag for shrunken.
#' It creates three Venn diagrams: one for significant genes, one for upregulated genes, and one for downregulated genes.
#' The Venn diagrams are then plotted side by side.
#'
#' @param DE A list of DE results from DESeq2, edgeR, and limma.
#' @param shrunken A boolean flag indicating whether to use shrunken significant expression. Default is TRUE.
#'
#' @return A ggplot object of the Venn diagrams.
#' @export
#'
#' @examples
#' \dontrun{
#' DE <- list(DESeq2 = DESeq2_results, edgeR = edgeR_results, limma = limma_results)
#' venns <- getVennMethodsDE(DE, shrunken = TRUE)
#' print(venns)
#' }
#' @importFrom dplyr filter
#' @importFrom ggVennDiagram ggVennDiagram
#' @importFrom ggplot2 ggtitle
#' @importFrom rlang .data
getVennMethodsDE <- function(DE, shrunken = TRUE) {
  # Check if the 'shrunken' variable is TRUE
  if (shrunken) {
    # If TRUE, set 'DESeq2_sig_exp' to "shrunken_sig_exp"
    DESeq2_sig_exp <- "shrunken_sig_exp"
  } else {
    # If FALSE, set 'DESeq2_sig_exp' to "sig_exp"
    DESeq2_sig_exp <- "sig_exp"
  }

  # Create a Venn diagram comparing significant genes (either upregulated or downregulated) identified by DESeq2, edgeR, and limma
  v1 <- ggVennDiagram::ggVennDiagram(
    list(
      DESeq2 = dplyr::filter(DE$DESeq2, !!sym(DESeq2_sig_exp) %in% c("UP", "DOWN")) |> rownames(),
      edgeR = dplyr::filter(DE$edgeR, .data$sig_exp %in% c("UP", "DOWN")) |> rownames(),
      limma = dplyr::filter(DE$limma, .data$sig_exp %in% c("UP", "DOWN")) |> rownames()
    )
  ) + ggplot2::ggtitle("DESeq2, edgeR, and limma significant genes")  # Add a title to the plot
  # Create a Venn diagram comparing upregulated genes identified by DESeq2, edgeR, and limma
  v2 <- ggVennDiagram::ggVennDiagram(
    list(
      DESeq2 = dplyr::filter(DE$DESeq2, !!sym(DESeq2_sig_exp) == "UP") |> rownames(),
      edgeR = dplyr::filter(DE$edgeR, .data$sig_exp == "UP") |> rownames(),
      limma = dplyr::filter(DE$limma, .data$sig_exp == "UP") |> rownames()
    )
  ) + ggplot2::ggtitle("DESeq2, edgeR, and limma upregulated genes")  # Add a title to the plot
  # Create a Venn diagram comparing downregulated genes identified by DESeq2, edgeR, and limma
  v3 <- ggVennDiagram::ggVennDiagram(
    list(
      DESeq2 = dplyr::filter(DE$DESeq2, !!sym(DESeq2_sig_exp) == "DOWN") |> rownames(),
      edgeR = dplyr::filter(DE$edgeR, .data$sig_exp == "DOWN") |> rownames(),
      limma = dplyr::filter(DE$limma, .data$sig_exp == "DOWN") |> rownames()
    )
  ) + ggplot2::ggtitle("DESeq2, edgeR, and limma downregulated genes")  # Add a title to the plot

  # Combine the three Venn diagrams into a single plot with three columns
  venns <- cowplot::plot_grid(v1, v2, v3, ncol = 3)
  # Return the combined plot
  venns
}

#' Heatmap Plot
#'
#' This function takes a DESeqDataSet, a DESeqResults object, a tx2gene data frame, a title, a boolean flag for z-score, and a score column.
#' It calculates the variance stabilizing transformation of the counts, selects the top and bottom genes based on the score column,
#' joins the transformed counts with the gene symbols, calculates the mean of the counts for genes sharing the same symbol,
#' and creates a heatmap of the counts.
#'
#' @param dds A DESeqDataSet containing the count data.
#' @param DEResults A DESeqResults object containing the DE results.
#' @param tx2gene A data frame containing the gene symbols and Ensembl IDs.
#' @param title A character string for the title of the heatmap. Default is "".
#' @param zscore A boolean flag indicating whether to z-score the counts. Default is TRUE.
#' @param scoreCol A character string for the score column. Default is "shrunken_log2FoldChange".
#' @param idCol A character string for the column name to use as the rownames in the heatmap. Default is "rownames".
#' @param symbolCol A character string for the column name to use as the gene symbols in the heatmap. Default is "SYMBOL".
#' @param annotationCols A character vector specifying the column names in the colData of the DESeqDataSet to use for the column
#' annotations in the heatmap. Default is c("Condition", "Sex").
#' @param remove A character vector specifying the sample names to remove from the DESeqDataSet. Default is NULL.
#' @param show An integer specifying the number of genes to show in each side of the heatmap. Default is 25.
#' @param ... Additional arguments to be passed to the pheatmap function.
#'
#' @return A heatmap of the counts.
#' @export
#'
#' @examples
#' \dontrun{
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
#' DEResults <- DESeq2::DESeq(dds)
#' tx2gene <- data.frame(GENESYMBOL = c("BRCA1", "BRCA2", "TP53"), 
#'  GENEID = c("ENSG00000012048", "ENSG00000139618", "ENSG00000141510"))
#' heatmap <- heatmapPlot(dds, DEResults, tx2gene, title = "Heatmap", 
#'  zscore = TRUE, scoreCol = "shrunken_log2FoldChange")
#' print(heatmap)
#' }
#' @importFrom ComplexHeatmap pheatmap
#' @importFrom DESeq2 vst
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr arrange desc pull tibble select rename filter group_by summarize_all ungroup
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom SummarizedExperiment assay colData
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
heatmapPlot <- function(
  dds, DEResults, tx2gene, title = "", zscore = TRUE, scoreCol = "shrunken_log2FoldChange", idCol = "rownames", symbolCol = "SYMBOL", 
  annotationCols = c("Condition", "Sex"), remove = NULL, show = 25, ...
) {
  # Apply variance stabilizing transformation (VST) to the DESeqDataSet
  vsd <- DESeq2::vst(dds, blind = FALSE)
  metadata <- SummarizedExperiment::colData(dds)
  if (!is.null(remove)) {
    # remove samples specified in the remove vector
    vsd <- vsd[, !colnames(vsd) %in% remove]
    metadata <- metadata[!rownames(metadata) %in% remove, ]
  }
  # Convert the VST data to a data frame
  vsd <- SummarizedExperiment::assay(vsd) |> as.data.frame()

  # Arrange the DEResults data frame by the score column in descending order
  DEResults <- dplyr::arrange(DEResults, dplyr::desc(!!sym(scoreCol)))
  # Filter out rows where the symbol column is NA
  genes <- dplyr::filter(DEResults, !is.na(!!sym(symbolCol)))
  # If there are more than 50 genes, select the top 25 and bottom 25 based on the score column
  top <- dplyr::arrange(genes, dplyr::desc(!!sym(scoreCol))) |> head(show)
  bottom <- dplyr::arrange(genes, !!sym(scoreCol)) |> head(show)
  genes <- base::rbind(top, bottom)
  

  # Handle the ID column assignment based on user input
  if (idCol == "rownames") {
    genes <- tibble::rownames_to_column(genes, var = "ENSEMBL")
  } else {
    genes$ENSEMBL <- genes[[idCol]]
  }

  # Select and join data frames for plotting
  genes <- dplyr::select(genes, "ENSEMBL", !!sym(symbolCol))
  data <- dplyr::inner_join(
    vsd |> tibble::rownames_to_column("ENSEMBL"),
    genes,
    by = "ENSEMBL"
  ) |>
    dplyr::select(-c("ENSEMBL")) |>
    dplyr::arrange(match(!!sym(symbolCol), DEResults[[symbolCol]])) |>
    tibble::column_to_rownames(var = symbolCol)

  # Create a data frame of column annotations based on user-specified columns
  # Use syms() instead of sym() for multiple column names
  annotation_col <- metadata |>
      as.data.frame() |>
      dplyr::select(!!!rlang::syms(annotationCols))

  # Define colors for the column annotations dynamically
  annotation_colors <- lapply(annotationCols, function(col) {
      unique_values <- unique(metadata[[col]])
      color_pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(8, length(unique_values)), "Set1"))(length(unique_values))
      names(color_pal) <- unique_values
      color_pal
  })
  names(annotation_colors) <- annotationCols

  # Scaling data if z-score is requested
  if (zscore) {
    data <- base::t(base::scale(base::t(data)))
  }

  # Generate heatmap
  data |>
    as.matrix() |>
    ComplexHeatmap::pheatmap(
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      angle_col = "0",
      main = title,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      scale = if (zscore) "row" else "none",
      name = if (zscore) "Z-score" else "VST",
      color = if (zscore) {
        rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
      } else {
        grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
      },
      ...
    )
}

#' Over Representation Bar Plot
#'
#' This function takes a data frame, column names for name, score, and count, a title, a filter, and colors.
#' It selects the specified columns, pivots the data frame longer, arranges the rows by descending p-value,
#' converts the name column to a factor with levels in the order of appearance, and converts the p-value to -log10 scale.
#' If a filter is provided, it filters the rows by the name column.
#' It then creates a bar plot of the p-values and a bar plot of the counts, and combines the two plots side by side.
#'
#' @param data A data frame containing the data.
#' @param nameCol A character string for the name column.
#' @param scoreCol A character string for the score column.
#' @param countCol A character string for the count column.
#' @param title A character string for the title of the plots. Default is "".
#' @param filter A character vector for filtering the rows by the name column. Default is NULL.
#' @param colors A character vector for the colors of the bars. Default is "blue".
#'
#' @return A ggplot object of the combined plots.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(Pathway = c("Pathway1", "Pathway2", "Pathway3"), 
#'  Pvalue = c(0.01, 0.05, 0.001), Count = c(10, 20, 30))
#' plot <- plotOverRepBar(data, "Pathway", "Pvalue", "Count", 
#'  title = "Over Representation Bar Plot", filter = c("Pathway1", "Pathway2"), colors = c("blue"))
#' print(plot)
#' }
#' @importFrom dplyr select arrange desc mutate filter
#' @importFrom ggplot2 ggplot aes geom_col labs theme_bw scale_fill_manual theme scale_x_continuous element_blank element_text
#' @importFrom tidyr pivot_longer
#' @importFrom patchwork plot_annotation
#' @importFrom scales pretty_breaks
#' @importFrom rlang .data
plotOverRepBar <- function(data, nameCol, scoreCol, countCol, title = "", filter, colors = c("blue")) {
  # Select the specified columns from the data frame and pivot the data so that each row is a unique combination of name and count
  # The remaining columns are transformed into two new columns: 'comparison' and 'pvalue'
  # The rows are then arranged in descending order of pvalue, and the pvalue column is transformed to -log10 scale
  data <- dplyr::select(data, !!sym(nameCol), !!sym(scoreCol), !!sym(countCol)) |>
      tidyr::pivot_longer(cols = -c(!!sym(nameCol), !!sym(countCol)), names_to = "comparison", values_to = "pvalue") |>
      dplyr::arrange(dplyr::desc(.data$pvalue)) |>
      dplyr::mutate(
        !!sym(nameCol) := factor(!!sym(nameCol), levels = unique(!!sym(nameCol))),
        pvalue = -log10(.data$pvalue)
      )

  # If a filter is provided, filter the data to include only the specified names
  if (!missing(filter)) {
    data <- dplyr::filter(data, !!sym(nameCol) %in% filter)
  }

  # Create a bar plot of pvalue versus name, with bars colored according to the 'colors' vector
  pPlot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$pvalue, y = !!sym(nameCol))) +
      ggplot2::geom_col(position = "dodge", fill = colors) +
      ggplot2::labs(x = "-log10(Pvalue)", y = "Pathway") +
      ggplot2::theme_bw()
  # Create a bar plot of count versus name, with gray bars
  countPlot <- ggplot2::ggplot(data, ggplot2::aes(x = !!sym(countCol), y = !!sym(nameCol))) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "# of genes") +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_manual(values = "gray") +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ) +
      # x axis display only whole numbers
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6))

  # Combine the pvalue and count plots side by side
  barPlots <- pPlot + countPlot
  # Add a title to the combined plot and center it
  combinedPlots <- barPlots + patchwork::plot_annotation(title = title, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  # Return the combined plot
  combinedPlots
}

#' Plots the number of reads per sample
#'
#' This function takes a list of log files from Salmon (a tool for transcript quantification from RNA-seq data),
#' extracts the total number of mapped reads from each log file, and plots a bar chart of the number of reads per sample.
#'
#' @param salmonFilesLog A list of file paths to the Salmon log files.
#'
#' @return A ggplot2 object representing the bar chart of the number of reads per sample.
#' @export
#'
#' @examples
#' \dontrun{
#' readsPerSamplePlot(c("path/to/log1", "path/to/log2"))
#' }
#' @importFrom ggplot2 aes geom_bar labs coord_flip ggtitle
#' @importFrom rlang .data
readsPerSamplePlot <- function(salmonFilesLog) {
  # get line from log file with number of reads
  numReads <- function(file) {
    readLines(file)[grep("Total # of mapped reads", readLines(file))]
  }
  numReadsTxt <- lapply(salmonFilesLog, numReads)
  numReads <- as.numeric(gsub("\\D", "", unlist(numReadsTxt)))
  names(numReads) <- names(numReadsTxt)
  df <- data.frame(Sample = names(numReads), NumReads = numReads)
  # plot number of reads per sample using ggplot2
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Sample, y = .data$NumReads)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(x = "Sample", y = "Number of reads") +
    ggplot2::coord_flip() +
    ggplot2::ggtitle("Number of reads per sample")
}

#' Plots the number of genes per sample
#'
#' This function takes a DESeqDataSet object that contains the transcript quantifications for each sample,
#' calculates the number of genes detected in each sample, and plots a bar chart of the number of genes per sample.
#'
#' @param ddsTxi An unfiltered DESeqDataSet object that contains the gene quantifications for each sample.
#'
#' @return A ggplot2 object representing the bar chart of the number of genes per sample.
#' @export
#'
#' @examples
#' \dontrun{
#' genesPerSamplePlot(ddsTxi)
#' }
#' @importFrom ggplot2 aes geom_bar labs coord_flip ggtitle
#' @importFrom rlang .data
genesPerSamplePlot <- function(ddsTxi) {
  counts <- t(counts(ddsTxi)) |> as.data.frame()
  numGenes <- rowSums(counts > 0)
  numGenes <- data.frame(Sample = rownames(counts), NumGenes = numGenes)
  ggplot2::ggplot(numGenes, ggplot2::aes(x = .data$Sample, y = .data$NumGenes)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(x = "Sample", y = "Number of genes") +
    ggplot2::coord_flip() +
    ggplot2::ggtitle("Number of genes detected per sample")
}

#' Plots the number of transcripts per sample
#'
#' This function takes a DESeqDataSet object that contains the transcript quantifications for each sample,
#' calculates the number of transcripts detected in each sample, and plots a bar chart of the number of transcripts per sample.
#'
#' @param ddsTxi An unfiltered DESeqDataSet object that contains the transcript quantifications for each sample.
#'
#' @return A ggplot2 object representing the bar chart of the number of transcripts per sample.
#' @export
#'
#' @examples
#' \dontrun{
#' transcriptsPerSamplePlot(ddsTxi)
#' }
#' @importFrom ggplot2 aes geom_bar labs
#' @importFrom rlang .data
transcriptsPerSamplePlot <- function(ddsTxi) {
  counts <- t(counts(ddsTxi)) |> as.data.frame()
  numTranscripts <- rowSums(counts > 0)
  numTranscripts <- data.frame(Sample = rownames(counts), NumTranscripts = numTranscripts)
  ggplot2::ggplot(numTranscripts, ggplot2::aes(x = .data$Sample, y = .data$NumTranscripts)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(x = "Sample", y = "Number of transcripts") +
    ggplot2::coord_flip() +
    ggplot2::ggtitle("Number of transcripts detected per sample")
}

#' Plots the distribution of Transcripts Per Million (TPM) for each sample
#'
#' This function takes a tximport object that contains the transcript quantifications for each sample,
#' calculates the TPM for each gene in each sample, and plots a density plot of the TPM values.
#'
#' @param txi A tximport object that contains the transcript quantifications for each sample.
#'
#' @return A ggplot2 object representing the density plot of the TPM values.
#' @export
#'
#' @examples
#' \dontrun{
#' tpmDistributionPlot(txi)
#' }
#' @importFrom dplyr mutate
#' @importFrom ggplot2 aes geom_density ggtitle
#' @importFrom tidyr gather
tpmDistributionPlot <- function(txi) {
  tpm <- t(txi$abundance) |> as.data.frame()
  tpm <- dplyr::mutate(tpm, Sample = rownames(tpm)) |>
    tidyr::gather("Gene", "TPM", -"Sample")
  ggplot2::ggplot(tpm, ggplot2::aes(x = log10(.data$TPM), fill = .data$Sample)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::ggtitle("TPM distribution")
}

homerHeatmapPlot <- function(dds, genes, tx2gene, title = "", zscore = FALSE) {
  vsd <- vst(dds, blind = FALSE)
  vsd <- assay(vsd) |> as.data.frame()
  data <- dplyr::left_join(
    vsd |> tibble::rownames_to_column("ENSEMBL"),
    tx2gene |> dplyr::select("GENESYMBOL", "GENEID") |> dplyr::rename(ENSEMBL = GENEID),
    by = "ENSEMBL",
    multiple = "first"
  ) |>
    select(-ENSEMBL) |>
    dplyr::filter(GENESYMBOL %in% genes) |>
    dplyr::group_by(GENESYMBOL) |>
    dplyr::summarize_all(mean) |>
    dplyr::ungroup() |>
    dplyr::arrange(match(GENESYMBOL, genes)) |>
    tibble::column_to_rownames("GENESYMBOL")
  # add annotation to columns in the heatmap
  annotation_col <- colData(dds) |>
    as.data.frame() |>
    dplyr::select(Condition, Sex)
  annotation_colors <- list(
    Condition = c("Room_Air" = "#E1AD25", "Oxygen" = "#00806B"),
    Sex = c("Male" = "#96CECE", "Female" = "#94609F")
  )
  if (zscore) {
    data <- t(scale(t(data)))
  }
  data |>
    as.matrix() |>
    ComplexHeatmap::pheatmap(
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      angle_col = "0",
      main = title,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      breaks = seq(0, 10, length.out = 256),
      scale = "none",
      name = if (zscore) "Z-score" else "VST",
      color = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
    )
}

plotPathwayHeatmap <- function(gse) {
  gse <- dplyr::filter(gse, Input == "all") |>
    dplyr::select(c("Sex", "ID", "NES")) |>
    tidyr::pivot_wider(names_from = "ID", values_from = "NES") |>
    dplyr::mutate(across(everything(), ~ replace(., is.na(.), 0))) |>
    tibble::column_to_rownames("Sex")

  ComplexHeatmap::Heatmap(
    gse,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    name = "NES",
    column_title = paste0(ncol(gse), " Pathways"),
    column_title_side = "bottom",
    left_annotation = ComplexHeatmap::rowAnnotation(
      Sex = c("Female", "Male"),
      col = list(Sex = c("Female" = "#94609F", "Male" = "#96CECE"))
    ),
    col = circlize::colorRamp2(c(-4, 0, 4), c("#E1AD25", "white", "#00806B")),
    border = TRUE
  )
}

plotOverRepVenn <- function(data) {
  data <- dplyr::filter(data, Input == "all")
  male <- dplyr::filter(data, Sex == "Male") |> pull(ID)
  female <- dplyr::filter(data, Sex == "Female") |> pull(ID)
  ggVennDiagram::ggVennDiagram(list(Male = male, Female = female))
}

