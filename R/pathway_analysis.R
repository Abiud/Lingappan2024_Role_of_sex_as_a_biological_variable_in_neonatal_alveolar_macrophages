#' Get Pathways GSE
#'
#' This function performs gene set enrichment analysis (GSEA) on a given set of genes using the clusterProfiler and ReactomePA packages.
#' It supports analysis for human and mouse genes.
#'
#' @param res A data frame containing the gene data. It should have columns for the statistics and gene names.
#' @param statCol A character string specifying the column name in res for the statistics. Default is "stat".
#' @param geneNameCol A character string specifying the column name in res for the gene names. If "rownames",
#' the row names of res will be used as the gene names. Default is "rownames".
#' @param organism A character string specifying the organism. It should be "human" or "mouse".
#' @return A list of data frames containing the GSEA results for the GO, KEGG, Reactome, and Hallmark pathways.
#' @examples
#' \dontrun{
#' GetPathwaysGSE(res, "stat", "rownames", "human")
#' }
#' @export
#' @importFrom clusterProfiler gseGO gseKEGG GSEA
#' @importFrom ReactomePA gsePathway
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
GetPathwaysGSE <- function(res, statCol = "stat", geneNameCol = "rownames", organism) {
  orgDBs <- list(
    mouse = org.Mm.eg.db::org.Mm.eg.db,
    human = org.Hs.eg.db::org.Hs.eg.db
  )
  geneList <- res[, statCol]
  if (geneNameCol == "rownames") {
    names(geneList) <- rownames(res)
  } else {
    names(geneList) <- as.character(res[, geneNameCol])
  }
  geneList <- sort(geneList, decreasing = TRUE)
  if (all(grepl("^ENSG", names(geneList))) || all(grepl("^ENSMUSG", names(geneList)))) {
    keyType <- "ENSEMBL"
  } else {
    keyType <- "SYMBOL"
  }
  geneList <- TranslateToENTREZ(geneList, keyType = keyType, orgDB = orgDBs[[organism]])

  GOGSE <- data.frame()
  KeggGSE <- data.frame()
  ReactGSE <- data.frame()
  HallmarkGSE <- data.frame()

  tryCatch(
    {
      GOGSE <- clusterProfiler::gseGO(
        geneList = geneList,
        OrgDb = orgDBs[[organism]],
        ont = "BP",
        minGSSize = 5,
        verbose = FALSE,
        eps = 1e-300
      )
    },
    error = function(e) {
      message("Error in GOGSE: ", e$message)
    }
  )
  tryCatch(
    {
      KeggGSE <- clusterProfiler::gseKEGG(
        geneList = geneList,
        organism = switch(organism,
          mouse = "mmu",
          human = "hsa"
        ),
        keyType = "ncbi-geneid",
        minGSSize = 5,
        verbose = FALSE,
        eps = 1e-300
      )
    },
    error = function(e) {
      message("Error in KeggGSE: ", e$message)
    }
  )
  tryCatch(
    {
      ReactGSE <- ReactomePA::gsePathway(
        geneList = geneList,
        organism = organism,
        minGSSize = 5,
        verbose = FALSE,
        eps = 1e-300
      )
    },
    error = function(e) {
      message("Error in ReactGSE: ", e$message)
    }
  )

  m_t2g <- GetHallmarkDB(organism)
  tryCatch(
    {
      HallmarkGSE <- clusterProfiler::GSEA(
        geneList,
        TERM2GENE = m_t2g,
        minGSSize = 5,
        verbose = FALSE,
        eps = 1e-300
      )
    },
    error = function(e) {
      message("Error in HallmarkGSE: ", e$message)
    }
  )
  return(SetReadablePathways(list(GO = GOGSE, Kegg = KeggGSE, Reactome = ReactGSE, Hallmark = HallmarkGSE), orgDB = orgDBs[[organism]]))
}

#' Get Overrepresented Pathways
#'
#' This function performs overrepresentation analysis (ORA) on a given set of genes using the clusterProfiler, ReactomePA, and org.Hs.eg.db packages.
#' It supports analysis for human and mouse genes.
#'
#' @param res A data frame containing the gene data. It should have columns for the statistics and gene names.
#' @param statCol A character string specifying the column name in res for the statistics. Default is "stat".
#' @param geneNameCol A character string specifying the column name in res for the gene names. If "rownames", the row names
#' of res will be used as the gene names. Default is "rownames".
#' @param organism A character string specifying the organism. It should be "human" or "mouse".
#' @return A list of data frames containing the ORA results for the GO, KEGG, Reactome, and Hallmark pathways.
#' @examples
#' \dontrun{
#' GetPathwaysOverRep(res, "stat", "rownames")
#' }
#' @export
#' @importFrom clusterProfiler enrichGO enrichKEGG enricher
#' @importFrom ReactomePA enrichPathway
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
GetPathwaysOverRep <- function(res, statCol = "stat", geneNameCol = "rownames", organism) {
  orgDBs <- list(
    "mouse" = org.Mm.eg.db::org.Mm.eg.db,
    "human" = org.Hs.eg.db::org.Hs.eg.db
  )
  geneList <- res[, statCol]
  if (geneNameCol == "rownames") {
    names(geneList) <- rownames(res)
  } else {
    names(geneList) <- as.character(res[, geneNameCol])
  }
  geneList <- sort(geneList, decreasing = TRUE)
  if (all(grepl("^ENSG", names(geneList))) || all(grepl("^ENSMUSG", names(geneList)))) {
    keyType <- "ENSEMBL"
  } else {
    keyType <- "SYMBOL"
  }
  geneList <- TranslateToENTREZ(geneList, keyType = keyType, orgDB = orgDBs[[organism]])

  GOOverRep <- data.frame()
  KeggOverRep <- data.frame()
  ReactOverRep <- data.frame()
  HallmarkOverRep <- data.frame()

  tryCatch(
    {
      GOOverRep <- enrichGO(
        gene = names(geneList),
        OrgDb = orgDBs[[organism]],
        ont = "BP",
        minGSSize = 5
      )
    },
    error = function(e) {
      message("Error in GOOverRep: ", e$message)
    }
  )
  tryCatch(
    {
      KeggOverRep <- enrichKEGG(
        gene = names(geneList),
        organism = switch(organism,
          mouse = "mmu",
          human = "hsa"
        ),
        minGSSize = 5
      )
    },
    error = function(e) {
      message("Error in KeggOverRep: ", e$message)
    }
  )
  tryCatch(
    {
      ReactOverRep <- ReactomePA::enrichPathway(
        gene = names(geneList),
        organism = organism,
        readable = TRUE,
        minGSSize = 5
      )
    },
    error = function(e) {
      message("Error in ReactOverRep: ", e$message)
    }
  )

  m_t2g <- GetHallmarkDB(organism)
  tryCatch(
    {
      HallmarkOverRep <- enricher(
        names(geneList),
        TERM2GENE = m_t2g,
        minGSSize = 5
      )
    },
    error = function(e) {
      message("Error in HallmarkOverRep: ", e$message)
    }
  )
  return(SetReadablePathways(list(GO = GOOverRep, Kegg = KeggOverRep, Reactome = ReactOverRep, Hallmark = HallmarkOverRep), orgDB = orgDBs[[organism]]))
}

#' Get Hallmark Database
#'
#' This function retrieves the Hallmark gene set database for a given organism using the msigdbr package.
#' It supports analysis for human and mouse genes.
#'
#' @param organism A character string specifying the organism. It should be "human" or "mouse".
#' @return A data frame containing the Hallmark gene set database for the specified organism.
#' @examples
#' \dontrun{
#' GetHallmarkDB("human")
#' }
#' @export
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr select
#' @importFrom rlang .data
GetHallmarkDB <- function(organism) {
  msigdbr(species = switch(organism,
    mouse = "Mus musculus", # nolint: indentation_linter.
    human = "Homo sapiens"
  ), category = "H") %>%
    dplyr::select(.data$gs_name, .data$entrez_gene)
}

#' Translate To ENTREZ
#'
#' This function translates a list of gene identifiers to ENTREZ IDs using a specified key type and organism database.
#'
#' @param geneList A named numeric vector where the names are gene identifiers and the values are the corresponding statistics.
#' @param keyType A character string specifying the type of the gene identifiers in geneList. Default is "ENSEMBL".
#' @param orgDB An organism database object from the AnnotationDbi package.
#' @return A named numeric vector where the names are ENTREZ IDs and the values are the corresponding statistics.
#' @examples
#' \dontrun{
#' TranslateToENTREZ(geneList, "ENSEMBL", org.Hs.eg.db)
#' }
#' @export
#' @importFrom clusterProfiler bitr
#' @importFrom dplyr filter
#' @importFrom rlang .data
TranslateToENTREZ <- function(geneList, keyType = "ENSEMBL", orgDB) {
  geneNames <- names(geneList)
  tryCatch(
    {
      # geneNames <- AnnotationDbi::select(orgDB, geneNames, "ENTREZID", keyType)
      geneNames <- clusterProfiler::bitr(geneNames, fromType = keyType, toType = "ENTREZID", OrgDb = orgDB)
    },
    error = function(e) {
      message("Error in TranslateToENTREZ: ", e$message)
    }
  )
  if ("ENTREZID" %in% colnames(geneNames)) {
    geneNames <- geneNames |> dplyr::filter(!is.na(.data$ENTREZID))
    names(geneList) <- geneNames$ENTREZID[match(names(geneList), geneNames[[keyType]])]
    geneList <- geneList[!is.na(names(geneList))]
    geneList <- sort(geneList, decreasing = TRUE)
    return(geneList)
  } else {
    return(c())
  }
}

#' Set Readable Pathways
#'
#' This function sets the gene symbols in the results of pathway analysis to be readable using a specified key type and organism database.
#'
#' @param pathways A list of data frames containing the results of pathway analysis. The data frames should have a column for the gene symbols.
#' @param orgDB An organism database object from the AnnotationDbi package.
#' @param keyType A character string specifying the type of the gene symbols in the data frames. Default is "ENTREZID".
#' @return A list of data frames containing the results of pathway analysis with readable gene symbols.
#' @examples
#' \dontrun{
#' SetReadablePathways(pathways, org.Hs.eg.db, "ENTREZID")
#' }
#' @export
#' @importFrom clusterProfiler setReadable
SetReadablePathways <- function(pathways, orgDB, keyType = "ENTREZID") {
  for (i in names(pathways)) {
    if (typeof(pathways[[i]]) == "S4") {
      pathways[[i]] <- setReadable(pathways[[i]], OrgDb = orgDB, keyType = keyType)
    } else {
      pathways[[i]] <- data.frame()
    }
  }
  return(pathways)
}

#' Get Gene Set Enrichment Analysis (GSEA)
#'
#' This function performs gene set enrichment analysis (GSEA) on a given set of genes using the clusterProfiler and ReactomePA packages.
#' It supports analysis for human and mouse genes.
#'
#' @param res A data frame containing the gene data. It should have columns for the statistics and gene names.
#' @param maxGeneNum A numeric value specifying the maximum number of genes to include in the analysis. Default is 30000.
#' @param statCol A character string specifying the column name in res for the statistics. Default is "stat".
#' @param sigExpCol A character string specifying the column name in res for the significance expression. Default is "sig_exp".
#' @param organism A character string specifying the organism. It should be "human" or "mouse".
#' @return A list containing three elements: all - a list of data frames containing the GSEA results for all genes, sig - a
#' list of data frames containing the GSEA results for significant genes, summary - a data frame containing the summary of the GSEA results.
#' @examples
#' \dontrun{
#' GetGSEA(res, maxGeneNum = 30000, statCol = "stat", sigExpCol = "sig_exp", organism = "human")
#' }
#' @export
#' @importFrom clusterProfiler gseGO gseKEGG GSEA
#' @importFrom ReactomePA gsePathway
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom dplyr filter rename relocate
#' @importFrom rlang sym .data
#' @importFrom purrr map
#' @importFrom utils head tail
GetGSEA <- function(res, maxGeneNum = 30000, statCol = "stat", sigExpCol = "sig_exp", organism) {
  if (nrow(res) > maxGeneNum) {
    res <- rbind(utils::head(res, 15000), utils::tail(res, 15000))
  }
  allPathways <- GetPathwaysGSE(res, statCol = "stat", geneNameCol = "rownames", organism = organism)
  sigPathways <- GetPathwaysGSE(
    res |> dplyr::filter(!!sym(sigExpCol) != "NO"),
    statCol = "stat", geneNameCol = "rownames", organism = organism
  )
  gsePathways <- merge_nested_df(list(all = allPathways, sig = sigPathways))
  if (nrow(gsePathways) == 0) {
    stop("No pathways found for GSEA analysis.")
  }
  gsePathways <- dplyr::rename(gsePathways, Input = .data$Level1, DB = .data$Level2) |>
    dplyr::relocate(.data$DB, .data$Input, .before = "ID")
  return(list(all = allPathways, sig = sigPathways, summary = gsePathways))
}

#' Get Overrepresentation Analysis
#'
#' This function performs overrepresentation analysis (ORA) on a given set of genes using the clusterProfiler and ReactomePA packages.
#' It supports analysis for human and mouse genes.
#'
#' @param res A data frame containing the gene data. It should have columns for the statistics and gene names.
#' @param statCol A character string specifying the column name in res for the statistics. Default is "stat".
#' @param sigExpCol A character string specifying the column name in res for the significance expression. Default is "sig_exp".
#' @param organism A character string specifying the organism. It should be "human" or "mouse".
#' @return A list containing four elements: all - a list of data frames containing the ORA results for all genes,
#' up - a list of data frames containing the ORA results for up-regulated genes, down - a list of data frames containing
#' the ORA results for down-regulated genes, summary - a data frame containing the summary of the ORA results.
#' @examples
#' \dontrun{
#' GetOverRepAnalysis(res, statCol = "stat", sigExpCol = "sig_exp", organism = "human")
#' }
#' @export
#' @importFrom clusterProfiler enrichGO enrichKEGG enricher
#' @importFrom ReactomePA enrichPathway
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom dplyr filter rename relocate
#' @importFrom rlang sym .data
#' @importFrom purrr map
GetOverRepAnalysis <- function(res, statCol = "stat", sigExpCol = "sig_exp", organism) {
  allPathways <- GetPathwaysOverRep(
    res |> dplyr::filter(!!sym(sigExpCol) != "NO"),
    statCol = "stat", geneNameCol = "rownames", organism = organism
  )
  upPathways <- GetPathwaysOverRep(
    res |> dplyr::filter(!!sym(sigExpCol) == "UP"),
    statCol = "stat", geneNameCol = "rownames", organism = organism
  )
  downPathways <- GetPathwaysOverRep(
    res |> dplyr::filter(!!sym(sigExpCol) == "DOWN"),
    statCol = "stat", geneNameCol = "rownames", organism = organism
  )
  overRepPathways <- merge_nested_df(
    list(all = allPathways, up = upPathways, down = downPathways),
    columnsToDrop = c("category", "subcategory")
  )
  if (nrow(overRepPathways) == 0) {
    stop("No pathways found for OverRep analysis.")
  }
  overRepPathways <- dplyr::rename(overRepPathways, Input = .data$Level1, DB = .data$Level2) |>
    dplyr::relocate(.data$DB, .data$Input, .before = "ID")
  return(list(all = allPathways, up = upPathways, down = downPathways, summary = overRepPathways))
}

#' Get Over Representation Data Frame
#'
#' This function takes the output of an over-representation analysis and transforms it into a tidy data frame.
#' It handles missing values, converts certain columns to character, reshapes the data frame, and renames the columns.
#'
#' @param overRep A data frame containing the results of an over-representation analysis.
#'
#' @return A tidy data frame containing the over-representation analysis results.
#' @export
#' @importFrom dplyr case_when mutate_at select rename_all where contains vars all_of
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom stringr str_remove
#' @importFrom rlang .data
GetOverRepDf <- function(overRep) {
  if (nrow(overRep) == 0) {
    return(data.frame())
  }
  for (col in c("Sex", "Comparison", "Input")) {
    if (!col %in% colnames(overRep)) {
      stop(paste("The input dataframe does not contain the required column:", col))
    }
  }
  # relocate columns GeneRatio:Count to the end
  cols <- c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count", "geneID")
  overRep <- overRep |>
    dplyr::mutate(Sex = dplyr::case_when(is.na(.data$Sex) ~ "All", .default = .data$Sex)) |>
    dplyr::mutate_at(cols, as.character) |>
    tidyr::pivot_longer(cols = dplyr::all_of(cols), names_to = "metric", values_to = "value") |>
    dplyr::mutate(comparison_metric = paste(.data$Comparison, .data$metric, sep = "_")) |>
    dplyr::select(-c("Comparison", "metric"))
  colsBeforeWide <- colnames(overRep)
  overRep <- tidyr::pivot_wider(overRep, names_from = "comparison_metric", values_from = "value")
  colsToSelect <- setdiff(colnames(overRep), colsBeforeWide)
  overRep <- overRep |>
    tidyr::pivot_longer(cols = dplyr::all_of(colsToSelect), names_to = "metric", values_to = "value") |>
    dplyr::mutate(sex_metric = paste(.data$Sex, .data$metric, sep = "_")) |>
    dplyr::select(-c("Sex", "metric"))
  colsBeforeWide <- colnames(overRep)
  overRep <- tidyr::pivot_wider(overRep, names_from = "sex_metric", values_from = "value")
  colsToSelect <- setdiff(colnames(overRep), colsBeforeWide)
  overRep <- tidyr::pivot_longer(overRep, cols = dplyr::all_of(colsToSelect), names_to = "metric", values_to = "value") |>
    dplyr::mutate(input_metric = paste(.data$Input, .data$metric, sep = "_")) |>
    dplyr::select(-c("Input", "metric")) |>
    tidyr::pivot_wider(names_from = "input_metric", values_from = "value") |>
    dplyr::mutate_at(
      dplyr::vars(
        dplyr::contains("pvalue"), dplyr::contains("p.adjust"), dplyr::contains("qvalue"), dplyr::contains("Count"),
      ), as.numeric
    ) |>
    dplyr::select(-dplyr::where(~ all(is.na(.)))) |>
    dplyr::select(
      c("DB", "ID", "Description"), dplyr::contains("p.adjust"), dplyr::contains("Count"), dplyr::contains("geneID"),
      dplyr::contains("pvalue"), dplyr::contains("qvalue"), dplyr::contains("GeneRatio"), dplyr::contains("BgRatio")
    ) |>
    dplyr::rename_all(~ stringr::str_remove(.x, "byExpSex_")) |>
    dplyr::rename_all(~ stringr::str_remove(.x, "All_"))
  return(overRep)
}

#' Get GSE Data Frame
#'
#' This function takes a GSE data frame and performs various transformations on it.
#' It checks for required columns, handles missing values, converts certain columns to character,
#' reshapes the data frame, and renames the columns.
#'
#' @param gse A data frame containing GSE data.
#'
#' @return A transformed GSE data frame.
#' @export
#' @importFrom dplyr case_when mutate_at select rename_all where contains vars all_of
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom stringr str_remove
#' @importFrom rlang .data
GetGSEDf <- function(gse) {
  if (nrow(gse) == 0) {
    return(data.frame())
  }
  for (col in c("Sex", "Comparison")) {
    if (!col %in% colnames(gse)) {
      stop(paste("The input dataframe does not contain the required column:", col))
    }
  }
  if ("Input" %in% colnames(gse)) {
    stop("gse dataframe should not contain column named 'Input'")
  }
  cols <- c("setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
  gse <- gse |>
    dplyr::mutate(Sex = dplyr::case_when(is.na(.data$Sex) ~ "All", .default = .data$Sex)) |>
    dplyr::mutate_at(cols, as.character) |>
    tidyr::pivot_longer(cols = dplyr::all_of(cols), names_to = "metric", values_to = "value") |>
    dplyr::mutate(comparison_metric = paste(.data$Comparison, .data$metric, sep = "_")) |>
    dplyr::select(-c("Comparison", "metric"))
  colsBeforeWide <- colnames(gse)
  gse <- tidyr::pivot_wider(gse, names_from = "comparison_metric", values_from = "value")
  colsToSelect <- setdiff(colnames(gse), colsBeforeWide)
  gse <- tidyr::pivot_longer(gse, cols = dplyr::all_of(colsToSelect), names_to = "metric", values_to = "value") |>
    dplyr::mutate(sex_metric = paste(.data$Sex, .data$metric, sep = "_")) |>
    dplyr::select(-c("Sex", "metric")) |>
    tidyr::pivot_wider(names_from = "sex_metric", values_from = "value") |>
    dplyr::rename_all(~ stringr::str_remove(.x, "All_")) |>
    dplyr::mutate_at(
      dplyr::vars(
        dplyr::contains("setSize"), dplyr::contains("enrichmentScore"), dplyr::contains("NES"), dplyr::contains("pvalue"),
        dplyr::contains("p.adjust"), dplyr::contains("qvalue"), dplyr::contains("rank")
      ), as.numeric
    ) |>
    dplyr::select(-dplyr::where(~ all(is.na(.)))) |>
    dplyr::select(
      c("DB", "ID", "Description"), dplyr::contains("pvalue"), dplyr::contains("NES"), dplyr::contains("setSize"), dplyr::contains("core_enrichment"),
      dplyr::contains("enrichmentScore"), dplyr::contains("p.adjust"), dplyr::contains("qvalue"), dplyr::contains("rank"), dplyr::contains("leading_edge")
    ) |>
    dplyr::rename_all(~ stringr::str_remove(.x, "byExpSex_"))
  return(gse)
}
