##' Gene enrichment using Enrichr
##'
##' Gene enrichment using Enrichr, slighthly modified by Aurora Maurizio.
##' @title Gene enrichment using Enrichr
##' @param genes (Required). Character vector of Entrez gene symbols as input. A data.frame
##' of gene symbols in first column is also acceptable, optionally a score denoting the
##' degree of membership between 0 and 1 in the second column.
##' @param databases (Required). Character vector of databases to search.
##' See https://maayanlab.cloud/Enrichr/ for available databases.
##' @param background (Optional). Character vector of Entrez gene symbols to be used as
##' background. A data.frame of gene symbols in first column is also acceptable.
##' Default is \code{"NULL"}. Enrichment analysis with background genes is only available
##' on the main site (Enrichr). Also, it is using a different API service (Speedrichr),
##' hence it is a little slower to complete and return the results.
##' @param include_overlap (Optional). Download database in GMT format to include 'Overlap'
##' in the resulting data.frame when analysing with a background. Default is \code{"FALSE"}.
##' @return Returns a list of data.frame of enrichment terms, p-values, ...
##' @author Wajid Jawaid \email{wj241@alumni.cam.ac.uk}
##' @importFrom httr POST
##' @importFrom httr use_proxy
##' @importFrom rjson fromJSON
##' @importFrom utils read.table
##' @export
##' @examples
##' data(input) # Load example input genes
##' data(background) # Load example background genes
##' dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
##'          "GO_Biological_Process_2023")
##' if (getOption("enrichR.live")) {
##'   enriched1 <- enrichr(input, dbs)
##'   print(head(enriched1[[1]]))
##' }

enrichr <- function(genes, databases) {
  if (length(genes) < 1) {
    stop("No genes have been given")
  }

  base.address <- getOption("enrichR.base.address")
  getEnrichr(url = base.address)

  if (!getOption("enrichR.live")) {
    stop("Enrichr website is unreachable")
  }

  if (is.null(databases)) {
    stop("No databases have been provided")
  }

  if (is.vector(genes) & !all(genes == "") & length(genes) != 0) {
    temp <- POST(url=paste0(getOption("enrichR.base.address"), "enrich"),
                 body=list(list=paste(genes, collapse="\n")))
  } else if (is.data.frame(genes)) {
    temp <- POST(url=paste0(getOption("enrichR.base.address"), "enrich"),
                 body=list(list=paste(paste(genes[,1], genes[,2], sep=","),
                                      collapse="\n")))
  } else {
    warning("genes must be a non-empty vector of gene names or a data.frame with genes and score.")
  }

  dbs <- as.list(databases)
  result <- lapply(dbs, function(x) {
    r <- getEnrichr(url = paste0(base.address, "export"), query = list(file = "API", backgroundType = x))
    if (!getOption("enrichR.live")) stop("Enrichr website is unreachable")
    r <- gsub("&#39;", "'", intToUtf8(r$content))
    tc <- textConnection(r)
    r <- read.table(tc, sep = "\t", header = TRUE, quote = "", comment.char = "")
    close(tc)
    return(r)
  })

  names(result) <- databases
  return(result)
}



#' Enrichment function
#'
#' Perform enrichment Analysis of RNA-Seq Data
#'
#' @param dfList Dataframes list
#' @param enrich.databases Vector of EnrichR databases to consult
#' @param p_adj Double. Adjusted pvalue threshold for the enrichment
#' @param logFC Double. Fold change threshold for the enrichment
#' @param save.results Logical. If TRUE saves input gene lists and enrichment results.
#' @return A list of enrichment tables for upregulated and downregulated genes
#' in the different enrichr databases
#' @examples
#' \dontrun{
#' df1 <- data.frame(GeneID  = c("MEST", "CDK1", "PCLAF", "BIRC5"),
#'                   baseMean = c(13490.22, 10490.23, 8888.33, 750.33),
#'                   log2FoldChange = c(5.78, 6.76, -7.78, -8.78),
#'                   padj = c(2.28e-143, 2.18e-115, 2.18e-45, 0.006),
#'                   row.names = c("MEST", "CDK1", "PCLAF", "BIRC5"))
#' df2 <- data.frame(GeneID  = c("MEST", "CDK1", "PCLAF", "BIRC5"),
#'                   baseMean = c(13490.22, 10490.23, 8888.33, 750.33),
#'                   log2FoldChange = c(5.78, 6.76, -7.78, -8.78),
#'                   padj = c(2.28e-143, 2.18e-115, 2.18e-45, 0.006),
#'                   row.names = c("MEST", "CDK1", "PCLAF", "BIRC5"))
#' dfList <- list(df1 = df1, df2 = df2)
#' test <- Enrichment(dfList, enrich.databases = c("GO_Cellular_Component_2021"),
#'                    save.results = FALSE)}
#' @family functional-annotation functions
#' @seealso \url{https://maayanlab.cloud/Enrichr/} for additional information about enrichR.
#' @importFrom enrichR listEnrichrDbs enrichr setEnrichrSite
#' @importFrom openxlsx write.xlsx
#' @importFrom utils write.table
#' @export


Enrichment <- function(dfList, enrich.databases  = c("GO_Biological_Process_2021",
                                                     "GO_Cellular_Component_2021",
                                                     "GO_Molecular_Function_2021",
                                                     "KEGG_2021_Human",
                                                     "MSigDB_Hallmark_2020",
                                                     "WikiPathways_2016",
                                                     "BioCarta_2016",
                                                     "Jensen_TISSUES",
                                                     "Jensen_COMPARTMENTS",
                                                     "Jensen_DISEASES"),
                       p_adj = 0.05, logFC = 1,
                       save.results = FALSE) {


  websiteLive <- getOption("enrichR.live", default = FALSE)

  if (websiteLive) {
    setEnrichrSite("Enrichr") # Human genes
    db <- listEnrichrDbs()
  } else {
    stop("enrichR website can not be reached at the moment. Please,
          check your internet connection and retry later.")
  }

  enrichr.list <- list()

  if (length(setdiff(enrich.databases, db$libraryName)) > 0) {
    warning(setdiff(enrich.databases, db$libraryName), " is not an enrichR geneset and will be removed.\n")
    enrich.databases <- intersect(enrich.databases, db$libraryName)
  }

  if (length(enrich.databases) == 0) {
    stop("Please provide at least one valid enrich.database.")
  }

  for (i in names(dfList)) {
    df_obj <- dfList[[i]]
    signif <- (df_obj[df_obj$padj <= p_adj, ])
    number_of_sig_genes  <- nrow(signif)

    message(i, " ", number_of_sig_genes, " significant genes\n")

    if (number_of_sig_genes == 0) {
      stop("no significant genes found. Enrichment can't be performed.")
    }

    neg <- nrow(signif[signif$log2FoldChange < logFC, ])

    message(i, " ", neg, " negative fold change\n")

    neg_list <- rownames(signif[signif$log2FoldChange < logFC, ])

    if (length(neg_list) == 0) {
      warning("There are no significantly downregulated genes in ", i)
  } else {
    if (save.results) {
      dir.create("enrichR/", showWarnings = FALSE, recursive = TRUE)
      write.table(neg_list, paste("./enrichR/FDRdown_", i,
                                  ".txt", sep = ""), quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }

    }

    pos  <- nrow(signif[signif$log2FoldChange > logFC, ])
    message(i, " ", pos, " positive fold change\n")

    pos_list  <- rownames(signif[signif$log2FoldChange > logFC, ])

    if (length(pos_list) == 0) {
      warning("There are no significantly upregulated genes in ", i)
    } else {
      if (save.results) {
      dir.create("enrichR/", showWarnings = FALSE, recursive = TRUE)
      write.table(pos_list, paste("./enrichR/FDRup_", i,
                                  ".txt", sep = ""), quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      }

    }

    enrichr.list[[i]] <- lapply(list(pos_list, neg_list), function(x) {
      enrichr(genes = x, databases = enrich.databases)
    })
    names(enrichr.list[[i]]) <-  c("fdr_up", "fdr_down")

  }

  if (save.results) {
    dir.create("enrichR/", showWarnings = FALSE, recursive = TRUE)
    for (i in names(dfList)) {
      for (j in c("fdr_up", "fdr_down")){
        filename <- paste("./enrichR/", i, j, ".xlsx", sep = "")
        if (!is.null(enrichr.list[[i]][[j]])) {
          write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
        }
      }
    }
  }

  return(enrichr.list)
}
