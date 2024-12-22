##' Set Enrichr Website
##'
##' Set Enrichr Website
##' @title Set Enrichr Website
##' @param site site requested
##' @return Changes Enrichr Website connection
##' @author Alexander Blume
##' @export
setEnrichrSite <- function(site) {
  enrichR.base.address <- "https://maayanlab.cloud/"
  site <- gsub(enrichR.base.address, "", site)
  matched <- grep(paste0("^",site),
                  getOption("enrichR.sites"),
                  ignore.case = TRUE, value = FALSE)
  if( length(matched) == 0 ) {
    message("Given website does not match available sites: ", site)
    message("Choose from:\n",
            paste("-",getOption("enrichR.sites"), "\n"))
  } else if (length(matched) > 1) {
    message("Given website matches multiple options: ", site)
    message(paste("-", getOption("enrichR.sites")[matched], "\n"),)
  } else {
    site <- getOption("enrichR.sites")[matched]

    options(enrichR.base.address = paste0(enrichR.base.address,site,"/"))
    message("Connection changed to ",paste0(enrichR.base.address,site,"/"))
    getEnrichr(url = paste0(enrichR.base.address,"datasetStatistics"))
  }
}


##' Helper function
##'
##' Helper function for HTTP methods GET
##' @title Check if EnrichR is online
##' @param url (Required). URL address to check
##' @param ... (Optional). Additional parameters to pass to GET
##' @return Invisible response object from GET request if successful, or NULL on failure
##' @importFrom httr GET
##' @importFrom httr status_code
##' @importFrom httr http_status
getEnrichr <- function(url, ...) {
  options(enrichR.live = FALSE) # Default to offline
  tryCatch({
    # Make the GET request
    response <- GET(url = url, ...)
    # Check the status code
    if (status_code(response) == 200) {
      options(enrichR.live = TRUE) # Set to live if the request succeeds
      invisible(response)
    } else {
      message(http_status(response)$message)
      NULL
    }
  },
  warning = function(w) {
    message("Warning: ", w$message)
    NULL
  },
  error = function(e) {
    message("Error: ", e$message)
    NULL
  })
}


##' Gene enrichment using Enrichr
##'
##' Gene enrichment using Enrichr, slighthly modified by Aurora Maurizio.
##' @title Gene enrichment using Enrichr
##' @param genes (Required). Character vector of Entrez gene symbols as input. A data.frame
##' of gene symbols in first column is also acceptable, optionally a score denoting the
##' degree of membership between 0 and 1 in the second column.
##' @param databases (Required). Character vector of databases to search.
##' See https://maayanlab.cloud/Enrichr/ for available databases.
##' @return Returns a list of data.frame of enrichment terms, p-values, ...
##' @author Wajid Jawaid \email{wj241@alumni.cam.ac.uk}
##' @importFrom httr POST
##' @importFrom httr use_proxy
##' @importFrom rjson fromJSON
##' @importFrom utils read.table
##' @export
##' @examples
##' \dontrun{
##' GeneID  = c("MEST", "CDK1", "PCLAF", "BIRC5")
##' dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
##'          "GO_Biological_Process_2023")
##' enriched1 <- enrichr(GeneID, dbs)
##' print(head(enriched1[[1]]))}


enrichr <- function(genes, databases) {
  if (length(genes) < 1) {
    stop("No genes have been given")
  }

  base.address <- "https://maayanlab.cloud/Enrichr/"
  enrichR.base.address <- "https://maayanlab.cloud/Enrichr/"
  getEnrichr(url = base.address)

  if (!getOption("enrichR.live")) {
    stop("Enrichr website is unreachable")
  }

  if (is.null(databases)) {
    stop("No databases have been provided")
  }

  if (is.vector(genes) & !all(genes == "") & length(genes) != 0) {
    temp <- POST(url=paste0(enrichR.base.address, "enrich"),
                 body=list(list=paste(genes, collapse="\n")))
  } else if (is.data.frame(genes)) {
    temp <- POST(url=paste0(enrichR.base.address, "enrich"),
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

##' Set Enrichr Website
##'
##' Set Enrichr Website
##' @title Set Enrichr Website
##' @param site site requested
##' @return Changes Enrichr Website connection
##' @author Alexander Blume
##' @export
setEnrichrSite <- function(site) {
  enrichR.sites.base.address = "https://maayanlab.cloud/"
  site <- gsub(enrichR.sites.base.address, "", site)
  matched <- grep(paste0("^",site),
                  getOption("enrichR.sites"),
                  ignore.case = TRUE, value = FALSE)
  if( length(matched) == 0 ) {
    message("Given website does not match available sites: ", site)
    message("Choose from:\n",
            paste("-",getOption("enrichR.sites"), "\n"))
  } else if (length(matched) > 1) {
    message("Given website matches multiple options: ", site)
    message(paste("-", getOption("enrichR.sites")[matched], "\n"),)
  } else {
    site <- getOption("enrichR.sites")[matched]

    options(enrichR.base.address = paste0("https://maayanlab.cloud/Enrichr/",site,"/"))
    message("Connection changed to ",paste0("https://maayanlab.cloud/Enrichr/",site,"/"))
    getEnrichr(url = paste0(getOption("enrichR.base.address"),"datasetStatistics"))
  }
}


##' Look up available databases on Enrichr
##'
##' Look up available databases on Enrichr
##' @title Look up available databases on Enrichr
##' @return A data.frame of available Enrichr databases
##' @author Wajid Jawaid \email{wj241@alumni.cam.ac.uk}
##' @importFrom rjson fromJSON
##' @export
##' @examples
##' dbs <- listEnrichrDbs()
listEnrichrDbs <- function() {
  dfSAF <- getOption("stringsAsFactors", FALSE)
  options(stringsAsFactors = FALSE)
  dbs <- getEnrichr(url = "https://maayanlab.cloud/Enrichr/datasetStatistics")
  if (!getOption("enrichR.live")) return()
  dbs <- intToUtf8(dbs$content)
  dbs <- fromJSON(dbs)
  dbs <- lapply(dbs$statistics, function(x) do.call(cbind.data.frame, x))
  dbs <- do.call(rbind.data.frame, dbs)
  options(stringsAsFactors = dfSAF)
  dbs
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
