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
#' dfList = list(df1 = df1, df2 = df2)
#' test = Enrichment(dfList, save.results = FALSE)
#' @family functional-annotation functions
#' @seealso \url{https://maayanlab.cloud/Enrichr/} for additional information about enrichR.
#' @importFrom enrichR listEnrichrDbs enrichr
#' @importFrom openxlsx write.xlsx
#' @importFrom utils write.table
#' @export


Enrichment <- function(dfList ,enrich.databases  = c("GO_Biological_Process_2021",
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

  #import::here(enrichR)
  #import::here(openxlsx)

  enrichr.list <- list()
  # -------------------------
  # enrichment Parameters
  # -------------------------

  #db = enrichR::listEnrichrDbs()
  db = listEnrichrDbs()
  if (length(setdiff(enrich.databases, db$libraryName))>0) {
    warning(paste(setdiff(enrich.databases, db$libraryName), "is not an enrichR geneset and will be removed.\n"))
    enrich.databases = intersect(enrich.databases, db$libraryName)
  }

  if (length(enrich.databases)==0) {
    stop("Please provide at least one valid enrich.database.")
  }


  for (i in names(dfList)) {
      df_obj <- dfList[[i]]
      signif <- (df_obj[df_obj$padj <= p_adj, ])
  number_of_sig_genes  <- nrow(signif)
  #print(head(signif))
  cat(i, number_of_sig_genes, "significant genes\n")

  if (number_of_sig_genes == 0 ) {
    stop("no significant genes found. Enrichment can't be performed.")
  }


  neg <- nrow(signif[signif$log2FoldChange < logFC, ])
  cat(i, neg, "negative fold change\n")

  neg_list <-  rownames(signif[signif$log2FoldChange < logFC, ])
  if (length(neg_list)==0) {
    warning(paste("There are no significantly downregulated genes in",i))
  } else {
    if (save.results) {
      dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
      write.table(neg_list, paste('./enrichR/FDRdown_',i,
                                  '.txt', sep =''), quote = F,
                  row.names = F, col.names = F)
    }

  }



  pos  <- nrow(signif[signif$log2FoldChange > logFC, ])
  cat(i, pos, "positive fold change\n")

  pos_list  <- rownames(signif[signif$log2FoldChange > logFC, ])

  if (length(pos_list)==0) {
    warning(paste("There are no significantly upregulated genes in",i))
  } else {
    if (save.results) {
      dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
      write.table(pos_list, paste('./enrichR/FDRup_',i,
                                  '.txt', sep =''), quote = F,
                  row.names = F, col.names = F)
    }

  }



  enrichr.list[[i]] <- lapply(list(pos_list,neg_list),function(x) {
    #enrichR::enrichr(genes = x, databases = enrich.databases)
    enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down")

  }

  if (save.results) {
    dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
    for (i in names(dfList)) {
      for (j in c("fdr_up","fdr_down")){
        filename = paste("./enrichR/",i,j,".xlsx", sep="")
        if(!is.null(enrichr.list[[i]][[j]])) {
          #openxlsx::write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
          write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
        }
      }}
  }


  return(enrichr.list)
  }

