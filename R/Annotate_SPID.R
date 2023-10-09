#' Format a string using placeholders - function from hypeR
#'
#' @param string A an unformatted string with placeholders
#' @param ... Variables to format placeholders with
#' @return A formatted string
#'
#' @examples
#' \dontrun{
#' format_str("Format with {1} and {2}", "x", "y")
#' }
#'
#' @keywords internal
.format_str <- function(string, ...) {
  args <- list(...)
  for (i in seq_along(args)) {
    pattern <- paste("\\{", i, "}", sep="")
    replacement <- args[[i]]
    string <- gsub(pattern, replacement, string)
  }
  return(string)
}




#' Get url base for species-specific enrichr libraries  - function from hypeR
#'
#' @param db A species
#' @return A url
#'
#' @keywords internal
enrichr_urls <- function(db=c("Enrichr", "YeastEnrichr", "FlyEnrichr", "WormEnrichr", "FishEnrichr")) {
  switch(match.arg(db),
         "Enrichr"      = "http://maayanlab.cloud/Enrichr/{1}",
         "YeastEnrichr" = "http://maayanlab.cloud/YeastEnrichr/{1}",
         "FlyEnrichr"   = "http://maayanlab.cloud/FlyEnrichr/{1}",
         "WormEnrichr"  = "http://maayanlab.cloud/WormEnrichr/{1}",
         "FishEnrichr"  = "http://maayanlab.cloud/FishEnrichr/{1}"
  )
}

#' Connect to the enrichr web application  - function from hypeR
#'
#' @param endpoint The url endpoint to connect to
#' @param db A species
#' @return A web response
#'
#' @importFrom httr GET http_status
#'
#' @keywords internal
enrichr_connect <- function(endpoint, db=c("Enrichr", "YeastEnrichr", "FlyEnrichr", "WormEnrichr", "FishEnrichr")) {
  url <- enrichr_urls(db)
  response <- httr::GET(.format_str(url, endpoint))
  if (!http_status(response)$category == "Success") {
    stop(http_status(response)$message)
  }
  return(response)
}



#' Download data from enrichr in the form of a named list - function from hypeR
#'
#' @param genesets A name corresponding to available genesets
#' @param db A species
#' @return A list of genesets
#'
#' @examples
#' ATLAS <- enrichr_download("Human_Gene_Atlas")
#'
#' @importFrom httr content
#'
#' @export
enrichr_download <- function(genesets, db=c("Enrichr", "YeastEnrichr", "FlyEnrichr", "WormEnrichr", "FishEnrichr")) {
  response <- enrichr_connect(.format_str("geneSetLibrary?mode=text&libraryName={1}", genesets), db)
  data <- content(response, "text")
  split <- strsplit(data, split="\n")[[1]]
  #genesets <- sapply(split, function(x) strsplit(x, "\t")[[1]])
  for (i in seq_along(split)) {
    genesets[i] <- unlist(strsplit(split[[i]][1], "\t"))[1]
  }
  names(genesets) <- unlist(lapply(genesets, function(x) x[1]))
  lapply(genesets, function(x) {
    genes <- x[3:length(x)]
    genes <- genes[genes != ""]
    unique(genes)
  })
}






#' Annotate_SPID
#'
#' Annotate Surface Protein Coding genes according to EnrichR Libraries
#'
#' @param DGE Data.frame containing annotated DEG list, as the output of DGE or Gene2SProtein function.
#' @param enrich.database String containing the EnrichR databases you would like to consult. Default: WikiPathway_2021_Human.
#' @param output_tsv Logical. If \code{TRUE}, outputs a tsv file with the results. By default, FALSE.
#' @return A dataframe with surface protein coding \code{DEGs} annotation.
#' @examples
#' # Deseq2 output sample
#' DGE = data.frame(GeneID = c("DLK1", "TOP2A"),
#'                  Mean_CPM_T = c(5.92, 9.91),
#'                  Mean_CPM_C = c(0.04, 0.03),
#'                  log2FoldChange = c(10.22, 8.42),
#'                  lfcSE = c(0.80, 0.48),
#'                  stat = c(12.68, 17.69),
#'                  pvalue = c(7.30135e-37, 4.37011e-70),
#'                  padj = c(1.49936e-35, 1.12976e-67))
#' library(enrichR)
#' annotated_DGE = Annotate_SPID(DGE, "WikiPathway_2021_Human")
#'
#' # Output of Gene2SProtein function
#' GeneNames = c("CIITA", "EPCAM", "DLK1", "CD24", "TOP2A")
#' SurfaceProteins_df = Gene2SProtein(GeneNames, input_type = "gene_name")
#' annotated_SP = Annotate_SPID(SurfaceProteins_df, "GO_Biological_Process_2021")
#' @section Warning:
#' Be sure that enrich.database exist.
#' @family functional-annotation functions
#' @seealso \code{\link{DGE}} function for DGE,
#' and \code{\link{Gene2SProtein}} function for Gene2SProtein analysis
#' @importFrom enrichR listEnrichrDbs
#' @importFrom assertr col_concat
#' @importFrom tidyr separate_rows
#' @importFrom magrittr %>%
#' @importFrom utils write.table
#' @importFrom dplyr group_by summarise ungroup
#'
#' @export



Annotate_SPID <- function(DGE,
                          enrich.database  = "WikiPathway_2021_Human",
                          output_tsv = FALSE) {

  if (is.null(dim(DGE))) {
    stop("DGE is empty")
  }

  if (length(enrich.database) != 1) {
    stop("You can select only one enrich.database at a time.")
  }

  #db = enrichR::listEnrichrDbs()
  db = listEnrichrDbs()
  if (!(enrich.database %in% db$libraryName)) {
    stop(enrich.database, "is not a valid enrichR geneset.")
  }
  #annotation_table = hypeR::enrichr_download(enrich.database)
  annotation_table <- enrichr_download(enrich.database)

  # here gives a warning, but we can safely ignore it
  suppressWarnings({annotation_table <- as.data.frame(do.call(rbind, annotation_table))})

  colNames <- colnames(annotation_table) # could be any number of column names here
  annotation_table['test'] <- col_concat(annotation_table, sep = " ")
  annotation_table['GeneID'] <- trimws(annotation_table$test, which = c("both")) #remove whitespaces
  annotation_table$term <- row.names(annotation_table)
  annotation_table_sub <- annotation_table[,c("term","GeneID")]

  #exploded = unique(tidyr::separate_rows(annotation_table_sub, GeneID, sep = " ", convert = FALSE))
  exploded <- unique(separate_rows(annotation_table_sub, GeneID, sep = " ", convert = FALSE))
  #group by gene.. that is: one gene x row with associated all the descriptions (space separated)
  grouped  <-
    exploded %>%
    #dplyr::group_by(GeneID) %>%
    #dplyr::summarise(temp = toString(term)) %>%
    #dplyr::ungroup()
    group_by(GeneID) %>%
    summarise(temp = toString(term)) %>%
    ungroup()

  colnames(grouped) <- c("GeneID",enrich.database)

  merged <- merge(DGE, grouped, by="GeneID", all.x = TRUE)


  if (output_tsv) {
    write.table(merged, paste(enrich.database, "_SP_annotation.tsv", sep ="_"), quote = FALSE, sep = "\t")
  }
  GeneID <- term <- NULL
  return(merged)
}


