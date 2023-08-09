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
#' @importFrom hypeR enrichr_download
#' @importFrom assertr col_concat
#' @importFrom tidyr separate_rows
#' @importFrom magrittr %>%
#' @importFrom utils write.table
#' @importFrom dplyr group_by summarise ungroup
#'
#' @export

Annotate_SPID <- function(DGE,
                          enrich.database  = "WikiPathway_2021_Human",
                          output_tsv = F) {
  #import::here(enrichR)
  #import::here(hypeR)
  #import::here(assertr)
  #import::here(tidyr)
  #import::here(dplyr)
  #import::here(magrittr,"%>%")

  #enrichR::listEnrichrSites()

  if (is.null(dim(DGE))) {
    stop("DGE is empty")
  }

  if (length(enrich.database) != 1) {
    stop("You can select only one enrich.database at a time.")
  }

  #db = enrichR::listEnrichrDbs()
  db = listEnrichrDbs()
  if (!(enrich.database %in% db$libraryName)) {
    stop(paste(enrich.database, "is not a valid enrichR geneset."))
  }
  #annotation_table = hypeR::enrichr_download(enrich.database)
  annotation_table = enrichr_download(enrich.database)

  # here gives a warning, but we can safely ignore it
  suppressWarnings({annotation_table = as.data.frame(do.call(rbind, annotation_table))})

  colNames = colnames(annotation_table) # could be any number of column names here
  #annotation_table['test'] = assertr::col_concat(annotation_table, sep = " ")
  annotation_table['test'] = col_concat(annotation_table, sep = " ")
  annotation_table['GeneID'] <- trimws(annotation_table$test, which = c("both")) #remove whitespaces
  annotation_table$term = row.names(annotation_table)
  annotation_table_sub=annotation_table[,c("term","GeneID")]

  #exploded = unique(tidyr::separate_rows(annotation_table_sub, GeneID, sep = " ", convert = FALSE))
  exploded = unique(separate_rows(annotation_table_sub, GeneID, sep = " ", convert = FALSE))
  #group by gene.. that is: one gene x row with associated all the descriptions (space separated)
  grouped  <-
    exploded %>%
    #dplyr::group_by(GeneID) %>%
    #dplyr::summarise(temp = toString(term)) %>%
    #dplyr::ungroup()
    group_by(GeneID) %>%
    summarise(temp = toString(term)) %>%
    ungroup()

  colnames(grouped)=c("GeneID",enrich.database)

  merged <- merge(DGE, grouped, by="GeneID", all.x = TRUE)


  if (output_tsv) {
    write.table(merged, paste(enrich.database, "_SP_annotation.tsv", sep ="_"), quote = F, sep = "\t")
  }
  GeneID <- term <- NULL
  return(merged)
}


