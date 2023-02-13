#' Annotate_SPID
#'
#' Annotate Surface Protein Coding genes according to EnrichR Libraries
#'
#' @param DGE Dataframe containing SP annotated DEG list
#' @param enrich.database String containg the EnrichR database you would like to consult. Default WikiPathway_2021_Human.
#' @param output_tsv Logical. If \code{TRUE}, outputs a tsv file with the results. By default, FALSE.
#' @return A dataframe with surface protein coding \code{DEGs} annotation
#' @examples
#' annotated = Annotate_SPID(DGE, "OMIM_Disease")
#' @section Warning:
#' Bla bla bla
#' @family aggregate functions
#' @seealso \code{\link{hello}} for DGE, and \code{\link{hello}} for Gene2SProtein analysis
#' @export

Annotate_SPID <- function(DGE, enrich.database  = "WikiPathway_2021_Human", output_tsv = F) {
  # options(warn=-1)
  import::here(hypeR)
  import::here(assertr)
  import::here(MASS)
  import::here(tidyr)
  import::here(dplyr)
  import::here(magrittr,"%>%")

  annoname = enrich.database
  annotation_table = hypeR::enrichr_download(enrich.database)
  annotation_table = as.data.frame(do.call(rbind, annotation_table)) #here gives a worning ignore it
  colNames = colnames(annotation_table) # could be any number of column names here
  annotation_table['test']=assertr::col_concat(annotation_table, sep = " ")
  annotation_table['GeneID'] <- trimws(annotation_table$test, which = c("both")) #remove whitespaces
  annotation_table$term = row.names(annotation_table)
  annotation_table_sub=annotation_table[,c("term","GeneID")]
  colnames(annotation_table_sub)
  exploded = unique(tidyr::separate_rows(annotation_table_sub, GeneID, sep = " ", convert = FALSE))
  #group by gene.. that is: one gene x row with associated all the descriptions (space separated)
  grouped  <-
    exploded %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarise(temp = toString(term)) %>%
    dplyr::ungroup()

  nrow(grouped)
  colnames(grouped)=c("GeneID",annoname)

  merged <- merge(DGE, grouped, by="GeneID", all.x = TRUE)

  if (output_tsv) {
    write.table(merged, paste(annoname, "_SP_annotation.tsv", sep ="_"), quote = F, sep = "\t")
  }
  return(merged)
}


