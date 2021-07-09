#' Detect Surface Proteins from a list of genes
#'
#' @param genes array of gene names.
#' @param output_tsv output a tsv file with the results - boolean variable.
#' @return A dataframe of surface proteins filtered from the \code{genes} array.
#' @examples
#' GeneNames = c("CIITA", "EPCAM", "DLK1", "CD24", "CDCP1", "LYVE1", "ABCD1", "VAMP1")
#' SurfaceProteins_df = Gene2SProtein(GeneNames, output_tsv = FALSE)
#' @export
Gene2SProtein <- function(genes, output_tsv) {
  surfaceome_table_url='https://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx'
  import::here(rio, import)
  ST = import(file = surfaceome_table_url,
              which = 1,
              readxl = FALSE,
              startRow = 2)
  # filter NA value in GeneID and UniProt.name
  ST = ST[!is.na(ST$UniProt.name),]
  ST = ST[!is.na(ST$GeneID),]
  row.names(ST) <- ST$UniProt.name
  proteins = ST[ST$UniProt.gene %in% genes,]
  proteins = proteins[!is.na(proteins$SS),]
  print(table(proteins$Surfaceome.Label))
  surface.proteins = proteins[proteins$Surfaceome.Label=='surface',]
  if (output_tsv) {
    write.table(surface.proteins, "surfaceProteins.tsv", quote = F, sep = "\t")
  }
  return(surface.proteins)
}
