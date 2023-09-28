#' Gene2SProtein function
#'
#' Detect Surface Proteins from a vector of genes.
#' The surface proteins are identified according to the in silico human surfaceome database,
#' proposed by Bernd Wollscheid group and available at \url{https://wlab.ethz.ch/surfaceome}
#'
#' @param genes A vector of genes.
#' @param input_type The gene identification type:
#' \code{gene_name}, \code{ensembl}, \code{entrez} or \code{uniProt_name}.
#' By default: \code{gene_name}.
#' @param output_tsv Logical. If \code{TRUE}, outputs a tsv file with the results. By default, FALSE.
#' @param output_filename Name of the tsv output file. Default is surfaceProteins.tsv.
#' @param Surfy_version  The version of surfy dataframe you wish to use. Choose between \code{log} or \code{newest}.
#' By default use the most recent \code{log} version.
#' If a log dataframe does not exist the \code{newest} is downloaded from \url{https://wlab.ethz.ch/surfaceome}.
#' @return A data frame with filtered surface proteins from the \code{genes} array.
#' The dataframe contains also addition information obtained from surfy.
#' @examples
#'  # from gene name IDs to Surface proteins
#'  GeneNames <- c("CIITA", "EPCAM", "DLK1", "CD24", "CDCP1", "LYVE1", "ABCD1", "VAMP1")
#'  SurfaceProteins_df <- Gene2SProtein(GeneNames, input_type = "gene_name")
#'
#'  # from ensembl IDs to Surface proteins
#'  Ensembl <- c("ENSG00000178343", "ENSG00000176895", "ENSG00000162419", "ENSG00000170776",
#'              "ENSG00000092529", "ENSG00000135926", "ENSG00000152595", "ENSG00000121577",
#'              "ENSG00000186094", "ENSG00000126773", "ENSG00000198918", "ENSG00000167378",
#'              "ENSG00000095574", "ENSG00000140678", "ENSG00000262484", "ENSG00000133739",
#'              "ENSG00000172469", "ENSG00000112992", "ENSG00000148343", "ENSG00000138593")
#'  SurfaceProteins_df <- Gene2SProtein(Ensembl, input_type = "ensembl",
#'                                    output_tsv = FALSE, Surfy_version = "new")
#' @section Warning:
#' The surfy database is interrogated using the gene identification type of your preference
#' between \code{gene_name}, \code{ensembl}, \code{entrez} or \code{uniProt_name}. Note that
#' you might loose some match due to difference gene version IDs.
#' @importFrom openxlsx read.xlsx write.xlsx
#' @importFrom utils write.table
#' @import BiocFileCache
#' @seealso \code{\link{DGE}} for DGE analysis,
#' \url{https://wlab.ethz.ch/surfaceome} for info on Surfy
#' @export

Gene2SProtein <- function(genes,
                          input_type = "gene_name",
                          output_tsv = FALSE,
                          output_filename = "surfaceProteins.tsv",
                          Surfy_version = "log") {

  now <- Sys.Date()
  log <- sort(list.files("log/", pattern = "*.xlsx"))

  # ---- surfy database versioning ----

  if (!Surfy_version %in% c("log", "new")) {
    stop("The specified database version is not available. \n  Choose between log or new")
  }

  if (length(log)>0 & Surfy_version == "log") {
    ST <- read.xlsx(xlsxFile = paste0(".log/", log[length(log)]),
                    startRow = 1)
  } else {
    dir.create(".log/", recursive = TRUE, showWarnings = FALSE)
    surfaceome_table_url='https://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx'


    bfc <- BiocFileCache(ask = FALSE)
    path <- bfcrpath(bfc, surfaceome_table_url)

    # fix the reading function
    ST <- read.xlsx(xlsxFile = path,
                    sheet = 1,
                    startRow = 2)


    write.xlsx(ST, paste0(".log/","table_S3_surfaceome_",now,".xlsx"), overwrite = TRUE)
    # -------------------------------------------------------------------------
  }

  # ---- input identification type ----

  ## check gene identification type
  if (!input_type %in% c("gene_name", "entrez", "ensembl", "uniProt_name")) {
    stop("The specified input type is not available. \n Choose between gene_name, entrez, ensembl or uniProt_name")
  }
  type <- if(input_type == "gene_name"){"UniProt.gene"
  } else if(input_type == "entrez"){"GeneID"
  } else if(input_type == "ensembl"){"Ensembl.gene"
  } else if (input_type == "uniProt_name"){"UniProt.name"}

  # ---- filter out NA value from database ----
  # GeneID and UniProt.name

  ST <- ST[!is.na(ST$Surfaceome.Label),]
  ST <- ST[!is.na(ST$UniProt.name),]
  ST <- ST[!is.na(ST[,type]),]
  row.names(ST) <- ST$UniProt.name

  # ---- filter the database for the input genes and checks ----
  proteins <- ST[ST[,type] %in% genes,]

  # check size proteins data.frame
  if (dim(proteins)[1] == 0) {

    warning("The input genes ",
            "do not have any match ",
            "in the surfaceome database. \n ",
            "Check gene alias alias and input type! ")
    surface.proteins <- data.frame(matrix(nrow = 0, ncol = length(colnames(ST))))
    colnames(surface.proteins)  <- colnames(ST)
  } else {
    # -------- filter surface proteins and checks --------
    surface.proteins <- proteins[proteins$Surfaceome.Label=='surface',]

    if (dim(surface.proteins)[1] == 0) {

      warning("No surface proteins were found among your list of genes")
    } else {

      message(dim(surface.proteins)[1], " out of ", length(genes), " genes ",
              "have a matching surface protein")
    }

  }

  # rename columns for label consistency
  surface.proteins$entrezID <- surface.proteins$GeneID
  surface.proteins$GeneID <- surface.proteins$UniProt.gene

  # -------- tsv --------
  if (output_tsv) {
    write.table(surface.proteins, output_filename, quote = FALSE, sep = "\t")
  }

  openxlsx <- read.xlsx <- write.xlsx <- write.table <- NULL
  return(surface.proteins)
}
