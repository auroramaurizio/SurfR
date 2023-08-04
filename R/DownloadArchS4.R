#' DownloadArchS4 function
#'
#' Download count matrix from \url{https://maayanlab.cloud/archs4/},
#' given a vector of input GEO Sample accessions numbers (GSM).
#'
#' @param GSM Vector with the GSM ids of the samples to consider.
#' @param species Specify the specie of yuor GSM samples. Either human or mouse.
#' @param print_tsv Logical. If \code{TRUE}, outputs a tsv file with the count matrix. By default, FALSE.
#' @param filename Name of the tsv output file. Default is matrix.tsv.
#' @return A count matrix with gene on the row and GSM ID on the column.
#' @section Warning:
#' Add disclaim on GEO data curation.
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/geo} for info on GSM.
#' \url{https://maayanlab.cloud/archs4/} for info on ArchS4.
#' @family public-data functions
#' @export
#'
DownloadArchS4 <- function(GSM, species, print_tsv = FALSE, filename = NULL) {
  import::here(rhdf5, h5read)
  import::here(rhdf5, H5close)


  set.seed(42)
  options(timeout=600)
  matrixh5_url=paste0('https://mssm-seq-matrix.s3.amazonaws.com/',species,'_matrix.h5')

  series = h5read(matrixh5_url,  name = "meta/Sample_series_id", s3 = TRUE)
  samples = h5read(matrixh5_url,  name = "meta/Sample_geo_accession", s3 = TRUE)
  genes = h5read(matrixh5_url,  name = "meta/genes", s3 = TRUE)


  if (length(samples %in% GSM) == 0) {
    stop("The defined GSM ids do not have any match in ArchS4 database. \n We suggest to contact ArchS4 curator to add them.")
  }
  sample_locations = which(samples %in% GSM)

  # extract gene expression from compressed data
  expression = h5read(matrixh5_url, "data/expression",
                      index=list(1:length(genes), sample_locations), s3 = TRUE)
  H5close()

  rownames(expression) = genes
  colnames(expression) = samples[sample_locations]

  #-------- Print file --------
  if (print_tsv) {
    write.table(expression, file=filename, sep="\t", quote=FALSE)
  }

  return(as.data.frame(expression))
}
