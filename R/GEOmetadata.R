#' GEOmetadata function
#'
#' Download metadata from \url{https://www.ncbi.nlm.nih.gov/geo},
#' given an input GEO accession series.
#'
#' @param GSE The GSE series ID.
#' @param GPL The GPL series number.
#' Required only if the chosen GSE series ID include data from multiple sequencing platform
#' @return A dataframe with all the available characteristics in GEO metadata \code{genes} array.
#' @examples
#' metadata_GSE123345 = GEOmetadata(GSE = "GSE123345")
#' metadata_GSE123372_GPL17021 = GEOmetadata(GSE = "GSE123372", GPL = c("GPL17021"))
#' @section Warning:
#' If the GEO accession series has more than 1 sequencing platforms you need to
#' specify the GPL series numbers.
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/geo} for info on GEO repository
#' @export
GEOmetadata <- function(GSE, GPL = "") {
  import::here(stringr, str_remove_all)
  import::here(stringr, str_sub)
  import::here(stringr, str_split)
  import::here(stringr, str_replace_all)

  if (GPL != "") { file = paste(GSE,"-",GPL,"_series_matrix.txt.gz", sep='')
  } else { file = paste(GSE,"_series_matrix.txt.gz", sep='') }

  system(paste("wget https://ftp.ncbi.nlm.nih.gov/geo/series/",str_sub(GSE, start = 1, end = -4),"nnn/",
               GSE,"/matrix/",file, sep =''))

  metaHEAD = read.table(file, check.names = FALSE, sep = "\t", skipNul = TRUE, fill = TRUE)
  N_head = grep(pattern = "Sample_title", metaHEAD$V1)
  meta = read.table(file, skip = N_head-1 , check.names = FALSE, sep = "\t", skipNul = TRUE, fill = TRUE)
  GSMid = grep(pattern = "Sample_geo_accession", x = meta$V1)
  title = grep(pattern = "Sample_title", x = meta$V1)
  characteristic2select = grep(pattern = "Sample_characteristics", x = meta$V1)
  row2select = c(title, GSMid, characteristic2select)
  metadata_2_format = t(meta[row2select,])

  parentesis = grep(pattern = "\\(", metadata_2_format["V2",])
  metadata_2_format[,parentesis] = str_remove_all(string = metadata_2_format[,parentesis], pattern = "\\(")
  metadata_2_format[,parentesis] = str_remove_all(string = metadata_2_format[,parentesis], pattern = "\\)")

  metadata = as.data.frame(metadata_2_format[2:dim(metadata_2_format)[1],])
  metadata$GSE = GSE
  metadata = metadata[,c("GSE", colnames(metadata)[1:length(colnames(metadata))-1])]
  cnames = c("GSE","title","GSM")

  for (i in colnames(metadata)[4:length(colnames(metadata))]) {
    sep = unlist(str_split(metadata[,i], pattern = ": "))
    name = sep[1]
    info = sep[grep(pattern = name, sep, invert = T)]
    metadata[,i]  = info
    cnames = c(cnames, name)
  }
  colnames(metadata) = cnames
  row.names(metadata) = metadata$GSM
  colnames(metadata) = str_replace_all(string=colnames(metadata), pattern=" ", replacement= "_")


  stringr <- str_remove_all <- str_sub <- str_split <- str_split <- str_replace_all <- read.table <- NULL
  file.remove(file)
  return(metadata)

}
