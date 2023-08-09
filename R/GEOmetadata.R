#' GEOmetadata function
#'
#' Download metadata from \url{https://www.ncbi.nlm.nih.gov/geo},
#' given an input GEO accession series.
#'
#' @param GSE The GSE series ID.
#' @param GPL The GPL series numbers.
#' Required only if the chosen GSE series ID include data from multiple sequencing platforms.
#' @return A dataframe with all the available characteristics in GEO metadata \code{genes} array.
#' @examples
#' # only one sequencing platform
#' mGSE133671 = GEOmetadata(GSE = "GSE133671")
#' # multiple sequencing platforms
#' mGSE59483 = GEOmetadata("GSE59483", GPL = c("GPL11154", "GPL15520"))
#' @section Warning:
#' If the GEO accession series has more than 1 sequencing platforms you need to
#' specify the GPL series numbers.
#' Add disclaim on GEO data curation!
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/geo} for info on GEO repository
#' @family public-data functions
#' @importFrom stringr str_remove_all str_sub str_replace_all str_split
#' @export
GEOmetadata <- function(GSE, GPL = "") {
  #import::here(stringr, str_remove_all)
  #import::here(stringr, str_sub)
  #import::here(stringr, str_split)
  #import::here(stringr, str_replace_all)


  metafromfile <- function(file) {
    metaHEAD = read.table(file, check.names = FALSE, sep = "\t", skipNul = TRUE, fill = TRUE)
    N_head = grep(pattern = "Sample_title", metaHEAD$V1)
    meta = read.table(file, skip = N_head-1 , check.names = FALSE, sep = "\t", skipNul = TRUE, fill = TRUE)
    GSMid    = grep(pattern = "Sample_geo_accession",   x = meta$V1)
    title    = grep(pattern = "Sample_title",           x = meta$V1)
    GPLid    = grep(pattern = "Sample_platform_id",     x = meta$V1)
    c2select = grep(pattern = "Sample_characteristics", x = meta$V1)
    row2select = c(title, GSMid, GPLid, c2select)
    metadata_2_format = t(meta[row2select,])

    parentesis = grep(pattern = "\\(", metadata_2_format["V2",])
    metadata_2_format[,parentesis] = str_remove_all(string = metadata_2_format[,parentesis], pattern = "\\(")
    metadata_2_format[,parentesis] = str_remove_all(string = metadata_2_format[,parentesis], pattern = "\\)")

    metadata = as.data.frame(metadata_2_format[2:dim(metadata_2_format)[1],])
    metadata$GSE = GSE
    metadata = metadata[,c("GSE", colnames(metadata)[1:length(colnames(metadata))-1])]
    cnames = c("GSE","title","GSM", "GPL")

    # replace empty spaces with NA
    metadata[metadata==''] <- NA


    for (i in 5:length(colnames(metadata))) {

      item = unlist(lapply(str_split(metadata[,i], pattern = ": "), `[[`, 1))
      # check if the characteristic is the same
      if (length(unique(item[!is.na(item)]))==1) {
        name = unique(item[!is.na(item)])
      } else {
        name = "multi_characteristic"
        warning(paste0(i,"th ",
                      "characteristic column in ", GSE,
                      " metadata presents multiple variables: ",
                      paste(unique(item[!is.na(item)]), collapse = ", ")))
      }
      # check for NA
      metadata[,i][is.na(metadata[,i])] <- paste(name, ": NA")
      metadata[,i]  = unlist(lapply(str_split(metadata[,i], pattern = ": "), `[[`, 2))
      cnames = c(cnames, name)
      }

    colnames(metadata) = cnames
    row.names(metadata) = metadata$GSM
    colnames(metadata) = str_replace_all(string=colnames(metadata), pattern=" ", replacement= "_")
    return(metadata)
  }

  if (length(GPL) == 1 && GPL == "") {
    # only one sequencing platform
    file = paste(GSE,"_series_matrix.txt.gz", sep='')
    system(paste("wget https://ftp.ncbi.nlm.nih.gov/geo/series/",str_sub(GSE, start = 1, end = -4),"nnn/",
                 GSE,"/matrix/",file, sep =''))

    metadata = metafromfile(file)
    file.remove(file)

  } else {
    #multiple sequencing platforms
    metadata_i = list()
    for (gpl in GPL) {
      file = paste(GSE,"-",gpl,"_series_matrix.txt.gz", sep='')
      system(paste("wget https://ftp.ncbi.nlm.nih.gov/geo/series/",str_sub(GSE, start = 1, end = -4),"nnn/",
                   GSE,"/matrix/",file, sep =''))
      metadata_i[[gpl]] =  metafromfile(file)
      file.remove(file)
      }
    metadata = Reduce(f = rbind, metadata_i)
    }

  stringr <- str_remove_all <- str_sub <- str_split <- str_replace_all <- read.table <- NULL
  return(metadata)

}

