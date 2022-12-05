#' plotPCA function
#'
#' Plot PCA highlighting one or two data features
#'
#' @param matrix Count matrix with gene on the row and sample ID on the column.
#' @param metadata Sample metadata
#' @param nTOP number of top genes to use for principal components,
#' selected by highest row variance
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param color.by Name of one or more metadata columns to color cells by.
#' @param pt.shape If NULL, all points are circles \(default\).
#' You can specify any metadata attribute allowing for both different colors
#' and different shapes on cells.
#' @param cols.use Vector of colors, each color corresponds to an identity class.
#' By default, ggplot assigns colors.
#' @param label Logical. If \code{TRUE} adds samples label. Default = FALSE.
#' @param new.label If NULL, use the sample names as in metadata.
#' Otherwise you can specify new labels.
#' @param label.size Sets size of labels
#' @return PCA plot objec created by ggplot2,
#' which can be assigned and further customized.
#' @examples
#'
#' plotPCA(matrix = count_GSE133671,
#'         metadata = meta_GSE133671,
#'         color = "GSE",
#'         shape = "condition")
#' @export
plotPCA <- function(matrix,
                    metadata,
                    nTOP = 500,
                    dims = c(1,2),
                    color.by = NULL,
                    pt.shape = NULL,
                    cols.use = NULL,
                    label = FALSE,
                    new.label = NULL,
                    label.size = 5) {


  import::here(ggplot2)

  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }

  pca
  return(pca)
}
