#' SVenn
#'
#' Venn diagram of common surface proteins overexpressed among up to 7 different studies
#'
#' @param S_list A list of a maximum of 7 surface protein sets detected in different studies.
#' @param cols.use Vector of colors, each color corresponds to a study.
#' By default, ggplot assigns colors.
#' @param opacity Degree of opacity for the color(s) specified with cols.use (less opacity, more transparency).
#' @param output_intersectionFile logical. If TRUE (default) write an xlsx output of protein in the intersections.
#' @param filename Name of the output file with the intersections.
#' @return venn plot of common genes.
#' @examples
#' S_list <- list(SP1 = c("EPCAM", "CD24",  "DLK1",  "CDCP1", "LYVE1"),
#'               SP2 = c("DLK1", "EPCAM", "EGFR", "UPK1A", "UPK2"))
#' SP <- SVenn(S_list, output_intersectionFile = FALSE)
#' @family plot functions
#' @importFrom venn venn
#' @importFrom openxlsx write.xlsx
#' @importFrom scales hue_pal
#' @seealso \code{\link{Gene2SProtein}} for detection of Surface proteins from a list of genes.
#' @export


SVenn <- function(S_list, cols.use = NULL,
                  opacity = 0.5,
                  output_intersectionFile = TRUE,
                  filename = "intersection.xlsx") {

  if (length(S_list) > 7) {
    stop("This function can plot Venn diagram with up to 7 sets")
  }

  if (is.null(x = cols.use)) {
    cols.use <- hue_pal()(length(x = names(S_list)))
  } else if (length(cols.use) < length(x = names(S_list))) {

    stop("you have", length(x = names(S_list)),
         "unique elements and supplied only", length(cols.use), "color \n")
  }


    venn(S_list,
         opacity = opacity,
         box = FALSE,
         zcolor = cols.use,
         ilcs = 1.5,
         sncs = 1.5,
         ggplot = FALSE)

    # Table down pathway intersection
    if (output_intersectionFile) {
      list_intersection <- attr(x = venn(S_list,
                                         intersections = TRUE,
                                         opacity = opacity,
                                         box = FALSE,
                                         zcolor = cols.use,
                                         ilcs = 1.5,
                                         sncs = 1.5,
                                         ggplot = FALSE),
                                         "intersections")

      df <- data.frame()
      for (intersection in names(list_intersection)) {
        entry <- data.frame(Intersection = intersection,
                            SurfaceProteins = list_intersection[[intersection]])
        df <- rbind(df, entry)
      }
      write.xlsx(df, filename, asTable = TRUE, overwrite = TRUE)
    }
}
