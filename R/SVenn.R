#' SVenn
#'
#' Venn diagram of common surface proteins overexpressed among different studies
#'
#' @param S_list A list of surface proteins detected in different studies.
#' @param cols.use Vector of colors, each color corresponds to a study.
#' By default, ggplot assigns colors.
#' @param main Plot title. Default = Venn.
#' @param opacity Degree of opacity for the color(s) specified with cols.use (less opacity, more transparency).
#' @param output_intersection logical. If TRUE write an output of protein in the intersections
#' @param filename Name of the output file with the intersections
#' @examples
#' SVenn(S_list, main = "my Venn")
#' @return a
#' @family aggregate functions
#' @seealso \code{\link{hello}} for counts data and metadata download, and \code{\link{hello}} for Gene2SProtein analysis
#' @export
#'

SVenn <- function(S_List,cols.use=NULL,
                  main = "Venn",
                  opacity = 0.5,
                  output_intersection = FALSE,
                  filename = "intersection.xlsx") {
  import::venn
  import::ggplot2

  if (is.null(x = cols.use)) {
    cols.use = scales::hue_pal()(length(x = names(S_list)) )
  } else if (length(cols.use) < length(x = names(S_list))) {
    stop(paste("you have",length(x = names(S_list)),
               "unique elements and supplied only",length(cols.use),"color \n"))
  }

  SP <- venn::venn(S_list,
                   opacity = opacity,
                   box = FALSE,
                   ilab = TRUE,
                   zcolor = cols.use,
                   ilcs = 1.5,
                   sncs = 1.5,
                   ggplot = FALSE)

  SP <- SP + ggplot2::ggtitle(main)


  # ------ Table down pathway intersection --------
  if (output_intersection) {
    a = venn::venn(S_list, intersections = T)
    list_intersection = attr(x = a, "intersections")

    df = data.frame()
    for (intersection in names(list_intersection)) {
      entry = data.frame(Intersection = intersection,
                         SurfaceProteins = list_intersection[[intersection]])
      df = rbind(df, entry)
    }
    write.xlsx(df, filename, asTable = T)
  }


}
