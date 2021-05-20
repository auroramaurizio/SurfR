#' Prova import package function
#'
#' @param x a vector or a data.frame.
#' @param mean "geo" for geometric mean, "har" for harmonic mean.
#' @return The geometric or harmonic mean of \code{x}.
#' @examples
#' gm = prova(trees$Volume, mean = "geo")
#' hm = prova(trees$Volume, mean = "har")
#' @export
prova <- function(x, mean) {
  import::here(psych, geometric.mean, harmonic.mean)
  if (mean == "geo") {
    a = geometric.mean(x)
  } else if (mean == "har") {
    warning("are you sure you want to do harmonic mean?!")
    a = harmonic.mean(x)
  } else { stop("mean parameters is wrong!")}
  return(a)
}
