#' metaRNAseq function
#'
#' Perform Meta-Analysis of RNA-Seq Data
#'
#' @param ind_deg List of indipendent named DEG dataframes with p-values to be combined.
#' @param test_statistic p-value combination technique (inverse normal or Fisher):
#' \code{fishercomb}, \code{invnorm}.
#' By default: \code{fishercomb}.
#' @param BHth Benjamini Hochberg threshold.
#' @param nrep Vector of numbers of replicates used in each study to calculate the previous one-sided p-values.
#' @param adjpval.t threshold to represent as binary the Meta-Analysis output adjpval.
#' @param plot Logical. If TRUE plot histogram of pvalues.
#' By default, the False Discovery Rate is controlled at 0.05.
#' @return A list with \code{DEindices} of DEG at the chosen Benjamini Hochberg threshold, and
#' \code{TestStatistic}, \code{rawpval}, \code{adjpval}, \code{binaryadjpval} vectors for differential expression in the meta-analysis.
#' @examples
#' # Deseq2 output samples
#' DGE1 <- data.frame(GeneID = c("DLK1", "EPCAM"),
#'                  Mean_CPM_T = c(5.92, 9.91),
#'                  Mean_CPM_C = c(0.04, 0.03),
#'                  log2FoldChange = c(10.22, 8.42),
#'                  lfcSE = c(0.80, 0.48),
#'                  stat = c(12.68, 17.69),
#'                  pvalue = c(7.30135e-37, 4.37011e-70),
#'                  padj = c(1.49936e-35, 1.12976e-67),
#'                  row.names = c("DLK1", "EPCAM"))
#' DGE2 <- data.frame(GeneID = c("DLK1", "EPCAM"),
#'                  Mean_CPM_T = c(3.92, 8.91),
#'                  Mean_CPM_C = c(0.04, 0.03),
#'                  log2FoldChange = c(7.22, 5.81),
#'                  lfcSE = c(0.80, 0.48),
#'                  stat = c(12.68, 17.69),
#'                  pvalue = c(7.30135e-37, 4.37011e-70),
#'                  padj = c(1.49936e-35, 1.12976e-67),
#'                  row.names = c("DLK1", "EPCAM"))
#' # input list
#' ind_deg <- list(DEG1_df = DGE1, DEG2_df = DGE2)
#' # perform meta-analysis
#' comb_pval_df <- metaRNAseq(ind_deg, test_statistic = "invnorm", BHth = 0.05, nrep = c(2,2))
#' @family meta-analysis functions
#' @seealso \code{\link{DGE}} for DGE analysis,
#' and \url{https://cran.r-project.org/web/packages/metaRNASeq/vignettes/metaRNASeq.pdf}
#' for metaRNASeq package info
#' @importFrom metaRNASeq fishercomb invnorm
#' @importFrom gridExtra marrangeGrob
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics hist
#' @importFrom ggplot2 ggsave
#' @export


metaRNAseq <- function(ind_deg,
                       test_statistic = "fishercomb",
                       BHth = 0.05,
                       adjpval.t = 0.05,
                       nrep = NULL,
                       plot = FALSE) {

  # Check if ind_deg is a list of at least two data.frame
  if (!is.list(ind_deg)) {
    stop("ind_deg is not a list. Please provide a list of at least two data.frames")
  }
  if (length(ind_deg) < 2) {
    stop("ind_deg contains ", length(ind_deg), " data.frame. Please provide a list of at least two data.frames")
  }

  common_genes <- Reduce(intersect, lapply(ind_deg, rownames))
  if (length(common_genes) == 0) {
    stop(" your DGE data.frames do not have any shared gene")
  }

  # if test_statistic is invnorm check that nrep has the same number of elements than ind_deg
  if (test_statistic == "invnorm" && length(ind_deg) != length(nrep)) {
    stop("nrep must have the same number of elements of ind_deg.")
  }

  histp  <- list()
  rawpval <- list()

  for (i in seq_along(ind_deg)){
    if (!("pvalue" %in% colnames(ind_deg[[i]]))) {
      stop("The DGE dataframe ", names(ind_deg)[i], " must include a column named pvalue.")

    }
    ind_deg[[i]] <- ind_deg[[i]][common_genes, ]
    rawpval[[i]] <- ind_deg[[i]][["pvalue"]]
    if (plot) {
      pdf(file = paste(names(ind_deg[i]), "_raw_pval_hist.pdf", sep = "", collapse = NULL))
      hist(rawpval[[i]], breaks = 100, col = "grey", main = names(ind_deg[i]), xlab = "Raw p-values")
      dev.off()
    }
  }

  if (plot) {
    ggsave(filename = "raw_pval.pdf",
           plot = marrangeGrob(histp, nrow = 1, ncol = 1),
           device = "pdf")
  }

  if (test_statistic == "fishercomb") {
    fish_comb <- fishercomb(rawpval, BHth = BHth)
    fish_comb$DEname <- common_genes
    fish_comb$binaryadjpval <- ifelse(fish_comb$adjpval <= adjpval.t, 1, 0)
    if (plot) {
      pdf(file = paste("fishercomb_pval_hist.pdf", sep = "", collapse = NULL))
      hist(fish_comb$rawpval, breaks = 100, col = "grey",
           main = names(ind_deg), xlab = "Raw p-values")
      dev.off()
    }
    return(fish_comb)
  } else if (test_statistic == "invnorm") {
    inv_norm <- invnorm(rawpval, nrep = nrep, BHth = BHth)
    inv_norm$DEname <- common_genes
    inv_norm$binaryadjpval <- ifelse(inv_norm$adjpval <= adjpval.t, 1, 0)
    if (plot) {
      pdf(file = paste("invnorm_pval_hist.pdf", sep = "", collapse = NULL))
      hist(inv_norm$rawpval, breaks = 100, col = "grey", main = names(ind_deg),
           xlab = "Raw p-values")
      dev.off()
    }
    return(inv_norm)
  } else {
    stop("The defined test_statistic could only be either fishercomb or invnorm.")
  }

}
