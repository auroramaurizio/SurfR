#' metaRNAseq function
#'
#' Perform Meta-Analysis of RNA-Seq Data
#'
#' @param ind_deg List of indipendent DEG dataframes with p-values to be combined.
#' @param test_statistic p-value combination technique (inverse normal or Fisher):
#' \code{fishercomb}, \code{invnorm}.
#' By default: \code{fishercomb}.
#' @param BHth Benjamini Hochberg threshold. 
#' By default, the False Discovery Rate is controlled at 5%.
#' @return A list with \code{DEindices} of DEG at the chosen Benjamini Hochberg threshold, and 
#' \code{TestStatistic}, \code{rawpval}, \code{adjpval} vector for differential expression in the meta-analysis.
#' @examples
#' ind_deg = list(DEG1_df, DEG2_df, DEG3_df)
#' comb_pval_df = Meta(ind_deg, test_statistic = "invnorm", BHth = 0.05)
#' @section Warning:
#' Bla bla bla
#' @family aggregate functions
#' @seealso \code{\link{hello}} for DGE analysis, and \code{\link{hello}} for metaRNASeq package info
#' @export


Meta <- function(ind_deg,
                 test_statistic = "fishercomb",
                 BHth = 0.05) {
  
  import::here(metaRNASeq)
  
  common_genes = row.names(Reduce(intersect, ind_deg ))
  
  rawpval <- list()
  nrep <-c()
  for (i in 1:length(ind_deg)){
    ind_deg[[i]] <- ind_deg [[i]][common_genes,]
    rawpval[[i]] <- ind_deg[[i]][["pvalue"]]
    nrep[i] <- length(ind_deg[[i]])
  }
  
  if (test_statistic == "fishercomb") {
  fish_comb <- fishercomb(rawpval, BHth = BHth)
  return(fish_comb)
  } else if (test_statistic == "invnorm"){
  # nrep = Vector of numbers of replicates used in each study to calculate the previous one-sided p-values.
  inv_norm <- invnorm(rawpval, nrep = nrep, BHth = BHth) 
  return(inv_norm)} }
  


