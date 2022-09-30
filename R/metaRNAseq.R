#' metaRNAseq function
#'
#' Perform Meta-Analysis of RNA-Seq Data
#'
#' @param ind_deg List of indipendent DEG dataframes with p-values to be combined.
#' @param test_statistic p-value combination technique (inverse normal or Fisher):
#' \code{fishercomb}, \code{invnorm}.
#' By default: \code{fishercomb}.
#' @param BHth Benjamini Hochberg threshold.
#' @param nrep Vector of numbers of replicates used in each study to calculate the previous one-sided p-values. 
#' By default, the False Discovery Rate is controlled at 5%.
#' @return A list with \code{DEindices} of DEG at the chosen Benjamini Hochberg threshold, and 
#' \code{TestStatistic}, \code{rawpval}, \code{adjpval} vector for differential expression in the meta-analysis.
#' @examples
#' ind_deg = list(DEG1_df, DEG2_df, DEG3_df)
#' comb_pval_df = Meta(ind_deg, test_statistic = "invnorm", BHth = 0.05, nrep = c(2,2,2))
#' @section Warning:
#' Bla bla bla
#' @family aggregate functions
#' @seealso \code{\link{hello}} for DGE analysis, and \code{\link{hello}} for metaRNASeq package info
#' @export
metaRNAseq <- function(ind_deg,
                 test_statistic = "fishercomb",
                 BHth = 0.05,
                 nrep = NULL) {
  
  import::here(metaRNASeq)
  
  common_genes = Reduce(intersect, lapply(ind_deg, rownames))
  
  histp  <- list()
  rawpval <- list()
  for (i in 1:length(ind_deg)){
    ind_deg[[i]] <- ind_deg [[i]][common_genes,]
    rawpval[[i]] <- ind_deg[[i]][["pvalue"]]
    pdf(file = paste(names(ind_deg[i]),"_raw_pval_hist.pdf",sep ="", collapse = NULL))
    hist(rawpval[[i]], breaks=100, col="grey", main= names(ind_deg[i]), xlab="Raw p-values")
    dev.off()
  }
  
  ggplot2::ggsave(filename = "raw_pval.pdf",
                  plot = gridExtra::marrangeGrob(histp, nrow = 1, ncol = 1), 
                  device = "pdf")
  
  if (test_statistic == "fishercomb") {
  fish_comb <- fishercomb(rawpval, BHth = BHth)
  fish_comb$DEname = common_genes
  fish_comb$DE.fishercomb=ifelse(fish_comb$adjpval<=0.05,1,0)
  pdf(file = paste("fishercomb_pval_hist.pdf",sep ="", collapse = NULL))
  hist(fish_comb$rawpval, breaks=100, col="grey", main= names(ind_deg), xlab="Raw p-values")
  dev.off()
  return(fish_comb)
  } else if (test_statistic == "invnorm"){
  inv_norm <- invnorm(rawpval, nrep = nrep, BHth = BHth)
  inv_norm$DEname = common_genes
  inv_norm$DE.inv_norm=ifelse(inv_norm$adjpval<=0.05,1,0) 
  pdf(file = paste("invnorm_pval_hist.pdf",sep ="", collapse = NULL))
  hist(inv_norm$rawpval, breaks=100, col="grey", main= names(ind_deg), xlab="Raw p-values")
  dev.off()
  return(inv_norm)}} 


