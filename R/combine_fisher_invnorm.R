#' combine_fisher_invnorm
#'
#' Combine Meta-Analysis results with individual DE tables
#'
#' @param ind_deg List of indipendent DEG dataframes with p-values to be combined.
#' @param invnorm inverse normal p-value combination technique dataframe (output of metaRNAseq)
#' @param fishercomb Fisher p-value combination technique dataframe (output of metaRNAseq)
#' @param adjpval threshold to represent as binary the Meta-Analysis output adjpval.
#' @param output_tsv logical. If TRUE, it outputs table with results. Default: TRUE
#' @param output_filename File name for the results file.
#' @return A dataframe with \code{DEindices} and \code{DEname} of DEG at the chosen Benjamini Hochberg threshold, and
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
#' # perform invnorm meta-analysis
#' invnorm <- metaRNAseq(ind_deg, test_statistic = "invnorm", BHth = 0.05, nrep = c(2,2))
#' # perform fishercomb meta-analysis
#' fishercomb <- metaRNAseq(ind_deg, test_statistic = "fishercomb", BHth = 0.05)
#' # combine results
#' comb_pval_df <- combine_fisher_invnorm(ind_deg,
#'                                       invnorm, fishercomb,
#'                                       adjpval = 0.05,
#'                                       output_tsv = FALSE)
#' @family meta-analysis functions
#' @seealso \code{\link{DGE}} function for DGE analysis,
#' and \url{https://cran.r-project.org/web/packages/metaRNASeq/vignettes/metaRNASeq.pdf}
#' for metaRNASeq package info
#' @importFrom dplyr relocate
#' @importFrom magrittr %>%
#' @importFrom utils write.table
#' @export

combine_fisher_invnorm <- function(ind_deg,
                                   invnorm,
                                   fishercomb,
                                   adjpval = 0.05,
                                   output_tsv = TRUE,
                                   output_filename = "combine_fisher_invnorm.tsv") {

  #import::here(dplyr)
  #import::here(magrittr,"%>%")

  # Check if ind_deg is a list of at least two data.frame

  if (!is.list(ind_deg)) {
    stop("ind_deg is not a list. Please provide a list of at least two data.frames")
  }
  if(length(ind_deg) <2 ) {
    stop("ind_deg contains", length(ind_deg), "data.frame. Please provide a list of at least 2 data.frames")
  }

  common_genes <- Reduce(intersect, lapply(ind_deg, rownames))
  # check if common_genes is empty
  if (length(common_genes) == 0) {
    stop(" your DGE data.frames do not have common genes")
  }

  DE <- data.frame(genes=common_genes)
  FC <- data.frame(genes=common_genes)

  for (i in seq_along(ind_deg)) {
    ind_deg[[i]] <- ind_deg [[i]][common_genes,]
    ind_deg[[i]][["binarypadj"]] <- ifelse(ind_deg[[i]][["padj"]]<=adjpval,1,0)
    FC[[paste(names(ind_deg[i]), "log2FC", sep ="_")]] <- ind_deg[[i]][["log2FoldChange"]]
    DE[[paste(names(ind_deg[i]), "binarypadj", sep ="_")]] <- ind_deg[[i]][["binarypadj"]]
  }

  FC <- FC[-c(1) ]
  DE <- DE[-c(1) ]

  signsFC <- mapply(FC, FUN=function(x) sign(x))
  sumsigns <- apply(signsFC,1,sum)
  commonsgnFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns),0)
  FC$signFC <- commonsgnFC

  comb <- cbind(FC,DE)

  comb$DE_fishercomb <- invnorm$binaryadjpval
  comb$DE_invnorm <- fishercomb$binaryadjpval
  comb$GeneID <- common_genes

  #comb <- comb %>%  dplyr::relocate(GeneID)
  comb <- comb %>% relocate(GeneID)



  # -------- tsv --------
  if (output_tsv) {
    write.table(comb, output_filename, quote = FALSE, sep = "\t")
  }

  GeneID <- NULL

  return(comb)

}





