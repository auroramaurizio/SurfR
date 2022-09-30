#' combine_fisher_invnorm
#'
#' Combine Meta-Analysis results with individual DE tables
#'
#' @param ind_deg List of indipendent DEG dataframes with p-values to be combined.
#' @param invnorm inverse normal p-value combination technique dataframe (output of metaRNAseq)
#' @param fishercomb Fisher p-value combination technique dataframe (output of metaRNAseq)
#' @param adjpval threshold to represent as binary the Meta-Analysis output adjpval.
#' @return A dataframe with \code{DEindices} and \code{DEname} of DEG at the chosen Benjamini Hochberg threshold, and 
#' \code{TestStatistic}, \code{rawpval}, \code{adjpval}, \code{binaryadjpval} vectors for differential expression in the meta-analysis.
#' @examples
#' ind_deg = list(DEG1_df, DEG2_df, DEG3_df)
#' names(ind_deg) = c("DEG1_df", "DEG2_df", "DEG3_df")
#' comb_pval_df = combine_fisher_invnorm(ind_deg, invnorm, fishercomb, adjpval = 0.05, output_tsv = T, output_filename = "combine_fisher_invnorm.tsv")
#' @section Warning:
#' Bla bla bla
#' @family aggregate functions
#' @seealso \code{\link{hello}} for DGE analysis, and \code{\link{hello}} for metaRNASeq package info
#' @export

combine_fisher_invnorm <- function(ind_deg,
                                   invnorm,
                                   fishercomb,
                                   adjpval = 0.05,
                                   output_tsv = T,
                                   output_filename = "combine_fisher_invnorm.tsv") {
  
                          #import::here(metaRNASeq)

                          DE = data.frame(genes=common_genes)
                          FC = data.frame(genes=common_genes)

                          for (i in 1:length(ind_deg)) {
                              ind_deg[[i]] <- ind_deg [[i]][common_genes,]
                              ind_deg[[i]][["binarypadj"]]=ifelse(ind_deg[[i]][["padj"]]<=adjpval,1,0)
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
                          
                            # -------- tsv --------
                          if (output_tsv) {
                          write.table(dgeResults, output_filename, quote = F, sep = "\t")
                          }
 
                          return(comb)}





