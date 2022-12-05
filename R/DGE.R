#' DGE function
#'
#' Perform Differential Gene Expression Analysis of RNA-Seq Data
#'
#' @param expression Dataframe with counts 
#' @param metadata Dataframe with sample metadata
#' @param TEST Character. sample name in metadata
#' @param CTRL Character. sample name in metadata  
#' @param FC_filt Dataframe with counts 
#' @param Nreplica  Double. Minimum number of replicates in each group
#' @param alpha Double. the significance cutoff used for optimizing the independent filtering (by default 0.1). 
#' If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param output_tsv Logical. If \code{TRUE}, outputs a tsv file with the results. By default, FALSE.
#' @param output_filename Name of the tsv output file. Default is DEGs.tsv.
#' @return A dataframe with \code{DEGs}  
#' @examples
#' DEG1_df = DGE(expression, metadata, TEST = "T", CTRL = "N", Nreplica = 45, alpha = 0.05, output_tsv = T, output_filename = "DEG1.tsv")
#' @section Warning:
#' Bla bla bla
#' @family aggregate functions
#' @seealso \code{\link{hello}} for counts data and metadata download, and \code{\link{hello}} for Gene2SProtein analysis
#' @export
DGE <- function(expression,
                metadata,
                Nreplica,
                TEST,
                CTRL,
                alpha = 0.1, 
                exp_filt_ctrl = 0,
                FC_filt = 0,
                output_tsv = T,
                output_filename = "DEGs.tsv") {
  
  import::here(DESeq2)
  import::here(edgeR)
  
  fCountsData <- expression[,match(metadata_ArchS4$samples, colnames(expression))] 
  # perform DGE
  design.formula = as.formula("~condition+series")
  dds <- DESeqDataSetFromMatrix(countData = fCountsData,
                                colData = metadata,
                                design = design.formula)
  
  min.samples = Nreplica
  keep <- rowSums(edgeR::cpm(counts(dds)) >= 1) >= min.samples
  dds <- dds[keep,]
  cpm = edgeR::cpm(counts(dds), log = F)
  dga <- DESeq(object = dds, 
               test = "Wald", 
               fitType = "parametric", 
               betaPrior = FALSE,
               minReplicatesForReplace = Inf)
  
  
  dgeResults <- results(dga, 
                      contrast = c("condition","T","N"), 
                      cooksCutoff          = Inf,
                      independentFiltering = TRUE,
                      alpha                = alpha,
                      pAdjustMethod        = "BH")
  
  dgeResults = dgeResults[order(dgeResults$pvalue, decreasing = F),]
  
  
  dgeResults$Mean_CPM_C = rowMeans(cpm[row.names(dgeResults), 
                                       row.names(metadata[metadata$condition==CTRL,])])
  dgeResults$Mean_CPM_T = rowMeans(cpm[row.names(dgeResults), 
                                       row.names(metadata[metadata$condition==TEST,])])
  

  
  dgeResults$GeneID = row.names(dgeResults)
  dgeResults = dgeResults[order(dgeResults$log2FoldChange, decreasing = T),
                          c('GeneID', 'Mean_CPM_T','Mean_CPM_C',
                            'log2FoldChange','lfcSE','stat','pvalue','padj')]
  
  # -------- tsv --------
  if (output_tsv) {
    write.table(dgeResults, output_filename, quote = F, sep = "\t")
  }
  
  return(dgeResults) }
  
  


