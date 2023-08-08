#' DGE function
#'
#' Perform Differential Gene Expression Analysis of RNA-Seq Data
#'
#' @param expression Dataframe with counts
#' @param metadata Dataframe with sample metadata
#' @param Nreplica  Double. Minimum number of replicates in each group
#' @param design Design formula for DGE
#' @param condition Column of the metadata ti use for DGE results
#' @param TEST Character. sample name in metadata
#' @param CTRL Character. sample name in metadata
#' @param FC_filt Dataframe with counts
#' @param alpha Double. the significance cutoff used for optimizing the independent filtering (by default 0.1).
#' If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param output_tsv Logical. If \code{TRUE}, outputs a tsv file with the results. By default, FALSE.
#' @param output_filename Name of the tsv output file. Default is DEGs.tsv.
#' @return A dataframe with \code{DEGs}
#' @examples
#' # Simulation of bulk RNA data
#' countData <- matrix(floor(runif(10000, min=0, max=101)),ncol=4)
#' colnames(countData) <- paste0("sample", 1:4)
#' rownames(countData) <- paste0("gene", 1:(10000/4))
#' metadata = data.frame(samplesID = paste0("sample", 1:4),
#'                      condition = factor(c("A","A","B","B")))
#' row.names(metadata) <- metadata$samplesID
#' # Perform DGE
#' DGEresults = DGE(expression = countData, metadata = metadata,
#'                  Nreplica = 2,
#'                  design = "~condition",condition = "condition",
#'                  TEST = "A", CTRL = "B")
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts results
#' @importFrom edgeR cpm
#' @importFrom utils write.table
#' @export

DGE <- function(expression,
                metadata,
                Nreplica,
                design = "~condition",
                condition = "condition",
                TEST,
                CTRL,
                alpha = 0.05,
                FC_filt = 0,
                output_tsv = FALSE,
                output_filename = "DEGs.tsv") {

  #import::here(DESeq2, DESeqDataSetFromMatrix)
  #import::here(DESeq2, DESeq)
  #import::here(DESeq2, counts)
  #import::here(DESeq2, results)

  #import::here(edgeR)


  # check the match between TEST CTRL and metadata
  if (length(c(CTRL, TEST) %in% metadata[,condition]) == 0) {
    stop(paste0("In your metadata the column ",condition, "does not contain ",  CTRL, " and ", TEST,"."))
  } else if (length(c(CTRL) %in% metadata[,condition]) == 0) {
    stop(paste0("In your metadata the column ",condition, "does not contain ",  CTRL, "."))
  } else if (length(c(TEST) %in% metadata[,condition]) == 0) {
    stop(paste0("In your metadata the column ",condition, "does not contain ",  TEST,"."))
  }

  if (!(setequal(row.names(metadata), colnames(expression)))) {
    stop("row.names of metadata must be equal to expression matrix colnames")
  }

  # perform DGE
  design.formula = as.formula(design)
  dds <- DESeqDataSetFromMatrix(countData = expression,
                                colData = metadata,
                                design = design.formula)

  min.samples = Nreplica
  #keep <- rowSums(edgeR::cpm(counts(dds)) >= 1) >= min.samples
  keep <- rowSums(cpm(counts(dds)) >= 1) >= min.samples
  dds <- dds[keep,]
  #cpm = edgeR::cpm(counts(dds), log = F)
  cpm = cpm(counts(dds), log = F)
  dga <- DESeq(object = dds,
               test = "Wald",
               fitType = "parametric",
               betaPrior = FALSE,
               minReplicatesForReplace = Inf)


  dgeResults <- results(dga,
                      contrast = c(condition,TEST,CTRL),
                      cooksCutoff          = Inf,
                      independentFiltering = TRUE,
                      alpha                = alpha,
                      pAdjustMethod        = "BH")

  dgeResults = dgeResults[order(dgeResults$pvalue, decreasing = F),]



  dgeResults$Mean_CPM_C = rowMeans(cpm[row.names(dgeResults),
                                       row.names(metadata[metadata[,condition] ==CTRL,])])
  dgeResults$Mean_CPM_T = rowMeans(cpm[row.names(dgeResults),
                                       row.names(metadata[metadata[,condition] ==TEST,])])



  dgeResults$GeneID = row.names(dgeResults)
  dgeResults = dgeResults[order(dgeResults$log2FoldChange, decreasing = T),
                          c('GeneID', 'Mean_CPM_T','Mean_CPM_C',
                            'log2FoldChange','lfcSE','stat','pvalue','padj')]

  # -------- tsv --------
  if (output_tsv) {
    write.table(dgeResults, output_filename, quote = F, sep = "\t")
  }

  return(dgeResults) }



