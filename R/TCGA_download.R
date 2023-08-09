#' TCGA_download function
#'
#' Downloads count matrix data from TCGA_download
#'
#' @param project, Character. A valid project from TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param whichcounts Character. Counts data form to use. Choose from: unstranded, stranded_first,stranded_second. By default, unstranded.
#' @param barcodes Character. A vector with names of the barcodes you want to download. If NULL (default) it downloads all the available barcodes in the project.
#' @param save.matrix Logical. If \code{TRUE}, outputs a tsv file with the Matrix. By default, FALSE.
#' @param save.metadata Logical. If \code{TRUE}, outputs a tsv file with the metadata. By default, FALSE.
#' @return A list containing the Matrix and the metadata.
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare
#' @importFrom SummarizedExperiment assay
#' @importFrom utils write.table
#' @import biomaRt
#' @examples
#' GBM_list_s1 <- TCGA_download(project="TCGA-GBM",
#'                              whichcounts = "unstranded",
#'                              save.matrix = FALSE, save.metadata = FALSE,
#'                              barcodes = c("TCGA-06-0878-01A-01R-1849-01"))
#' # remove downloaded data from TCGA
#' unlink('GDCdata', recursive = TRUE, force = TRUE)
#' file.remove("MANIFEST.txt")
#' @family public-data functions
#' @export
TCGA_download <- function(project,
                 whichcounts = "unstranded",
                 save.matrix = FALSE,
                 save.metadata = FALSE,
                 barcodes = NULL) {

      #import::here(TCGAbiolinks)
      #import::here(data.table)
      #import::here(biomaRt)
      #import::here(SummarizedExperiment)

      if (is.null(barcodes)) {
      #   query <- TCGAbiolinks::GDCquery(project,
         query <- GDCquery(project,
                                         data.category= "Transcriptome Profiling",
                                         data.type = "Gene Expression Quantification",
                                         sample.type = c("Solid Tissue Normal","Primary Tumor"),
                                         workflow.type="STAR - Counts")
       } else {
      #   query <- TCGAbiolinks::GDCquery(project,
         query <- GDCquery(project,
                                         data.category= "Transcriptome Profiling",
                                         data.type = "Gene Expression Quantification",
                                         sample.type = c("Solid Tissue Normal","Primary Tumor"),
                                         workflow.type="STAR - Counts",
                                         barcode = barcodes)
       }



       #tryCatch(TCGAbiolinks::GDCdownload(query, method = "client"), error = function(e)TCGAbiolinks::GDCdownload(query))
       #data <- TCGAbiolinks::GDCprepare(query, summarizedExperiment = TRUE)
       tryCatch(GDCdownload(query, method = "client"), error = function(e)GDCdownload(query))
       data <- GDCprepare(query, summarizedExperiment = TRUE)

       dir = 'TCGA/'

       # Save matrix
       #Matrix <- SummarizedExperiment::assay(data, whichcounts)
       Matrix <- assay(data, whichcounts)
       row.names(Matrix) <- data@rowRanges$gene_name

       # Save metadata
       metadata = as.data.frame(data@colData)
       #annotation_table = subset(metadata, select = -c(treatments, primary_site, disease_type))

       if (save.matrix ) {
         dir.create(dir, recursive = TRUE, showWarnings = FALSE)
         write.table(Matrix, file=paste0(dir, project,"_expression.tsv"), sep="\t", quote=FALSE)
       }

       if (save.metadata ) {
         dir.create(dir, recursive = TRUE, showWarnings = FALSE)
         write.table(metadata, file=paste0(dir, project,"_metadata.tsv",sep=""),sep="\t",quote=FALSE)
       }

       mat.met = list(Matrix, metadata)

       #print(TCGAbiolinks:::getProjectSummary(project))

       return(mat.met) }










