#' TCGA_download function
#'
#' Downloads count matrix data from TCGA_download
#'
#' @param project, Character. A valid project from TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param whichcounts Character. Counts data form to use. Choose from: unstranded, stranded_first,stranded_second. By default, unstranded.
#' @param save.matrix Logical. If \code{TRUE}, outputs a tsv file with the Matrix. By default, FALSE.
#' @param save.metdata Logical. If \code{TRUE}, outputs a tsv file with the metadata. By default, FALSE.
#' @return A list containing the Matrix and the metadata. 
#' @examples
#' TCGA_download(project="TCGA-UVM", whichcounts = "unstranded", save.matrix = F, save.metadata = T) 
#' @section Warning:
#' Bla bla bla
#' @family aggregate functions
#' @seealso \code{\link{hello}} for \code{\link{hello}} for DGE analysis.
#' @export


TCGA_download <- function(project, 
                 whichcounts = "unstranded",
                 save.matrix = F,
                 save.metadata = F) {
  
      import::here(TCGAbiolinks)
      import::here(data.table)
      import::here(biomaRt)


       query <- GDCquery(project, 
                         data.category="Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification",
                         sample.type = c("Solid Tissue Normal","Primary Tumor"),
                         workflow.type="STAR - Counts")
       
       tryCatch(GDCdownload(query, method = "client"),error = function(e)GDCdownload(query))

       data <- GDCprepare(query)

       dir = 'TCGA/'
       dir.create(dir)
       setwd(dir)

       
       # Save matrix
       Matrix <- assay(data, whichcounts)
       row.names(Matrix) <- data@rowRanges$gene_name

       # Save metadata
       metadata = as.data.frame(data@colData)
       #annotation_table = subset(metadata, select = -c(treatments, primary_site, disease_type))
    
       if (save.matrix ) {
       write.table(Matrix, file=paste(project,"_expression.tsv"), sep="\t", quote=FALSE)
       }
    
       if (save.metadata ) {
       write.table(metadata, file=paste(project,"_metadata.tsv",sep=""),sep="\t",quote=F)
       }
    
       mat.met = list(Matrix, metadata)
       
       #print(TCGAbiolinks:::getProjectSummary(project))
    
       return(mat.met) }
 








    
