#' Enrichment_barplot
#'
#' Barplot representing the top up-regulated or down-regulated significant pathways
#'
#' @param Enrich A list of enrichment tables for up and down-regulated genes in the different enrichR databases.
#' Output of Enrichment.R function for one DGE experiment.
#' @param enrich.databases Vector of EnrichR databases to consider. These databases must be present in the Enrich list.
#' @param p_adj Double. Minimum Adjusted pvalue threshold for the enrichment
#' @param num_term Double. Number of up-regulated and dw-regulated terms to represent
#' @param cond String. Title of the plot.
#' @param plot Logical. If TRUE save plot as pdf.
#' @return bar plot of significant pathways.
#' @examples
#' dbs <- c("GO_Biological_Process_2021","WikiPathways_2016", "MSigDB_Hallmark_2020")
#' dfList <- list()
#' if (requireNamespace("enrichR", quietly = TRUE)) {
#'     up_genes <- c("RUNX1", "DLK1", "TOP2A", "EPCAM", "GATA1", "KDR")
#'     dfList[["fdr_up"]] <- enrichR::enrichr(up_genes, dbs)
#'     dw_genes <- c("CD275", "COL1A1", "COL1A2","LUM", "SOX9")
#'     dfList[["fdr_down"]] <- enrichR::enrichr(dw_genes, dbs)
#'     # Plot for upregulated gene
#'     Enrichment_barplot(dfList,
#'                        enrich.databases = dbs,
#'                        p_adj = 0.01, num_term = 5, cond = "UP")
#'     # Plot for downregulated gene
#'     Enrichment_barplot(dfList,
#'                        enrich.databases = dbs,
#'                        p_adj = 0.01, num_term = 5, cond = "DOWN")
#' } else {
#'     print("example require enrichR package")
#' }
#' @family functional-annotation functions
#' @family plot functions
#' @import knitr
#' @importFrom stringr str_replace
#' @importFrom ggplot2 ggplot geom_bar aes scale_fill_gradient scale_x_discrete labs coord_flip
#' @importFrom utils head
#' @importFrom grDevices pdf dev.off
#' @export
#'


Enrichment_barplot <- function(Enrich,
                               enrich.databases  = c("GO_Biological_Process_2021",
                                                     "GO_Cellular_Component_2021",
                                                     "GO_Molecular_Function_2021"),
                               p_adj = 0.05, num_term = 10, cond = "UP",
                               plot = FALSE) {


     #import::here(stringr)
     #import::here(ggplot2)

     enrichR.table <- data.frame()

     # check if Enrich contains enrich.databases


     if (cond == "UP") {
       enrich_list <- Enrich[["fdr_up"]]

     } else if (cond == "DOWN") {
       enrich_list <- Enrich[["fdr_down"]]
     }


     if (length(setdiff(enrich.databases, names(enrich_list)))>0) {
       warning(setdiff(enrich.databases, names(enrich_list)), "not present in your Enrichment analysis.")
       enrich.databases <- intersect(enrich.databases, names(enrich_list))
     }

     if (length(enrich.databases) == 0) {
       stop("Selected Enrich databases not present in your enrichment analysis.")
     }

     for (dat in enrich.databases) {
     Table <- enrich_list[[dat]]
     enrichR.table <- rbind(enrichR.table, Table)
     }

     row.names(enrichR.table) <- enrichR.table$Term

     fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))

     p <- enrichR.table[enrichR.table$Adjusted.P.value < p_adj,"Term"]

     if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
     #pathways.dataframe <- data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value)
     pathways.dataframe <- data.frame(
       Pathway = p,
       gene.ratio = vapply(p, fx, FUN.VALUE = numeric(1)),
       p.value = enrichR.table[p, ]$P.value,
       p.value.adj = enrichR.table[p, ]$Adjusted.P.value
     )
     # Formatting the dataframe for the plot
     pathways.dataframe <- pathways.dataframe[order(pathways.dataframe$p.value.adj),]
     pathways.dataframe$Pathway.num <- as.factor(dim(pathways.dataframe)[1]:1) # as.factor is somehow necessary to write the pathway names
     }

     ###############################################
     # N Top significant pathways
     ###############################################

     top_sig <- head(pathways.dataframe[seq_len(num_term), ], num_term)


     ###############################################
     # Significant pathway barplot
     ###############################################

     top_sig$Log10Adj.P.value <- -log10(top_sig$p.value.adj)

     #remove the description in brackets (GO:0071346)
     top_sig$Pathway <- stringr::str_replace(top_sig$Pathway, "\\s*\\([^\\)]+\\)", " ")

     top_sig <- top_sig[order(top_sig$gene.ratio),]

     p_top_sig <- ggplot(top_sig, ggplot2::aes(x=Pathway, y=gene.ratio, fill=Log10Adj.P.value)) +
       geom_bar(stat="identity")+
       scale_fill_gradient(low = "blue", high = "red",  name = "-Log10(p.adj.value)") +
       scale_x_discrete(limits = top_sig$Pathway)+
       labs(title = cond) +
       coord_flip()

     if (plot) {
       pdf("Top_sig_Pathways_barplot.pdf", 10, 5)
       p_top_sig
       dev.off()
     }

     Pathway <- gene.ratio <- Log10Adj.P.value <- NULL
     return(p_top_sig)

}



