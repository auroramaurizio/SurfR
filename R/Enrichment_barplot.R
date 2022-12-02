#' Enrichment_barplot
#'
#' Barplot representing the top up-regulated or down-regulated significant pathways
#'
#' @param dfList A list of enrichment tables for up and down-regulated genes in the different enrichrdatabases. Output of Enrichment.R function.
#' @param enrich.databases Vector of EnrichR databases to consider 
#' @param p_adj Double. Minimum Adjusted pvalue threshold for the enrichment
#' @param num_term Double. Number of up-regulated and dw-regulated terms to represent
#' @param cond String. Title of the plot.
#' @examples
#' Enrichment_barplot(dfList, enrich.databases = c("GO_Biological_Process_2018","GO_Cellular_Component_2018"), p_adj = 0.01, num_term = 5, cond = "UP")
#' @section Warning:
#' Bla bla bla
#' @family aggregate functions
#' @seealso \code{\link{hello}} for counts data and metadata download, and \code{\link{hello}} for Gene2SProtein analysis
#' @export
#' 


Enrichment_barplot <- function(dfList, enrich.databases  = c("GO_Biological_Process_2021",
                                                     "GO_Cellular_Component_2021",
                                                     "GO_Molecular_Function_2021"), p_adj = 0.05, num_term = 10, cond = "UP") {
     
  
     import::here(stringr)
     import::here(openxlsx)
  
  
     enrichR.table = data.frame()

     for (dat in enrich.databases) {
     Table <- openxlsx::read.xlsx(xlsxFile=dfList, sheet=dat, startRow=1, colNames=T, rowNames=T,
                               detectDates=F, skipEmptyRows=T, skipEmptyCols=T, na.strings="NA", fillMergedCells=F)

     enrichR.table = rbind(enrichR.table, Table)
     }
     
     fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))

     p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < p_adj,])

     if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
     pathways.dataframe = data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value)
     # Formatting the dataframe for the plot
     pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
     pathways.dataframe$Pathway.num = as.factor(dim(pathways.dataframe)[1]:1) # as.factor is somehow necessary to write the pathway names
     }
     
     ###############################################
     # N Top significant pathways 
     ###############################################

     top_sig = head(pathways.dataframe[1:num_term,], num_term)
     
     ###############################################
     # Significant pathway barplot
     ###############################################
     
     top_sig$Log10Adj.P.value = -log10(top_sig$p.value.adj)
     
     #remove the description in brackets (GO:0071346)
     top_sig$Pathway = str_replace(top_sig$Pathway, "\\s*\\([^\\)]+\\)", " ")
     
     top_sig <- top_sig[order(top_sig$gene.ratio),] 
     
     
     p_top_sig <-ggplot(top_sig, aes(x=Pathway, y=gene.ratio, fill=Log10Adj.P.value)) +
       geom_bar(stat="identity")+
       scale_fill_gradient(low = "blue", high = "red",  name = "-Log10(p.adj.value)") + 
       scale_x_discrete(limits = top_sig$Pathway)+
       labs(title = cond) +
       coord_flip()
     
     pdf("Top_sig_Pathways_barplot.pdf", 10, 5)
     p_top_sig
     dev.off() 
     
}
     
    

