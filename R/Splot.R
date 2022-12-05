#' Splot function
#'
#' Plot a barplot with features of Surface Protein
#'
#' @param SurfaceProteins_df Output dataframe of Gene2SProtein function.
#' @param group.by Name of columns to plot. Default = Membranome.Almen.main-class.
#' @param cols.use Vector of colors, each color corresponds to an identity class.
#' By default, ggplot assigns colors.
#' @param main Plot title. Default = Almen main class.
#' @return plot objec created by ggplot2,
#' which can be assigned and further customized.
#' @examples
#'  GeneNames = c("CIITA", "EPCAM", "DLK1", "CD24", "CDCP1", "LYVE1", "ABCD1", "VAMP1")
#'  SurfaceProteins_df = Gene2SProtein(GeneNames, input_type = "gene_name")
#'  Splot(SurfaceProteins_df)
#' @export
Splot <- function(SurfaceProteins_df,
                  group.by = "Membranome.Almen.main-class",
                  cols.use = NULL,
                  main = "Almen main class") {
  import::here(ggplot2)
  import::here(scales)

  SurfaceProteins_df$v2plot = SurfaceProteins_df[,group.by]


  if (is.null(x = cols.use)) {
    cols.use = hue_pal()(length(x = levels(as.factor(SurfaceProteins_df$v2plot))) )
    if (length(which (is.na(SurfaceProteins_df$v2plot)) > 0)) {
      warning("NA value in your dataframe")
      # Add grey for NA protein
      cols.use = c(cols.use, "#CCCCCC")
    }
  } else if (length(cols.use) < length(x = unique(SurfaceProteins_df$v2plot))) {
    stop(paste("you have",length(x = unique(SurfaceProteins_df$v2plot)), "unique elements and supplied only",length(cols.use),"color \n",
               "Be carefull to NA value"))
  }

  plot <- ggplot(SurfaceProteins_df, aes(x=v2plot)) +
    geom_bar(color = 'black', fill = cols.use )+
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=22, vjust = 1, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position="top",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = group.by, y = "N", title = main)

  ggplot2 <- ggrepel <- geom_label_repel <- hue_pal <- aes <- color <- shape <- sampleNames <- NULL
  return(plot)
}

