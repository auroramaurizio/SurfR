#' plotPCA function
#'
#' Plot PCA highlighting one or two data features
#'
#' @param matrix Filtered count matrix in CPM or RPKM with gene on the row and sample ID on the column.
#' @param metadata Sample metadata, row.names must be samples names.
#' @param nTOP number of top genes to use for principal components,
#' selected by highest row variance
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param centering Logical. If \code{TRUE} center PCs
#' @param scaling Logical. If \code{TRUE} scales PCs
#' @param color.by Name of one or more metadata columns to color point by.
#' @param shape.by Name of one or more metadata columns to shape point by. If NULL, all points are circles \(default\).
#' @param pt.size Size of the points in the plot.
#' @param cols.use Vector of colors, each color corresponds to an identity class.
#' By default, ggplot assigns colors.
#' @param shape.use Vector of shape, each shape corresponds to an identity class.
#' @param main Plot title. Default = PCA.
#' @param label Logical. If \code{TRUE} adds samples label. Default = FALSE.
#' @param new.label If NULL, use the sample names as in metadata row.names.
#' Otherwise you can specify new labels.
#' @return PCA plot objec created by ggplot2,
#' which can be assigned and further customized.
#' @examples
#' # Simulation of bulk RNA data
#' countData <- matrix(floor(runif(10000, min=0, max=101)),ncol=4)
#' colnames(countData) <- paste0("sample", 1:4)
#' rownames(countData) <- paste0("gene", 1:(10000/4))
#' metadata = data.frame(samplesID = paste0("sample", 1:4),
#'                      condition = factor(c("A","A","B","B")),
#'                      therapy = factor(c("T1","T2","T1","T2")))
#' row.names(metadata) <- metadata$samplesID
#' library(edgeR)
#' plotPCA(matrix = cpm(countData),
#'         metadata = metadata,
#'         nTOP = 100,
#'         dims = c(1,2),
#'         color.by = "condition", shape.by = "therapy",
#'         label = FALSE, main = "PCA")
#' @family plot functions
#' @importFrom ggplot2 ggplot aes geom_point unit xlab ylab coord_fixed ggtitle theme theme_bw geom_hline geom_vline element_text element_rect scale_color_manual scale_shape_manual
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats prcomp
#' @importFrom scales hue_pal
#' @export
#'
#'


plotPCA <- function(matrix,
                    metadata,
                    nTOP = 500,
                    dims = c(1,2),
                    centering = TRUE,
                    scaling = TRUE,
                    color.by = NULL,
                    shape.by = NULL,
                    pt.size = 6,
                    cols.use = NULL,
                    shape.use = NULL,
                    main = "PCA",
                    label = FALSE,
                    new.label = NULL) {


  #import::here(ggplot2)
  #import::here(ggrepel)
  #import::here(stats, prcomp)

  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }

  matrix <- matrix[,match(row.names(metadata), colnames(matrix))]
  if (length(x = row.names(metadata)) != length(x = colnames(matrix))) {
    if (!all(row.names(metadata) == colnames(matrix))) {
      stop("matrix column and metadata rows don't match")
    }
  }


  # sort by variance and select topN
  vary <- apply(matrix,1,var)
  vary_s <- sort(vary, decreasing = TRUE)
  TOP_N <- names(vary_s[1:nTOP])
  mtx_TOP <- matrix[TOP_N,]


  # pca with prcomp
  pcx = dims[1]; pcy = dims[2]
  centering = centering
  scaling = scaling

  # PCA
  pca = prcomp(t(mtx_TOP), center=centering, scale=scaling)
  var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
  score = as.data.frame(pca$x)
  score

  # plot paramters
  xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
  ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
  cum = var[pcx]+var[pcy]
  names = rownames(pca$x)

  if (is.null(x = color.by)) {
    color.by = colnames(metadata)[1]
    warning("color.by option is NULL, colorig by first column of metadata")
  }

  score$color = as.factor(metadata[, color.by])

  if (is.null(x = shape.by)) {
    score$shape = ""
    shape.by = "shape"
    warning("shape.by option is NULL. All points will have same shape")
  } else {
    score$shape = as.factor(metadata[, shape.by])
  }

  if (is.null(x = cols.use)) {
    #cols.use = scales::hue_pal()(length(x = levels(x = score$color)))
    cols.use = hue_pal()(length(x = levels(x = score$color)))
  } else if (length(cols.use) < length(x = levels(x = score$color))) {
    stop(paste("you have",length(x = levels(x = score$color)), "factors and supplied only",length(cols.use),"color"))
  }

  if (is.null(x = shape.use)) {
    shape.use = c(16:25, 1:15)
  } else if (length(shape.use) < length(x = levels(x = score$shape))) {
    stop(paste("you have",length(x = levels(x = score$shape)), "factors and supplied only",length(shape.use),"shape"))
  }


  if (label) {
    if (is.null(x = new.label)) {
      score$sampleNames = row.names(metadata)
    } else {
      score$sampleNames = new.label
    }
    # pca = ggplot2::ggplot(score, ggplot2::aes(score[,pcx], y=score[,pcy], color=color, shape = shape, label = sampleNames)) +
    #   ggrepel::geom_label_repel(data= score, ggplot2::aes(x=score[,pcx], y=score[,pcy],
    #                                     color=color, label = sampleNames),
    #                    size = 6,  box.padding = ggplot2::unit(1, "lines"),
    #                    point.padding = ggplot2::unit(0.1, "lines"),
    #                    segment.color = 'grey50') +
    #   ggplot2::geom_point(size= pt.size)+
    #   ggplot2::xlab(xlab) +
    #   ggplot2::ylab(ylab) +
    #   ggplot2::coord_fixed() + ggplot2::ggtitle(main) +
    #   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    #   ggplot2::theme_bw(base_size = 16) + #+ theme(legend.position = "right") +
    #   ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    #   ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    #   ggplot2::theme(plot.title = ggplot2::element_text(color="black", size=16, face="bold.italic"),
    #         axis.text.x = ggplot2::element_text(angle = 0, face = "bold", color = "black", size=12, hjust =.5),
    #         axis.title.x = ggplot2::element_text(face = "bold", color = "black", size = 14),
    #         axis.text.y = ggplot2::element_text(angle = 0, face = "bold", color = "black", size=12),
    #         axis.title.y = ggplot2::element_text(face = "bold", color = "black", size = 14),
    #         legend.text = ggplot2::element_text(face = "bold", color = "black", size = 10),
    #         legend.position="right",
    #         panel.background = ggplot2::element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    #   ggplot2::scale_color_manual(values=cols.use, name = color.by) +
    #   ggplot2::scale_shape_manual(values = shape.use, name = shape.by)
    # } else {
    #pca = ggplot2::ggplot(score, ggplot2::aes(score[,pcx], y=score[,pcy], color=color, shape = shape)) +
    #   ggplot2::geom_point(size= pt.size)+
    #  ggplot2::xlab(xlab) +
    #   ggplot2::ylab(ylab) +
    #   ggplot2::coord_fixed() + ggplot2::ggtitle(main) +
    #   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    #   ggplot2::theme_bw(base_size = 16) + #+ theme(legend.position = "right") +
    #   ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    #   ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    #  ggplot2::theme(plot.title = ggplot2::element_text(color="black", size=16, face="bold.italic"),
    #         axis.text.x = ggplot2::element_text(angle = 0, face = "bold", color = "black", size=12, hjust =.5),
    #         axis.title.x = ggplot2::element_text(face = "bold", color = "black", size = 14),
    #         axis.text.y = ggplot2::element_text(angle = 0, face = "bold", color = "black", size=12),
    #         axis.title.y = ggplot2::element_text(face = "bold", color = "black", size = 14),
    #         legend.text = ggplot2::element_text(face = "bold", color = "black", size = 10),
    #         legend.position="right",
    #         panel.background = ggplot2::element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    #   ggplot2::scale_color_manual(values=cols.use, name = color.by) +
    #   ggplot2::scale_shape_manual(values = shape.use, name = shape.by)

    # }

    pca = ggplot(score, aes(score[,pcx], y=score[,pcy], color=color, shape = shape, label = sampleNames)) +
      geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy],
                                        color=color, label = sampleNames),
                       size = 6,  box.padding = ggplot2::unit(1, "lines"),
                       point.padding = ggplot2::unit(0.1, "lines"),
                       segment.color = 'grey50') +
      geom_point(size= pt.size)+
      xlab(xlab) +
      ylab(ylab) +
      coord_fixed() + ggtitle(main) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw(base_size = 16) + #+ theme(legend.position = "right") +
      geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
      geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
      theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
            axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =.5),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            legend.text = element_text(face = "bold", color = "black", size = 10),
            legend.position="right",
            panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
      scale_color_manual(values=cols.use, name = color.by) +
      scale_shape_manual(values = shape.use, name = shape.by)
  } else {
    pca = ggplot(score, aes(score[,pcx], y=score[,pcy], color=color, shape = shape)) +
      geom_point(size= pt.size)+
      xlab(xlab) +
      ylab(ylab) +
      coord_fixed() + ggtitle(main) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw(base_size = 16) + #+ theme(legend.position = "right") +
      geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
      geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
      theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
            axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=12, hjust =.5),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=12),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            legend.text = element_text(face = "bold", color = "black", size = 10),
            legend.position="right",
            panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
      scale_color_manual(values=cols.use, name = color.by) +
      scale_shape_manual(values = shape.use, name = shape.by)

  }



  ggplot2 <- ggrepel <- geom_label_repel <- stats <- prcomp <- hue_pal <- aes <- color <- shape <- sampleNames <- NULL
  return(pca)
}

