#' @title Plot expression of genes
#' @description This function allows to plot the expression of genes in individual color scales
#' and their corresponding "merge" visualization
#' @param object A Seurat object
#' @param gene1 Gene 1
#' @param gene2 Gene 2 (optional)
#' @param dims Dimensions to plot. Only two are valid.
#' @param reduction Dimensionality reduction (e.g. \code{pca}, \code{umap}, \code{tsne}, ...)
#' @param group Group variable in metadata
#' @param type Data type: \code{counts} or \code{data} slots
#' @param qclip Quantile value to clip gene expression values. This parameter reduces the effect of
#' outlier cells with high gene expression according to a quantile of the distribution (default 1: No clipping). All
#' cells for each gene expressing a value greater to the quantile value in the distribution are rescaled to
#' the corresponding value of that quantile
#' @param alpha Alpha level to adjust transparency of colors
#' @param size Size of points
#' @param bgColor Background color
#' @param returnGrid If \code{TRUE}, a grid is generated using cowplot. Otherwise, a list with the plots is
#' returned
#' @return A grid plot with individual genes and merge gene expression values
#' @importFrom cowplot plot_grid
#' @importFrom cowplot get_legend
#' @export
#' @author José Alquicira Hernández
#' @examples
#'
#' # Plot gene expression of CD8A and CD3D from log-normalized data over a UMAP plot
#'
#' plotMerge(pbmc, "CD8A", "CD3D", reduction = "umap", type = "data")
#'


plotExp <- function(object, gene1, gene2 = NULL, dims = c(1,2), reduction = "umap",
                    group = NULL, subset = NULL, type = "data", qclip = 1, alpha = 0.7, size = 0.3,
                    bgColor = "#171716", returnGrid = TRUE){

  if(!is(object, "Seurat")){
    stop("Input object must be of 'Seurat' class")
  }

  # Test reduction existence
  if(!reduction %in% names(slot(object, "reductions"))){
    stop("reduction data not found in object")
  }

  # Extract gene expression matrix
  expData <- GetAssayData(object, type)

  # Test existence of genes to plot
  if(!gene1 %in% rownames(expData)){
    stop(paste0("gene '", gene1, "' is not present in data"))
  }

  if(!is.null(gene2)){
    if(!gene2 %in% rownames(expData)){
      stop(paste0("gene '", gene2, "' is not present in data"))
    }
  }


  # Get groupping variable if provided

  if(!is.null(group)){

    metadata <- slot(object, "meta.data")

    if(group %in% colnames(metadata)){
      groupVar <- metadata[,group]
      isGroup <- TRUE
    }else{
      stop("Group variable not present in meta.data")
    }


  }else{
    isGroup <- FALSE
  }

  # Extract embeddings
  cellEmbeddings <- as.data.frame(Embeddings(slot(object, "reductions")[[reduction]]))

  if(!is.null(subset)){
    i <- rownames(cellEmbeddings) %in% subset
    if(!length(i)){
      stop("No cells were obtained after subsetting. Check provided cell ids")
    }
    cellEmbeddings <- cellEmbeddings[i,]
    expData <- expData[,i]
    if(isGroup) groupVar <- groupVar[i]
  }


  # Validate dimensions
  dimsInEmbeddings <- dims %in% seq_len(ncol(cellEmbeddings))

  if(!all(dimsInEmbeddings)){
    missingDim <- dims[which(!dimsInEmbeddings)]
    stop(paste("Dimension", missingDim, "does not exist in", reduction, "\n  "))
  }

  cellEmbeddings <- cellEmbeddings[,dims]

  # Extract gene expression values

  if(is.null(gene2)){
    expData <- expData[rownames(expData) %in% gene1, , drop = FALSE]
  }else{
    expData <- expData[rownames(expData) %in% c(gene1, gene2),]
  }

  # Clip outliers
  clipValues <- function(x, gene, qclip){
    clip <- quantile(x, qclip)
    if(clip){
      i <- x > clip
      x[i] <- clip
    }else{
      message(paste(qclip, "quantile for gene", gene , "is equal to 0. Clipping step will be ommited"))
    }
    x
  }

  expData <- as.data.frame(Matrix::t(expData))
  expData <- mapply(clipValues, x = expData, gene = colnames(expData), MoreArgs = list(qclip)) -> expDataGenes


  # Scale data to [0,1] range
  expData <- apply(expData, 2, function(x) x / max(x))

  red <- expData[,gene1]

  if(is.null(gene2)){
    blue <- 0
  }else{
    blue <-  expData[,gene2]
  }

  overlay <- rgb(green = 0, red = red, blue = blue)

  # Combine gene/color gradients with dimensionality reductiond data
  cellEmbeddings <- cbind(cellEmbeddings, overlay)
  cellEmbeddings$overlay <- as.character(cellEmbeddings$overlay)

  #cellEmbeddings <- cbind(cellEmbeddings, t(as.matrix(expDataGenes)))
  cellEmbeddings <- cbind(cellEmbeddings, expDataGenes)


  # Plot data
  dimNames <- names(cellEmbeddings)[1:2]

  # Add group variable if available
  if(isGroup){
    cellEmbeddings[,group] <- groupVar
  }


  # Sort embeddings by color
  i <- order(cellEmbeddings$overlay)
  cellEmbeddings <- cellEmbeddings[i,]


  ggplot(cellEmbeddings) +
    aes_string(dimNames[1], dimNames[2]) +
    geom_point(color = cellEmbeddings$overlay, alpha = alpha, size = size) +
    xlab(gsub("_", " ", dimNames[1])) +
    ylab(gsub("_", " ", dimNames[2])) +
    ggtitle("Merge") +
    theme_classic() +
    theme(panel.background = element_rect(fill = bgColor, colour = NA)) -> p

  if(isGroup){
    p <- p + facet_wrap(as.formula(paste("~", group)))
  }


  genes <- c(gene1, gene2)
  colors <- c("red", "blue")
  colors <- colors[seq_along(genes)]

  pGenes <- mapply(FUN = .plotGeneExp,
                   genes, colors,
                   embed = list(cellEmbeddings),
                   alpha = list(alpha),
                   size = list(size),
                   bgColor = list(bgColor),
                   SIMPLIFY = FALSE)

  if(is.null(gene2)){
    p <- pGenes[[1]]
    if(isGroup){
      p <- p + facet_wrap(as.formula(paste("~", group)))
    }
    return(p)
  }


  if(isGroup){
    legends <- mapply(function(x, lab){ x + labs(color = lab)}, pGenes, names(pGenes), SIMPLIFY = FALSE)
    legends <- lapply(legends, get_legend)
    legend <- plot_grid(plotlist = legends, ncol = 1)
    plot_grid(p, legend, nrow = 1, rel_widths = c(1, 0.15))
  }else{

    pGenes$merge <- p

    if(returnGrid){
      plot_grid(plotlist = pGenes, nrow = 1, rel_widths = c(1, 1, 4/5))
    }else{
      pGenes
    }
  }


}


.plotGeneExp <- function(gene, col, embed, alpha, size, bgColor){

  # Sort embeddings by gene expression values
  i <- order(embed[,gene])
  embed <- embed[i,]

  dimNames <- names(embed)[1:2]

  ggplot(embed) +
    aes_string(dimNames[1], dimNames[2], color = gene) +
    geom_point(alpha = alpha, size = size) +
    scale_color_gradient(low = "black", high = col) +
    ggtitle(gene) +
    xlab(gsub("_", " ", dimNames[1])) +
    ylab(gsub("_", " ", dimNames[2])) +
    theme_classic() +
    labs(color = "Expr") +
    theme(panel.background = element_rect(fill = bgColor, colour = NA))

}
