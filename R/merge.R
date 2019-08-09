#' @title Plot gene expression of genes and their overlap
#' @description This function allows to plot the expression of genes in individual color scales
#' and their corresponding "merge" value
#' @param object A Seurat object
#' @param gene1 Gene 1
#' @param gene2 Gene 2
#' @param type Data type: \code{counts} or \code{data} slots
#' @param reduction Dimensionality reduction (e.g. \code{pca}, \code{umap}, \code{tsne}, ...)
#' @param qclip Quantile value to clip gene expression values. This parameter reduces the effect of
#' outlier cells with high gene expression according to a quantile of the distribution (default 0.99). All
#' cells for each gene expressing a value greater to the quantile value in the distribution are rescaled to
#' the corresponding value of that quantile
#' @param alpha Alpha level to adjust transparency of colors
#' @param size Size of points
#' @return A grid plot with individual genes and merge gene expression values
#' @importFrom cowplot plot_grid
#' @export
#' @author José Alquicira Hernández
#' @examples
#'
#' # Plot gene expression of CD8A and CD3D from log-normalized data over a UMAP plot
#'
#' plotMerge(pbmc, "CD8A", "CD3D", reduction = "umap", type = "data")
#'


plotMerge <- function(object, gene1, gene2, type = "data", reduction = "umap", qclip = 0.95, alpha = 0.3){

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
  if(!gene2 %in% rownames(expData)){
    stop(paste0("gene '", gene2, "' is not present in data"))
  }

  # Extract embeddings
  cellEmbeddings <- as.data.frame(Embeddings(slot(object, "reductions")[[reduction]]))

  # Extract gene expression values
  expData <- expData[rownames(expData) %in% c(gene1, gene2),] -> expDataGenes


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
  expData <- mapply(clipValues, x = expData, gene = colnames(expData), MoreArgs = list(qclip))


  # Scale data to [0,1] range
  expData <- apply(expData, 2, function(x) x / max(x))

  overlay <- rgb(green = 0, red=expData[,gene1], blue = expData[,gene2])


  # Combine gene/color gradients with dimensionality reductiond data
  cellEmbeddings <- cbind(cellEmbeddings, overlay)
  cellEmbeddings$overlay <- as.character(cellEmbeddings$overlay)

  #cellEmbeddings <- cbind(cellEmbeddings, t(as.matrix(expDataGenes)))
  cellEmbeddings <- cbind(cellEmbeddings, expData)


  # Plot data
  dimNames <- names(cellEmbeddings)[1:2]

  # Sort embeddings by color
  i <- order(cellEmbeddings$overlay)
  cellEmbeddings <- cellEmbeddings[i,]

  ggplot(cellEmbeddings) +
    aes_string(dimNames[1], dimNames[2]) +
    geom_point(color = cellEmbeddings$overlay, alpha = alpha, size = 0.9) +
    xlab(gsub("_", " ", dimNames[1])) +
    ylab(gsub("_", " ", dimNames[2])) +
    ggtitle("Merge") +
    theme_classic() -> p


  pGenes <- mapply(FUN = .plotGeneExp,
                   c(gene1, gene2), c("red", "blue"),
                   embed = list(cellEmbeddings),
                   alpha = list(alpha),
                   SIMPLIFY = FALSE)

  pGenes$merge <- p

  plot_grid(plotlist = pGenes, nrow = 1, rel_widths = c(1, 1, 3/4))


}


.plotGeneExp <- function(gene, col, embed, alpha){

  # Sort embeddings by gene expression values
  i <- order(embed[,gene])
  embed <- embed[i,]

  dimNames <- names(embed)[1:2]

  ggplot(embed) +
    aes_string(dimNames[1], dimNames[2], color = gene) +
    geom_point(alpha = alpha, size = 0.9) +
    scale_color_gradient(low = "black", high = col) +
    ggtitle(gene) +
    xlab(gsub("_", " ", dimNames[1])) +
    ylab(gsub("_", " ", dimNames[2])) +
    theme_classic() +
    labs(color = "Expr")

}
