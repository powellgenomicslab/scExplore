#' @title Plot feature on low-dimensional space
#' @description This function allows to plot the values of a feature (gene or metadata column) on
#' a low-dimensional space
#' @param object A Seurat object
#' @param feature Feature to be plotted (gene expression or metadata)
#' @param dims Dimensions to plot. Only two are valid.
#' @param reduction Dimensionality reduction (e.g. \code{pca}, \code{umap}, \code{tsne}, ...)
#' @param group Group variable in metadata
#' @param type If gene expression is used, specify, rata type: \code{counts} or \code{data} slots
#' @param qclip Quantile value to clip gene expression values or continuous variables. This parameter reduces the effect of
#' outlier cells with high gene expression according to a quantile of the distribution (default 0.99). All
#' cells for each gene expressing a value greater to the quantile value in the distribution are rescaled to
#' the corresponding value of that quantile
#' @param label Label clusters when a categorical variable is used. Removes legend.
#' @param alpha Alpha level to adjust transparency of colors
#' @param size Size of points
#' @return A plot with feature of interest
#' @importFrom cowplot plot_grid
#' @importFrom ggrepel geom_label_repel
#' @importFrom cluster clara
#' @export
#' @author José Alquicira Hernández
#' @examples
#'
#' # Plot gene expression of CD8A and CD3D from log-normalized data over a UMAP plot
#'
#' plotFeature(pbmc, "seurat_clusters")
#'


plotFeature <- function(object, feature, dims = c(1,2), reduction = "umap", group = NULL, type = "data", qclip = 0.99, label = TRUE, alpha = 0.7, size = 1){

  if(!is(object, "Seurat")){
    stop("Input object must be of 'Seurat' class")
  }

  # Test reduction existence
  if(!reduction %in% names(slot(object, "reductions"))){
    stop("reduction data not found in object")
  }

  # Extract metadata
  metadata <- slot(object, "meta.data")

  # Test existence oF feature in metadata
  if(feature %in% colnames(metadata)){
    var <- metadata[,feature]
    featureType <- "Values"
  }else{
    expData <- GetAssayData(object, type)
    inExpData <- rownames(expData) == feature
    featureType <- "Expr"
    # Test existence of feature in gene expression data
    if(!any(inExpData)) stop("Feature not present in meta.data or expression data")
    var <- expData[inExpData,]
  }

  # Get groupping variable if provided

  if(!is.null(group)){

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
  cellEmbeddings <- cellEmbeddings[,dims]

  # Validate dimensions
  dimsInEmbeddings <- dims %in% seq_len(ncol(cellEmbeddings))

  if(!all(dimsInEmbeddings)){
    missingDim <- dims[which(!dimsInEmbeddings)]
    stop(paste("Dimension", missingDim, "does not exist in", reduction, "\n  "))
  }



  if(is(var, "numeric") | is(var, "integer")){

    # Clip outliers
    clipValues <- function(x, feature, qclip){
      clip <- quantile(x, qclip)
      if(clip){
        i <- x > clip
        x[i] <- clip
      }else{
        message(paste(qclip, "quantile for feature", feature , "is equal to 0. Clipping step will be ommited"))
      }
      x
    }

    var <- clipValues(var, feature = feature, qclip = qclip)

    # Combine feature with cell embeddings

    cellEmbeddings <- cbind(cellEmbeddings, var)
    names(cellEmbeddings)[3] <- feature

    # Plot data
    dimNames <- names(cellEmbeddings)[1:2]

    if(isGroup){
      cellEmbeddings[,group] <- groupVar
    }

    # Sort embeddings by color
    i <- order(cellEmbeddings[, feature])
    cellEmbeddings <- cellEmbeddings[i,]


    p <-
      ggplot(cellEmbeddings) +
      aes_string(dimNames[1], dimNames[2], color = feature) +
      geom_point(alpha = alpha, size = size) +
      xlab(gsub("_", " ", dimNames[1])) +
      ylab(gsub("_", " ", dimNames[2])) +
      ggtitle(feature) +
      theme_classic() +
      labs(color = featureType) +
      scale_color_viridis_c()


  }else{

    if(is(var, "character")) var <- as.factor(var)

    cellEmbeddings <- cbind(cellEmbeddings, var)
    names(cellEmbeddings)[3] <- feature
    dimNames <- names(cellEmbeddings)[1:2]

    if(is(var, "factor")){
      # Combine embeddings with variable

      nLevs <- length(levels(var))

      if(nLevs <= 16){
        pal <- pal16
      }else if(nLevs <= 21){
        pal <- pal21
      }else if(nLevs > 21){
        pal <- rainbow(n = nLevs)
      }

      if(label){
        cellEmbeddingsByClass <- split(cellEmbeddings, cellEmbeddings[, feature])

        centroids <- as.data.frame(Reduce(rbind,
                                          lapply(cellEmbeddingsByClass,
                                                 function(x){
                                                   x <- apply(x[,c(1,2)], 2, clara, k = 1)
                                                   as.data.frame(lapply(x, "[[", "medoids"))
                                                 }
                                          )))
        centroids[,feature] <- as.factor(names(cellEmbeddingsByClass))
      }

      if(isGroup){
        cellEmbeddings[,group] <- groupVar
      }

      p <-
        ggplot(cellEmbeddings) +
        aes_string(dimNames[1], dimNames[2], color = feature) +
        geom_point(alpha = alpha, size = size) +
        scale_color_manual(values = pal) +
        ggtitle(feature) +
        xlab(gsub("_", " ", dimNames[1])) +
        ylab(gsub("_", " ", dimNames[2])) +
        theme_classic() +
        labs(color = "")


      if(label){
        p <- p + geom_label_repel(aes_string(label = feature),
                                  color = "black",
                                  label.size = NA,
                                  fill = NA,
                                  data = centroids) + theme(legend.position = "none")
      }
    }

  }

  if(!is.null(group)){
    p <- p + facet_wrap(as.formula(paste("~", group)))
  }

  p


}
