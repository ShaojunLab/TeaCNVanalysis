#source("~/Library/Mobile Documents/com~apple~CloudDocs/code/R/ArchRHeatmap.R")
library(ComplexHeatmap)
#' Plot Peak-Gene links Heatmap from input matrix
#' This function plots side by side heatmaps of linked ATAC and Gene regions from signaficant links
#' @param mATAC matrix of ATAC peaks
#' @param mRNA matrix of RNA peaks, have the same rows and columns with mATAC
#' @param k An integer describing the number of k-means clusters to group peak-to-gene links prior to plotting heatmaps.
plotLinksHeatmap_pair <- function(
    mATAC = NULL, 
    mRNA = NULL,
    meta_links = NULL,
    colAno = NULL,
    k = 25,
    nPlot = 25000,
    limitsATAC = c(-2, 2),
    limitsRNA = c(-2, 2),
    groupBy = "Clusters",
    palGroup = NULL,
    palATAC = paletteContinuous("blueYellow"),
    palRNA = paletteContinuous("solarExtra")){
  nd <- list("ComplexHeatmap","dplyr","circlize","ggplot2","ggsci","S4Vectors")
  lapply(nd, require, character.only = TRUE)
  if(!is.null(meta_links)){
    if(nrow(meta_links) == 0){stop("No significant peak2genelinks found in meta_links!")}
  }
  # Normalization Matrices by row
  if(!is.null(mATAC)){
    mATAC <- .rowZscores(mATAC)
    tree_row=.cluster_mat(mATAC,distance = "correlation", method = "complete")	 #clust by row (features)
    clust=cutree(tree_row$clust, k)
    if(nrow(mATAC) > nPlot){
      nPK <- nPlot * table(clust) / length(clust) 
      splitK <- split(seq_len(nrow(mATAC)), clust)
      kDF <- lapply(seq_along(splitK), function(x){
        idx <- sample(splitK[[x]], floor(nPK[x]))
        k <- rep(x, length(idx))
        DataFrame(k = k, idx = idx)
      }) %>% Reduce("rbind", .)
    }else{
      kDF <- DataFrame(k = clust, idx = seq_len(nrow(mATAC)))
    }
   # match(rownames(kDF),colnames(mATAC))
    #order row
    rowOrder <- tree_row$clust$order
    colOrder <- colnames(mATAC)
    # Plot Heatmaps
    htATAC <- .ArchRHeatmap(
      mat = mATAC[rowOrder,colOrder],
      scale = FALSE,
      limits = limitsATAC,
      color = palATAC, 
      colData =  colAno,
      colorMap = NULL,
      clusterCols = FALSE,
      clusterRows = F,
      split = NULL,
      labelRows = FALSE,
      labelCols = FALSE,
      draw = FALSE,
      name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks")
    )
  }
 # "Constructing RNA Heatmap!"
  if(!is.null(mRNA)){
    mRNA <- .rowZscores(mRNA)
    if(!exists("rowOrder")){rowOrder <-rownames(mRNA)}
    if(!exists("colOrder")){colOrder <-colnames(mRNA)}
    htRNA <- .ArchRHeatmap(
      mat = mRNA[rowOrder,colOrder], 
      scale = FALSE,
      limits = limitsRNA,
      color = palRNA, 
      colData = colAno,
      colorMap = NULL,
      clusterCols = FALSE,
      clusterRows = F,
      split = NULL,
      labelRows = FALSE,
      labelCols = FALSE,
      draw = FALSE,
      name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
    )
  }
  if((!is.null(mRNA))&(!is.null(mATAC))){
    return(htATAC+htRNA)
  }else if(is.null(mRNA)&(!is.null(mATAC))){
    return(htATAC)
    }else{return(htRNA)}
 
}

#' @description normalization data to [-2,2]
.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  row_Sds <- matrixStats::rowSds(m)
  row_Sds <- ifelse(row_Sds==0,1e-6,row_Sds)
  z <- sweep(m - rowMeans(m), 1, row_Sds,`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

#' @description group on columns
.groupMeans <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}


#' Re-map a character vector of labels from an old set of labels to a new set of labels
#'
#' This function takes a character vector of labels and uses a set of old and new labels
#' to re-map from the old label set to the new label set.
#'
#' @param labels A character vector containing lables to map.
#' @param newLabels A character vector (same length as oldLabels) to map labels to from oldLabels.
#' @param oldLabels A character vector (same length as newLabels) to map labels from to newLabels
#' @export
mapLabels <- function(labels = NULL, newLabels = NULL, oldLabels = names(newLabels)){
  
  if(length(newLabels) != length(oldLabels)){
    stop("newLabels and oldLabels must be equal length!")
  }
  
  if(!requireNamespace("plyr", quietly = TRUE)){
    labels <- paste0(labels)
    oldLabels <- paste0(oldLabels)
    newLabels <- paste0(newLabels)
    labelsNew <- labels
    for(i in seq_along(oldLabels)){
      labelsNew[labels == oldLabels[i]] <- newLabels[i]
    }
    paste0(labelsNew)
  }else{
    paste0(plyr::mapvalues(x = labels, from = oldLabels, to = newLabels))
  }
  
}

#####
.cluster_mat = function(mat,distance="euclidean", method="ward.D2"){

  if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
    stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
    stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
  }
  if(distance[1] == "correlation"){
    d = as.dist(1 - cor(t(mat)))
  }else{
    if(class(distance) == "dist"){
      d = distance
    }else{
      d = dist(mat, method = distance)
    }
  }
  
  return(list(clust=hclust(d, method = method),dist=d))
}

  