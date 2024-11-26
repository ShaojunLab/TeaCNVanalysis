#source("fun_pseudoBulk_mat.R",chdir = T)
#source("colorPalette.R",chdir = T)
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(scales)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(dplyr)
  library(Seurat)
  library(Signac)
  library(ggthemes)
})
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
  library(patchwork)
}

cb_pattern <- c("#4dbbd5ff","#E64b35ff","#00a087ff","#a65628",'#fdbf6f',"#3c5488ff","#f39b7fff","#91d1c2ff","#8491b4ff","#7e6148ff","#b09c85ff","#7876b1ff","#377eb8","#4daf4a","#984ea3","#ff7f00","#f781bf","#999999",'#1f78b4','#b2df8a','#33a02c', '#fb9a99')
cols_type <- c("#B80C09","#0B4F6C","#D4B483","#4357AD","#48A9A6","#7699D4","#E4DFDA","#9448BC","#B9929F","#A6E1FA","#F0C808","#FF7F11")

#plot for count matrix 
#' @title plot4Count()
#' @description Data visualization for the # of reads(or peaks) per cell, and # of reads (or cells) per peak (column : cell)
plot4Count <- function(matrix,outdir="./",outfile="rawdataCheck.pdf",
                       addline=FALSE,
                       format = "pdf",
                       min_Nreads_per_cell=NULL,
                       min_Nfeature_per_cell=NULL,
                       min_Nreads_per_feature=NULL,
                       min_Ncell_per_feature=NULL){
  if(!is.matrix(matrix)){matrix <- as.matrix(matrix)}
  suppressMessages(library(pheatmap))
  suppressMessages(library(RColorBrewer))
  plotDir <- outdir
  if(!file.exists(plotDir)){dir.create(plotDir,recursive=T)}
  reads_per_cell = colSums(matrix)
  reads_per_peak = rowSums(matrix)
  peaks_per_cell = colSums(matrix>0)
  cells_per_peak = rowSums(matrix>0)
  if(format=="pdf"){
    pdf(file=paste0(plotDir,"/",outfile),width = 8,height = 8)
  }else{
    png(file=paste0(plotDir,"/",outfile),width = 8,height =8,units = 'in',res= 300)
  }
  par(mfrow=c(3,2))
  #layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  hist(log10(reads_per_cell),main='reads per cell',col='wheat',breaks = 20,prob = TRUE)
  lines(density(log10(reads_per_cell)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(reads_per_cell)), col = "blue",lty=2)
  text(x=log10(mean(reads_per_cell))+0.2,y=max(density(log10(reads_per_cell))$y-0.1),paste0("Mean: ",round(mean(reads_per_cell))),cex = .75,col="blue")
  
  if(addline){
    abline(v = log10(min_Nreads_per_cell), col = "red",lty=2)
  }
  hist(log10(peaks_per_cell), main='peaks per cell', col='wheat',breaks = 20,prob = TRUE)
  lines(density(log10(peaks_per_cell)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(peaks_per_cell)), col = "blue",lty=2)
  text(x=log10(mean(peaks_per_cell))+0.2,y=max(density(log10(peaks_per_cell))$y-0.1),paste0("Mean: ",round(mean(peaks_per_cell))),cex = .75,col="blue")
  if(addline){abline(v = log10(min_Nfeature_per_cell), col = "red",lty=2)}
  
  hist(log10(reads_per_peak),main='reads per peak',col='wheat',breaks=20,prob = TRUE)
  lines(density(log10(reads_per_peak)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(reads_per_peak)), col = "blue",lty=2)
  text(x=log10(mean(reads_per_peak)) + 0.2, y=max(density(log10(reads_per_peak))$y*0.6),paste0("Mean: ",round(mean(reads_per_peak))),cex = .75,col="blue")
  
  if(addline){abline(v = log10(min_Nreads_per_feature), col = "red",lty=2)}
  hist(log10(cells_per_peak),main='cells per peak',col='wheat',breaks=20,prob = TRUE)
  lines(density(log10(cells_per_peak)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(cells_per_peak)), col = "blue",lty=2)
  text(x=log10(mean(cells_per_peak)) + 0.2, y=max(density(log10(cells_per_peak))$y*0.6),paste0("Mean: ",round(mean(cells_per_peak))),cex = .75,col="blue")
  if(addline){abline(v = log10(min_Ncell_per_feature), col = "red",lty=2)}
  plot(reads_per_cell, peaks_per_cell, log='xy', col='wheat')
  
  dev.off()
}

#' Plot a Heatmap of Identified ATAC Marker Features
#' This function will plot a heatmap of a matrix of Feature-Cell
#' @param @param scaleRows A boolean value that indicates whether the heatmap should display row-wise z-scores instead of raw values.
#' @param filename a pdf file path and name for output.
plotPeakHeatmap <- function(object,features = NULL,assay=NULL,slot = "counts",scaleRows = TRUE,limits = c(-2, 2),pal = NULL,
                            group.by=NULL,group.colors = NULL,filename="heatmap_peaks.pdf",Wid=12,Hei=12,
                            showRowDendrogram = FALSE,clusterRows = FALSE,
                            heatmap_legend_side="bottom",annotation_legend_side = "bottom"){
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  
  cells <- colnames(x = object)
  data <- as.data.frame(x = as.matrix(x = GetAssayData(
    object = object,
    slot = slot)[features,, drop = FALSE]))
  group.by <- group.by %||% 'seurat_clusters'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  groups.use <- groups.use[order(groups.use[,group.by]),,drop=F]
  colOrder <- rownames(groups.use)
  
  #filter low expressed features
  data <- data[rowSums(data>0)>1,,drop=F]
  data <- as.matrix(data)
  if(nrow(data) == 0){
    stop("No Markers Found!")
  }
  message(sprintf("Identified %s markers!", nrow(data)))
  if(scaleRows){
    mat <- .rowZscores(data,limit=F)
  }else{
    mat <- data
  }
  
  
  mat <-mat[,colOrder]
  
  pal_col = paletteContinuous(set ="blueYellow", n = 100)
  #Prepare Limits if needed
  breaks <- NULL
  if (scaleRows & !is.null(limits)) {
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }else{
    limits <- c(round(min(mat),2), round(quantile(mat,0.98),2))
  }
  
  breaks <- seq(min(limits), max(limits), length.out = length(pal_col))
  pal_col <- circlize::colorRamp2(breaks, pal_col)
  
  default.colors <- c(hue_pal()(length(x = levels(x = groups.use[,group.by]))))
  if (!is.null(x = names(x = group.colors))) {
    cols <- unname(obj = group.colors[levels(x = groups.use[,group.by])])
  } else {
    cols <- group.colors[1:length(x = levels(x = groups.use[,group.by]))] %||% default.colors
  }
  
  colorMap=list()
  names(cols) = levels(groups.use[,group.by])
  colorMap[[group.by]] <- cols
  
  ht1Anno <- HeatmapAnnotation(
    df = groups.use,
    col = colorMap, 
    show_annotation_name = TRUE,
    gp = gpar(col = "NA"),
    annotation_legend_param =
      list(
        nrow = min(5, ceiling(nrow(unique(groups.use))/2))
      )
  )
  
  if(scaleRows){
    fileName <- paste0(assay," Z-Scores\n", length(unique(rownames(mat))), " features, ",ncol(mat)," cells")
    
  }else{
    fileName <- paste0(assay," counts\n", length(unique(rownames(mat))), " features, ",ncol(mat)," cells")
    
  }
  ht_tumor <- Heatmap(mat,
                      name = fileName,
                      show_column_names = F, 
                      show_row_names = T,
                      border = F,
                      row_title = NULL, #"%s",
                      row_title_gp = gpar(fontsize = 10),
                      heatmap_legend_param = list(
                        direction = "horizontal",
                        #title_position = "leftcenter-rot",
                        legend_height = unit(3, "cm")),
                      col=pal_col,
                      cluster_columns = F,
                      cluster_rows = clusterRows, 
                      show_row_dend = showRowDendrogram, 
                      clustering_method_rows = "ward.D2",
                      #Annotation
                      top_annotation = ht1Anno)
  
  pdf(filename,width = Wid,height = Hei)
  draw(ht_tumor, merge_legend = TRUE,heatmap_legend_side=heatmap_legend_side,annotation_legend_side = annotation_legend_side)
  dev.off()
  
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

scale_rows <- function(x, min = -2, max = 2, limit = FALSE) {
  if(!limit){min=0;max =1}
  apply(x, 1, function(row) {
    min_row <- min(row)
    max_row <- max(row)
    min + (max - min) * (row - min_row) / (max_row - min_row)
  })
}


#' Plot Peak-Gene links Heatmap
#' @param gplinks a 'GRanges' links from Links(object[["peaks"]])
#' @param score_cutoff coef.result>score_cutoff, default 0.05
#' @param cell_ano_names column names of object for annotation
#' @param cell_ano_colors a list of colors with names of cell_ano_names
plotLinksHeatmap <- function(object,assay.rna="SCT",assay.atac="peaks",
                             assay.cnv=NULL,
                             cells=NULL,Nsample=500,Nrows_show=10000,group_by=NULL,group.colors=NULL,
                             cell_ano_names=NULL,cell_ano_colors=NULL,
                             gplinks = NULL, pvalue_cutoff = 0.05,score_cutoff = 0.05,scale=TRUE,
                             limitsATAC = c(-2, 2),
                             limitsRNA = c(-2, 2),
                             palATAC = paletteContinuous("blueYellow"),
                             palRNA = paletteContinuous("solarExtra"),
                             filename="./heatmap_links.pdf",Wid=16,Hei=10
){
  nd <- list("ComplexHeatmap","dplyr","circlize","ggplot2","ggsci","fastcluster","Matrix","foreach","doParallel")
  lapply(nd, require, character.only = TRUE)
  res<- list()
  if(is.null(gplinks)){
    if(length(Links(obj[[assay.atac]]))>0){
      gplinks <- Links(obj[[assay.atac]])
      gplinks <- as.data.frame(gplinks)
    }else{stop(paste0("No links found in assays ",assay.atac,", run LinkPeaks() first."))}
  }
  gplinks_sig <- gplinks[gplinks$pvalue<pvalue_cutoff&gplinks$score>score_cutoff,,drop=F]
  if(nrow(gplinks_sig)<1){
    stop(paste0("No significant gene-peak links were found in links."))
  }
  gplinks_sig <- gplinks_sig[order(gplinks_sig$pvalue,gplinks_sig$score,decreasing = c(FALSE,TRUE)),]
  res$gplinks_sig <- gplinks_sig
  
  #down-sample links for plot
  if(Nrows_show<nrow(gplinks_sig)){
    sigLink <- gplinks_sig[1:Nrows_show,,drop=F]
  }
  
  if(length(object@neighbors)<1){
    object <- FindNeighbors(object, dims = 1:30,force.recalc=T,return.neighbor=TRUE,k.param=100,assay = assay.rna)
  }
  if(is.null(cells)){
    cells=colnames(object)
  }
  
  if(Nsample>= length(cells)){
    randomseeds=cells
  }else{
    randomseeds=sample(cells,size=Nsample)
  }
  res$randomseeds <- randomseeds
  pseudoExp <- gen_pseudo_mat(object,randomseeds,assay=assay.rna,expMatrix=obj2[[assay.rna]]@data)
  pseudoPeak <- gen_pseudo_mat(object,randomseeds,assay=assay.rna,expMatrix=obj2[[assay.atac]]@data)

  
  if(!is.null(assay.cnv)){
    pseudoCNV <- gen_pseudo_mat(object,randomseeds,assay=assay.rna,expMatrix=obj2[[assay.cnv]]@counts,rowValue="singlecellCNV")
    
    ##
    gene.filt <- intersect(rownames(pseudoCNV),sigLink$gene)
    sigLink <- sigLink[sigLink$gene %in% gene.filt,,drop=F]  ###Filtering genes!
    cnv_exp=pseudoCNV[gene.filt,,drop=F]
  }
  
  rna_exp=pseudoExp[unique(sigLink$gene),,drop=F]
  peak_exp=pseudoPeak[unique(sigLink$peak),,drop=F]
  #cell order
  show_pdata=unique(c(group_by,cell_ano_names))
  cell_ano <- object[[]][randomseeds,show_pdata]
  cell_ano <- cell_ano[order(cell_ano[,group_by],cell_ano[,cell_ano_names]),,drop=F]
  cell_odr <- rownames(cell_ano)
  #ordering for peak matrix
  peak_exp <- as.matrix(peak_exp)
  #colnames(peak_exp) <- gsub("\\.","-",colnames(peak_exp))
  peaks.keep  <- which(rownames(peak_exp)%in%sigLink$peak)
  if(length(peaks.keep)>0){
    peak_exp <- peak_exp[peaks.keep,]
    index1=na.omit(match(sigLink$peak,row.names(peak_exp)))
  }else{stop("No significant peak found in peak matrix. Check if peaks are identical.")}
  
  cell_indx1 <- match(cell_odr,colnames(peak_exp))
  peak_exp=as.matrix(peak_exp[index1,cell_indx1])
  
  #ordering gene expression matrix
  rna_exp <- as.matrix(rna_exp)
  #colnames(rna_exp) <- gsub("\\.","-",colnames(rna_exp))
  genes.keep  <- which(rownames(rna_exp)%in%sigLink$gene)
  #genesName.keep <- intersect(rownames(rna_exp)[genes.keep],rownames(cnv_exp))
  #genes.keep <-  which(rownames(rna_exp)%in%genesName.keep)
  if(length(genes.keep)>0){
    rna_exp <- rna_exp[genes.keep,]
    index2=na.omit(match(sigLink$gene,row.names(rna_exp)))
  }else{
    stop("No linked genes found in expression matrix. Check if gene names are identical.")
  }
  cell_indx2 <- match(cell_odr,colnames(rna_exp))
  rna_exp=as.matrix(rna_exp[index2,cell_indx2])
  
  if(!is.null(assay.cnv)){
    #ordering CNV matrix
    cnv_exp <- as.matrix(cnv_exp)  
    index3=match(sigLink$gene,rownames(cnv_exp))
    cell_indx3 <- match(cell_odr,colnames(cnv_exp))
    cnv_exp=as.matrix(cnv_exp[index3,cell_indx3])
    cnv_exp[is.na(cnv_exp)]<- 1   ####TBD
  }
  
  if(scale){
    mATAC <- .rowZscores(peak_exp,min = min(limitsATAC), max = max(limitsATAC),limit=T)
    mRNA <- .rowZscores(rna_exp,min = min(limitsRNA), max = max(limitsRNA),limit=T)
    if(!is.null(assay.cnv)){
      mCNV <- .rowZscores(cnv_exp,min = min(limitsRNA), max = max(limitsRNA),limit=T) 
    }
  }else{
    mATAC <-peak_exp
    # mATAC[mATAC>3] <- 3
    mRNA <- rna_exp
    # mRNA[mRNA>3] <- 3
    if(!is.null(assay.cnv)){
      mCNV <- cnv_exp
    }
  }
  
  sparse_mRNA <- as(mRNA, "sparseMatrix")
  # Compute hierarchical clustering on the distance matrix
  dist_matrix <- dist(sparse_mRNA)
  hc <- fastcluster::hclust(as.dist(dist_matrix), method = "ward.D2")
  rowOrder <- hc$order
  res$clust <- hc
  
  if(is.null(cell_ano_colors)){
    col_map_func <- function(dat){
      res = list()
      for(n in colnames(dat)){
        tmp_name = sort(unique(dat[, n]))
        tmp_col = pal_jco("default",alpha = 0.8)(length(tmp_name))
        names(tmp_col) = tmp_name
        res[[n]] = tmp_col
      }
      return(res)
    }
    cell_ano_colors <-  col_map_func(cell_ano)
  }
  
  default.colors <- c(hue_pal()(length(x = levels(x = cell_ano[,group_by]))))
  if (!is.null(x = names(x = group.colors))) {
    cols <- unname(obj = group.colors[levels(x = cell_ano[,group_by])])
  } else {
    cols <- group.colors[1:length(x = levels(x = cell_ano[,group_by]))] %||% default.colors
  }
  names(cols) = levels(cell_ano[,group_by])
  cell_ano_colors[[group_by]] <- cols
  
  
  ht1Anno <- HeatmapAnnotation(
    df = cell_ano,
    col = cell_ano_colors, 
    show_annotation_name = TRUE,
    gp = gpar(col = "NA"),
    annotation_legend_param =
      list(
        nrow = min(5, ceiling(nrow(cell_ano)/2))
      )
  )
  
  # Plot Heatmaps
  htATAC <- Heatmap(mATAC[rowOrder,cell_odr],
                   name = paste0("ATAC Z-Scores\n", nrow(gplinks_sig), " links, ",length(unique(gplinks_sig$peak))," peaks"),
                   show_column_names = F, 
                   show_row_names = F,
                   border = F,
                   row_title = NULL, #"%s",
                   row_title_gp = gpar(fontsize = 10),
                   heatmap_legend_param = list(
                     direction = "horizontal",
                     #title_position = "leftcenter-rot",
                     legend_height = unit(3, "cm")),
                   col=palATAC,
                   cluster_columns = F,
                   cluster_rows = F,
                   #Annotation
                   top_annotation = ht1Anno)
    
    

  htRNA <- Heatmap(mRNA[rowOrder,cell_odr],
                  name = paste0("RNA Z-Scores\n",nrow(gplinks_sig)," links, ", length(unique(gplinks_sig$gene)), " genes"),
                  show_column_names = F, 
                  show_row_names = F,
                  border = F,
                  row_title = NULL, #"%s",
                  row_title_gp = gpar(fontsize = 10),
                  heatmap_legend_param = list(
                    direction = "horizontal",
                    #title_position = "leftcenter-rot",
                    legend_height = unit(3, "cm")),
                  col=palRNA,
                  cluster_columns = F,
                  cluster_rows = F,
                  #Annotation
                  top_annotation = ht1Anno) 
  htList <- htATAC+htRNA
  if(!is.null(assay.cnv)){
    htCNV <- Heatmap(mCNV[rowOrder,cell_odr],
                     name = paste0("CNV Z-Scores\n",nrow(gplinks_sig[gplinks_sig$gene %in% rownames(pseudoCNV),,drop=F])," links, ", length(unique(intersect(rownames(pseudoCNV),gplinks_sig$gene))), " genes"),
                     show_column_names = F, 
                     show_row_names = F,
                     border = F,
                     row_title = NULL, #"%s",
                     row_title_gp = gpar(fontsize = 10),
                     heatmap_legend_param = list(
                       direction = "horizontal",
                       #title_position = "leftcenter-rot",
                       legend_height = unit(3, "cm")),
                     col=colorRampPalette(colors = c("darkblue", "white", "darkred"))(16),
                     cluster_columns = F,
                     cluster_rows = F,
                     #Annotation
                     top_annotation = ht1Anno) 
    htList <- htATAC+htRNA + htCNV
  }
  pdf(filename,width = Wid,height = Hei)
  draw(htList, merge_legend = TRUE,heatmap_legend_side="bottom",annotation_legend_side = "bottom")
  dev.off()
  return(res)
}


#' A ggplot-based dot plot wrapper function
#'
#' This function is edited.
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param color A numeric/categorical vector used to determine the coloration for each point.
#' @param discrete A boolean value indicating whether the supplied data is discrete (`TRUE`) or continuous (`FALSE`).
#' @param discreteSet The name of a custom palette from `colorPalettes` to use for categorical/discrete color.
#' This argument is only used if `discrete` is set to `TRUE`.
#' @param continuousSet The name of a custom palette from `colorPalettes` to use for numeric color.
#' This argument is only used if `discrete` is set to `FALSE`.
#' @param labelMeans A boolean value indicating whether the mean of each categorical/discrete color should be labeled.
#' @param pal A custom palette used to override discreteSet/continuousSet for coloring vector.
#' @param defaultColor The default color for points that do not have another color applied (i.e. `NA` values).
#' @param highlightPoints A integer vector describing which points to hightlight. The remainder of points will be colored light gray.
#' @param colorDensity A boolean value indicating whether the density of points on the plot should be indicated by color.
#' If `TRUE`, continuousSet is used as the color palette.
#' @param size The numeric size of the points to be plotted.
#' @param xlim A numeric vector of two values indicating the lower and upper bounds of the x-axis on the plot.
#' @param ylim A numeric vector of two values indicating the lower and upper bounds of the y-axis on the plot.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum and minimum
#' values if `xlim` and `ylim` are not provided. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param randomize A boolean value indicating whether to randomize the order of the points when plotting.
#' @param seed A numeric seed number for use in randomization.
#' @param colorTitle A title to be added to the legend if `color` is supplied.
#' @param colorOrder A vector that allows you to control the order of palette colors associated with the values in `color`.
#' For example if you have `color` as `c("a","b","c")` and want to have the first color selected from the palette be used for
#' "c", the second color for "b", and the third color for "a", you would supply the `colorOrder` as `c("c", "b", "a")`.
#' @param colorLimits A numeric vector of two values indicating the lower and upper bounds of colors if numeric. Values
#' beyond these limits are thresholded.
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param legendSize The size in inches to use for plotting the color legend.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param labelAsFactors A boolean indicating whether to label the `color` input as a numeric factor (`TRUE`) or with a character string (`FALSE`).
#' @param fgColor The foreground color of the plot.
#' @param bgColor The background color of the plot.
#' @param bgWidth The background relative width size of the halos in the labeling.
#' @param labelSize The numeric font size of labels.
#' @param addFit A string indicating the method to use for adding a fit/regression line to the plot (see `ggplot2::geom_smooth()` methods).
#' If set to `NULL`, no fit/regression line is added.
#' @param rastr A boolean value that indicates whether the plot should be rasterized using `ggrastr`. This does not rasterize
#' lines and labels, just the internal portions of the plot.
#' @param dpi The resolution in dots per inch to use for the plot.
#' @export
ggPoint <- function(
    x = NULL, 
    y = NULL, 
    color = NULL, 
    discrete = TRUE, 
    discreteSet = "stallion",
    continuousSet = "solarExtra", 
    labelMeans = TRUE,  
    pal = NULL, 
    defaultColor = "lightGrey",
    highlightPoints = NULL,
    colorDensity = FALSE,
    size = 1, 
    xlim = NULL, 
    ylim = NULL, 
    extend = 0.05, 
    xlabel = "x", 
    ylabel = "y", 
    title = "", 
    subtitle = "",
    randomize = FALSE, 
    seed = 1,
    colorTitle = NULL, 
    colorOrder = NULL, 
    colorLimits = NULL,
    alpha = 1, 
    baseSize = 10, 
    legendSize = 3,
    ratioYX = 1, 
    labelAsFactors = TRUE,
    fgColor = "black", 
    bgColor = "white", 
    bgWidth = 1,
    labelSize = 3,
    addFit = NULL, 
    rastr = FALSE, 
    dpi = 300,
    ...
){
  
  stopifnot(length(y) == length(x))
  if(length(x) < 5){
    stop("x must be at least length 5 to plot!")
  }
  
  if(randomize){
    set.seed(seed)
    idx <- sample(seq_along(x), length(x))
  }else{
    idx <- seq_along(x)
  }
  
  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    x <- x[include]
    y <- y[include]
    if(!is.null(color)){
      color <- color[include]
    }
  }
  
  if(is.null(xlim)){
    xlim <- range(df$x) %>% extendrange(f = extend)
  }
  
  if(is.null(ylim)){
    ylim <- range(df$y) %>% extendrange(f = extend)
  }
  
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  
  #Plot
  if (is.null(color) & !colorDensity) {
    
    p <- ggplot(df[idx,], aes(x = x, y = y)) + 
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
      xlab(xlabel) + ylab(ylabel) + 
      labs(title = title,
           subtitle = subtitle) +
      theme_4_plot(baseSize = baseSize)
    
    if(rastr){
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
    }else{
      p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
    }
    
  }else {
    
    if(colorDensity){
      
      discrete <- FALSE
      df <- .getDensity(x, y, n = 100, sample = NULL) #change
      df <- df[order(df$density), ,drop=FALSE]
      df$color <- df$density
      
      if(is.null(colorTitle)){
        colorTitle <- "density"
      }
      
    }else if(discrete){
      
      if(!is.null(highlightPoints)){
        if(length(highlightPoints) < length(color)){
          color[-highlightPoints] <- "Non.Highlighted"
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      color <- paste0(color)
      
      if(!is.null(colorOrder)){
        if(!all(color %in% colorOrder)){
          stop("Not all colors are in colorOrder!")
        }
      }else{
        colorOrder <- gtools::mixedsort(unique(color))
      }
      
      if(is.null(colorTitle)){
        colorTitle <- "color"
      }
      
      stopifnot(length(color) == nrow(df))
      df$color <- factor(color, levels = colorOrder)
      
      if(labelAsFactors){
        df$color <- factor(
          x = paste0(paste0(match(paste0(df$color), paste0(levels(df$color)))), "-", paste0(df$color)), 
          levels = paste0(seq_along(levels(df$color)), "-", levels(df$color))
        )
        if(!is.null(pal)){
          names(pal) <- paste0(levels(df$color))[match(names(pal), colorOrder)]
        }
        colorOrder <- paste0(levels(df$color))
      }
      
    }else{
      stopifnot(length(color) == nrow(df))
      if(!is.null(highlightPoints)){
        if(length(highlightPoints) < length(color)){
          color[-highlightPoints] <- NA
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      if(!is.null(colorLimits)){
        color[color < min(colorLimits)] <- min(colorLimits)
        color[color > max(colorLimits)] <- max(colorLimits)
      }
      df$color <- color
    }
    
    p <- ggplot(df[idx,], aes(x = x, y = y, color = color)) +  
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) + 
      xlab(xlabel) + ylab(ylabel) + 
      labs(title = title,
             subtitle = subtitle) + 
      theme_4_plot(baseSize = baseSize) +
      theme(legend.direction = "horizontal", legend.box.background = element_rect(color = NA)) +
      labs(color = colorTitle)
    
    if(rastr){
      
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha, 
        raster.width = min(par('fin')), 
        raster.height = (ratioYX * min(par('fin')))
      )
      
    }else{
      
      p <- p + geom_point(size = size, alpha = alpha)
      
    }
    
    if (discrete) {
      
      if (!is.null(pal)) {
        p <- p + scale_color_manual(values = pal)
      }else {
        pal <- paletteDiscrete(set = discreteSet, values = colorOrder)
        if(!is.null(highlightPoints)){
          pal[grep("Non.Highlighted", names(pal))] <- "lightgrey"
        }
        p <- p + scale_color_manual(values = pal) +
          guides(color = guide_legend(override.aes = list(size = legendSize, shape = 15)))
      }
      
      if (labelMeans) {
        
        dfMean <- split(df, df$color) %>% lapply(., function(x) {
          data.frame(x = median(x[, 1]), y = median(x[, 2]), color = x[1, 3])
        }) %>% Reduce("rbind", .)
        
        if(labelAsFactors){
          dfMean$label <- stringr::str_split(paste0(seq_len(nrow(dfMean))), pattern = "\\-", simplify=TRUE)[,1]
        }else{
          dfMean$label <- dfMean$color
        }
        dfMean$text <- stringr::str_split(dfMean$color, pattern = "-", simplify = TRUE)[,1]
        
        
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        xo <- bgWidth * diff(range(df$x)) / 300
        yo <- bgWidth * diff(range(df$y)) / 300
        for (i in theta) {
          p <- p + 
            geom_text(data = dfMean, 
                      aes_q(
                        x = bquote(x + .(cos(i) * xo)),
                        y = bquote(y + .(sin(i) * yo)),
                        label = ~text
                      ),
                      size = labelSize,
                      color = bgColor
            )
        }
        
        if(is.null(fgColor)){
          p <- p + geom_text(data = dfMean, aes(x = x, y = y, color = color, label = label), size = labelSize, show.legend = FALSE)
        }else{
          p <- p + geom_text(data = dfMean, aes(x = x, y = y, label = label), color = fgColor, size = labelSize, show.legend = FALSE) 
        }
        
      }
      
    }else{
      
      if (!is.null(pal)) {
        if(!is.null(colorLimits)){
          p <- p + scale_colour_gradientn(colors = pal, limits=colorLimits, na.value = "lightgrey")
        }else{
          p <- p + scale_colour_gradientn(colors = pal, na.value = "lightgrey")
        }
      }else {
        if(!is.null(colorLimits)){
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), limits=colorLimits, na.value = "lightgrey")
        }else{
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), na.value = "lightgrey")
        }
      }
    }
    
  }
  
  if (!is.null(addFit)) {
    result.spe <- cor.test(df$x, df$y,method="spearman")
    rho.spe <- round(result.spe$estimate,3)
    pvalue.spe <- ifelse(result.spe$p.value==0,0,format(result.spe$p.value, scientific = TRUE, digits = 3))
    result.pea <- cor.test(df$x, df$y,method="pearson")
    rho.pea <- round(result.pea$estimate,3)
    pvalue.pea <- ifelse(result.pea$p.value==0,0,format(result.pea$p.value, scientific = TRUE, digits = 3))

    
    p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, color = "grey") + 
      ggtitle(paste0(title, "\nPearson = ", rho.spe, ", p-value = ",pvalue.spe,
                     "\nSpearman = ", rho.pea, ", p-value = ",pvalue.pea))
  }
  
  p <- p + theme(legend.position = "bottom", legend.key = element_rect(size = 2))#, legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(0.1, 'cm'))
  
  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }
  
  return(p)
  
}
.getDensity <- function(x = NULL, y = NULL, n = 100, sample = NULL, densityMax = 0.95){
  df <- data.frame(x=x,y=y)
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  df$density <- dens$z[ii]
  df$density[df$density > quantile(unique(df$density),densityMax)] <- quantile(unique(df$density),densityMax) 
  if(!is.null(sample)){
    df <- df[sample(nrow(df), min(sample,nrow(df))),]
  }
  return(df)
}

theme_4_plot <- function(
    color = "black",
    textFamily = "sans",
    baseSize = 10, 
    baseLineSize = 0.5,
    baseRectSize = 0.5,
    plotMarginCm = 1,
    legendPosition = "bottom",
    legendTextSize = 5,
    axisTickCm = 0.1,
    xText90 = FALSE,
    yText90 = FALSE
){
  
  theme <- theme_bw() + theme(
    text = element_text(family = textFamily),
    axis.text = element_text(color = color, size = baseSize), 
    axis.title = element_text(color = color, size = baseSize),
    title = element_text(color = color, size = baseSize),
    plot.margin = unit(c(plotMarginCm, plotMarginCm, plotMarginCm, plotMarginCm), "cm"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = color, size = (4/3) * baseRectSize * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    axis.ticks.length = unit(axisTickCm, "cm"), 
    axis.ticks = element_line(color = color, size = baseLineSize * (4/3) * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(color = color, size = legendTextSize),
    legend.box.background = element_rect(color = NA),
    #legend.box.background = element_rect(fill = "transparent"),
    legend.position = legendPosition,
    strip.text = element_text(size = baseSize, color="black")#,
    #plot.background = element_rect(fill = "transparent", color = NA)
  )
  
  if(xText90){
    theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  if(yText90){
    theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90, vjust = 1))
  }
  
  return(theme)
  
}
paletteContinuous <- function(
    set = "solarExtra", 
    n = 256, 
    reverse = FALSE
){
  
  pal <- colorPalettes[[set]]
  palOut <- colorRampPalette(pal)(n)
  
  if(reverse){
    palOut <- rev(palOut)
  }
  
  return(palOut)
  
}

##Plot
#' @param object Seurat object
DensityScatter_v2 <- function(object,x = 'log10_nFrags', y = 'TSS.enrichment', log_x = FALSE, quantiles = TRUE
                           ){
  df <- object@meta.data[,c(x,y),drop=F]
  if(quantiles){
    quan_x <- as.numeric(round(quantile(df[,x],c(0.05,0.1,0.9,0.95)),2))
    quan_y <- as.numeric(round(quantile(df[,y],c(0.05,0.1,0.9,0.95)),2))
    subtitle = sprintf("%s: %d%%: %.2f  %d%%: %.2f  %d%%: %.2f  %d%%: %.2f\n%s: %d%%: %.2f  %d%%: %.2f  %d%%: %.2f  %d%%: %.2f",x, 5,quan_x[1],10,quan_x[2],90,quan_x[3],95,quan_x[4],y,5,quan_y[1],10,quan_y[2],90,quan_y[3],95,quan_y[4])
  }
  
  if(log_x){
    df[,x] <- log10(df[,x])
    xlabel <- paste0("Log ",x)
  }else{xlabel = x}

  ggPoint(
    x = df[,x], 
    y = df[,y], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = xlabel,
    ylabel = y,
    subtitle = subtitle,
    xlim = quantile(df[,x], probs = c(0.01,0.99)),
    ylim = c(0, quantile(df[,y], probs = 0.99))
  ) + geom_hline(yintercept = 3, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
}

# Function to change x labels of individual violin plots
VlnPlot_v2 <- function(object,group_by = NULL,
                       features.plot=c('nCount_ATAC', 'nFeature_ATAC',"log10_nFrags",'log10_nCount_ATAC',
                                       'TSS.enrichment', 'nucleosome_signal'), 
                       ident.include=NULL,
                       pt.size = 0,
                       xlab="",color_pattern=NULL,ncol=NULL) {
  
  # Main function
  main_function <- function(object = object,group_by=group_by, features.plot = features.plot,pt.size=pt.size, 
                            ident.include = ident.include, xlab = xlab,color_pattern=color_pattern,add_median_point=add_median_point) {
    pvln <- VlnPlot(object = object, features = features.plot, group.by = group_by, pt.size=pt.size,
            idents = ident.include,cols =color_pattern) + 
      theme(legend.position = 'none')+ 
      labs(x = xlab)
    
    return(pvln)
  }
  
  # Apply main function on all features
  p <- lapply(X = features.plot, object = object,group_by=group_by, pt.size=pt.size,ident.include = ident.include, xlab = xlab,color_pattern=color_pattern,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  if(is.null(ncol)){ncol=ceiling(sqrt(length(features.plot)))}
  cowplot::plot_grid(plotlist = p, ncol = ncol)
}


VlnPlotSelf <- function(obj,
                        features = NULL,lim = NULL,round_p = 0 ,
                        ncol = 1,color_pattern=NULL){
  require(Seurat);
  require(cowplot);
  require(ggplot2);
  
  # calculate median
  give.m = function(x) {
    la.df = data.frame(y = median(x) + max(x)/8,label = paste0(round(mean(x), round_p)));
    return(la.df)
  }
  if(is.null(color_pattern)){color_pattern <- cb_pattern}
  P = lapply(features,FUN = function(x){
    
    if(is.null(lim)){
      a = min(obj@meta.data[[x]])
      b = max(obj@meta.data[[x]])
      lim = c(a,b)
    }
    
    Seurat::VlnPlot(obj,features =x,cols = color_pattern,pt.size = 0,raster=FALSE) +
      stat_summary(fun = mean, geom = "point", col = "black") +
      stat_summary(fun.data = give.m,size = 4,position = 'identity',
                   geom = "text",
                   col = "black")+NoLegend()+ylab('') + ylim(lim[1],lim[2])
  })
  plot_grid_args <- c(P[c(1:length(features))],ncol = ncol)
  do.call(cowplot::plot_grid, plot_grid_args)
}

DimPlot.multi <- function(object,group_by=NULL,reduction = "umap",
                          color_list=NULL,label=FALSE,legend_position = "right"){
  Pl <- list()
  for(g in group_by){
    n_groups <- length(unique(object[[]][,g])) 
    nCol <- ceiling(n_groups/16)
    if(is.null(color_list)){
      color_by <- NULL
    }else{
      color_by <- color_list[[g]]
    }
    Pl[[g]] <- DimPlot(object, reduction = reduction, group.by = g, label = label, label.size = 4, repel = TRUE,cols =color_by,raster=FALSE) +
      ggtitle(g)+
      theme(legend.position = legend_position) &
      guides(color = guide_legend(override.aes = list(size=3),ncol = nCol))
  }
  cowplot::plot_grid(plotlist = Pl, ncol = ceiling(sqrt(length(group_by)))) 
}

Celltype.marker.featureplot <- function(obj,MainCellType_markers_plot=NULL,cellName,reduction = NULL,ncol = 4){
  if(is.null(MainCellType_markers_plot)){
    MainCellType_markers_plot <- list(immune_cell = 'PTPRC',
                                      T_cell = c('CD3D','CD3E','CD8A','CD8B','CD4'),
                                      B_cell = c('CD19','CD14','MS4A1','CD79A'),
                                      Monocyte_cell = c('CD68'),
                                      Epithelial_cell = c('EPCAM','KRT5','KRT18'),
                                      Endothelial_cell= c('PECAM1','VWF','PLVAP'),
                                      Fibroblast_cell =c("ACTA2","COL1A1",'COL1A2','DCN'), 
                                      Mucosa_cell = c('IGFBP2','TFF1'),
                                      Plasma = c("SDC1","CD38","MZB1","IGHA1"),
                                      Myeloid = c("CD68","CD14","APOE","LYZ","CD163"),
                                      Mast = c("ITGAX","KIT","TPSAB1"),
                                      DC = c("CLEC9A","CD1C","THBD"),
                                      pDC = c("ITGAX","LILRA4","UGCG","CLEC4C"),
                                      Nutrophils = c("CXCL8","S100A8","S100A9"))
  }
  p <- FeaturePlot(obj,features = c(MainCellType_markers_plot[[cellName]]),
                   reduction = reduction,pt.size = 0.1,
                   max.cutoff = 'q95',ncol = ncol) 
  return(p)   
}


#' @author Yanru Zhang
cellType.Plot <- function(obj,i = 'seurat_clusters', violinMarker = NULL,dotMarker = NULL,col_palette=NULL,assay.rna="RNA"){
  
  # function to obtain a 
  give.mat = function(marker,obj){
    suppressMessages(library(tidyr))
    suppressMessages(library(reshape2))
    Gene = marker$Gene
    Gene = intersect(Gene,rownames(obj))
    dat = as.matrix(obj[[assay.rna]]@data[Gene,])
    dat = as.data.frame(t(dat))
    dat[,'Clusters'] = obj@meta.data[,i]
    dat[,'CellId'] = colnames(obj)
    dat_long <- melt(dat,id.vars =c('Clusters','CellId'),variable.name = 'Gene' ,value.name='Expr')
    mat = merge(dat_long,marker,by.x = 'Gene',by.y = 'Gene')
    return(mat)}
  
  violinMarker$CellType = factor(violinMarker$CellType,levels = unique(violinMarker$CellType))
  mat = give.mat(obj = obj,marker = violinMarker)
  
  if(is.null(col_palette)){
    col_palette <- cb_pattern
  }

  if(length(col_palette)<length(unique(mat[,'Clusters']))){
    col_palette <- paletteDiscrete( values = length(unique(mat[,'Clusters'])))
  }
  P1 = ggplot(mat,aes(x = Clusters,y = Expr,fill = Clusters)) + 
    geom_violin(scale = "width")+
    theme_few()+
    scale_fill_manual(values = col_palette)+
    facet_grid( CellType+Gene ~. ,scales = 'free_y')+
    theme(strip.background=element_rect(colour="black", fill="white"),strip.text = element_text(size = 10))+
    NoLegend()
  
  Gene = unique(as.character(dotMarker$Gene)[nrow(dotMarker):1])
  Gene <- Gene[Gene%in%rownames(obj)]
  P2= DotPlot(object = obj,features = Gene) +
    coord_flip()+
    scale_color_gradientn(colours = c("#053061","#2166ac","#4393c3","#d1e5f0","#fddbc7",
                                      "#d6604d","#b2182b","#67001f"))+
    theme(axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = 0.5))
  P1+P2 
}

Marker.Gene.Heatmap <- function(obj,markerGene,top = 5){
  top5MG = markerGene %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
  new_obj = ScaleData(obj, features = top5MG$gene)
  DoHeatmap(new_obj, features = top5MG$gene,size =5) + scale_fill_gradientn(colors = c("#045a8d","white","#a50f15")) 
}
