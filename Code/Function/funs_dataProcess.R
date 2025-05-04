#Functions for data processing
suppressMessages({
  library(dplyr)
  library(Signac)
  library(ggplot2)
})


nCell_statistic <- function(obj,group_by,assay = "RNA", slot = "counts"){
  meta <- obj@meta.data
  df <- meta %>%
    # dplyr::group_by(!!sym(group_by)) %>%
    # dplyr::summarise(nCell=n())
  dplyr::count(!!sym(group_by))
  colnames(df)[2]<- "nCell"
  
  groups <- split(x = colnames(obj), f = obj@meta.data[,group_by])
  group_stats <- c()
  for(group_name in names(groups)) {
    cells <- groups[[group_name]]
    expr_matrix <- GetAssayData(obj, assay = assay, slot = slot)[, cells]
    unique_genes_count <- sum(rowSums(expr_matrix) > 0)
    total_counts <- sum(expr_matrix)
    n_g1 <- c(group_name,unique_genes_count, total_counts)
    group_stats <- rbind(group_stats,n_g1)
  }
  group_stats <- as.data.frame(group_stats)
  group_stats[,2] <- as.numeric(group_stats[,2])
  group_stats[,3] <- as.numeric(group_stats[,3])
  colnames(group_stats) <- c(group_by,"nFeature","nCount")
  df <- left_join(df,group_stats)
  
  return(df)
}

#' @param percent.mt =10 
#' @param nCount_ATAC.min = 500
#' @param nCount_ATAC.max = 15000
#' @param nFeature_ATAC.min = 300
#' @param nFeature_ATAC.max = 5000
#' @param removeDoubletByGene = T
#' @param removeCellByQuantile = F
#' @param max_percent default 0.95 如果设置按照百分比进行过滤，则需设置此参数
#' @param min_percent default 0.05 如果设置按照百分比进行过滤，则需设置此参数
#' @author Yanru zhang
Filter.Low.Quality.Cell.scMulti <-function(obj,
                                           nCount_ATAC.min = NULL,
                                           nCount_ATAC.max = NULL,
                                           nFeature_ATAC.min = NULL,
                                           nFeature_ATAC.max = NULL,
                                           nCount_RNA.min = NULL,
                                           nCount_RNA.max = NULL,
                                           nFeature_RNA.min = NULL,
                                           nFeature_RNA.max = NULL,
                                           percent.mt = NULL,
                                           removeCell_by_quantile = F,
                                           max.percent = 0.95,
                                           min.percent = 0.05,
                                           nucleosome_signal.max = NULL,
                                           TSS.enrichment.min = NULL
){
  if(removeCell_by_quantile){
    nCount_ATAC.min <- quantile(obj$nCount_ATAC,min.percent)
    nCount_ATAC.max <- quantile(obj$nCount_ATAC,max.percent)
    nFeature_ATAC.min <- quantile(obj$nFeature_ATAC,min.percent)
    nFeature_ATAC.max <- quantile(obj$nFeature_ATAC,max.percent)
  }
  
  #建立索引
  t.index = rep(TRUE, length=ncol(obj))
  
  #按照参数获得通过过滤条件后保留的细胞
  if (removeCell_by_quantile) {
    cat('Filter cells based on the provided percentage...\n')
  }
  if(!is.null(nCount_ATAC.min)){
    t.index = t.index & (obj$nCount_ATAC >= nCount_ATAC.min);
    cat('Filter out cells with too low nCount_ATAC,', sum(!(obj$nCount_ATAC >= nCount_ATAC.min)),"cells were filtered out; \n")
  }
  
  if(!is.null(nCount_ATAC.max)){
    t.index = t.index & (obj$nCount_ATAC <= nCount_ATAC.max);
    cat("Filter out cells with too high nCount_ATAC,", sum(!(obj$nCount_ATAC <= nCount_ATAC.max)),"cells were filtered out; \n")
  }
  
  if(!is.null(nFeature_ATAC.min)){
    t.index = t.index & (obj$nFeature_ATAC >= nFeature_ATAC.min);
    cat("Filter out cells with too low nFeature_ATAC,", sum(!(obj$nFeature_ATAC >= nFeature_ATAC.min)),"cells were filtered out; \n")
  }
  
  if(!is.null(nFeature_ATAC.max)){
    t.index = t.index & (obj$nFeature_ATAC <= nFeature_ATAC.max);
    cat("Filter out cells with too high nFeature_ATAC,",sum(!(obj$nFeature_ATAC <= nFeature_ATAC.max)),"cells were filtered out; \n")
  }
  if(!is.null(nCount_RNA.min)){
    t.index = t.index & (obj$nCount_RNA >= nCount_RNA.min);
    cat("Filter out cells with too low nCount_RNA,",sum(!(obj$nCount_RNA >= nCount_RNA.min)),"cells were filtered out; \n")
  }
  if(!is.null(nCount_RNA.max)){
    t.index = t.index & (obj$nCount_RNA <= nCount_RNA.max);
    cat("Filter out cells with too high nCount_RNA,",sum(!(obj$nCount_RNA <= nCount_RNA.max)),"cells were filtered out; \n")
  }
  if(!is.null(nFeature_RNA.min)){
    t.index = t.index & (obj$nFeature_RNA >= nFeature_RNA.min);
    cat("Filter out cells with too low nFeature_RNA,",sum(!(obj$nCount_RNA >= nCount_RNA.min)),"cells were filtered out; \n")
  }
  if(!is.null(nFeature_RNA.max)){
    t.index = t.index & (obj$nFeature_RNA <= nFeature_RNA.max);
    cat("Filter out cells with too high nFeature_RNA,",sum(!(obj$nFeature_RNA <= nFeature_RNA.max)),"cells were filtered out; \n")
  }
  if(!is.null(percent.mt)){
    t.index = t.index & (obj$percent.mt <= percent.mt);
    cat("Filter out cells with too high percent.mt,",sum(!(obj$percent.mt <= percent.mt)),"cells were filtered out; \n")
  }
  if(!is.null(nucleosome_signal.max)){
    t.index = t.index & (obj$nucleosome_signal <= nucleosome_signal.max);
    cat("Filter out cells with too high nucleosome_signal,",sum(!(obj$nucleosome_signal <= nFeature_ATAC.max)),"cells were filtered out; \n")
  }
  if(!is.null(TSS.enrichment.min)){
    t.index = t.index & (obj$TSS.enrichment >= TSS.enrichment.min);
    cat("Filter out cells with too low TSS.enrichment.min,",sum(!(obj$TSS.enrichment >= TSS.enrichment.min)),"cells were filtered out; \n")
  }
  #过滤细胞
  obj = subset(obj, cells = colnames(obj)[t.index])
  
  
  cat("Filtering done. Remained cells: ",ncol(obj),"\n")
  return(obj)
}

QC.statistic <- function(object,meta.cols=c("nCount_ATAC","nFeature_ATAC","log10_nFrags","TSS.enrichment","nucleosome_signal")){
  meta <- object@meta.data
  meta <- meta[,meta.cols,drop=F]
  
  percentile_fun <- function(x) {
    c(quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
    mean = mean(x, na.rm = TRUE))
  }
  stats_df <- meta %>% summarise_all(percentile_fun) %>%
    as.data.frame()
  row.names(stats_df) <- c("0%","25%","50%","75%","100%","Mean")
  stats_df <- t(stats_df)
  stats_df <- round(stats_df,2)
  return(stats_df)
}

runClustering.scATAC <- function(obj,assay = "ATAC",TopFeatures.cutoff='q5',
                                 algorithm = 3,
                                 dims=2:30,
                                 resolution = 0.8,
                                 features = NULL,
                                 features.exclude = NULL){
  DefaultAssay(obj) <- assay
  obj <- FindTopFeatures(obj, min.cutoff = TopFeatures.cutoff)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  if(!is.null(features)){
    features <- VariableFeatures(obj) 
    features <- features[!features %in%features.exclude ]
    obj[[assay]]@var.features <- features
  }
  obj <- RunUMAP(obj, reduction = 'lsi', dims =dims)
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = dims)
  obj <- FindClusters(object = obj,algorithm = algorithm, verbose = FALSE,resolution =resolution)
  return(obj)
}


getGeneActivity <- function(obj,assays="ATAC"){
  DefaultAssay(obj) = assays
  gene.activities <- GeneActivity(obj,assay=assays)
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  obj[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
  
  # obj <- NormalizeData(
  #   object = obj,
  #   assay = 'GeneActivity',
  #   normalization.method = 'LogNormalize',
  #   scale.factor = median(obj$GeneActivity,na.rm=TRUE))
  return(obj)
}

MainCellType_markers_enrich = list(T_cell_marker = c('CD3D', 'CD3E', 'CD8A', 'NKG7', 'MKI67','STMN1',
                                                     'CD4','CD8A',
                                                     'CCR7', 'SELL', 'IL7R', 'LEF1',
                                                     'FOXP3', 'TNFRSF4', 'TNFRSF18', 'IL2RA', 'CTLA4', 'BATF', 'RTKN2', 'IKZF2',
                                                     'GZMK', 'GZMA', 'GZMH', 'GZMB', 'PRF1', 'CX3CR1', 'KLRG1',
                                                     'PDCD1', 'TOX2', 'CXCL13', 'TNFRSF9', 'LAG3',
                                                     'GPR183', 'CCL5', 'GZMA', 'PRF1', 'CCR7', 'SELL', 'IL7R', 'S1PR1', 'ANXA1',
                                                     'ZNF683', 'CXCR6', 'ITGAE', 'ITGA1', 'CD69', 'IL7R', 'XCL1','PRDM1',
                                                     'KLRD1', 'FGFBP2', 'GNLY', 'GZMB', 'GZMH', 'FCGR3A', 'NKG7', 'CX3CR1', 'CD3E',
                                                     'IFI6', 'ISG15','TRDC','STMN1', 'MKI67', 'TYMS', 'TOP2A', 'TUBB',
                                                     'CXCR5', 'BCL6', 'ICA1', 'TOX', 'TOX2', 'IL6ST', 'ICOS','PDCD1','CD200',
                                                     'KLRB1','SLC4A10','NCR3','TRAV1-2'),
                                   B_cell_marker=c("CD19","CR2","MS4A1","BANK1","CD79A","CD79B", 'MZB1', 'HLA-DQA1', 'HLA-DQB1',
                                                   'MZB1', 'IGHG3','IGHG4','IGHA1',
                                                   'IL4R', 'TCL1A','FCER2', 'IGHD',
                                                   'TNFRSF13B', 'GPR183', 'MARCKS','CRIP1',
                                                   'ISG15', 'IFIT3', 'IFITM1', 'IFI44L', 'IFIT1',
                                                   'MARCKSL1', 'RGS13', 'SERPINA9', 'LRMP', 'LMO2', 'NEIL1',
                                                   'STMN1', 'MKI67', 'HMGB2','TUBA1B', 'HMGN2', 'TUBB', 'RRM2'),
                                   Myeloid_cell_marker=c('CD68','CD163','S100A8', 'S100A9', 'LYZ', 'CD14','TPSAB1',
                                                         'CD14','FCGR3B','ITGAM','CCR2','XCR1','LAMP3','CEACAM8', 'S100A8', 'S100A9',
                                                         'LYZ','CD68','CD163','S100A8','S100A9','ITGAX','CD14','FCGR3A','FCGR3B','CX3CR1','CCR2',
                                                         'CD68','CD14','HLA-DRA','HLA-DRB1','MMP9','CTSK','MRC1','CD4','CD8A','CD3E','CCL18','FCN1',
                                                         'CLEC9A','CD1C','THBD','ITGAX',
                                                         'CLEC4C','NRP1','IL3RA','LILRA4',
                                                         'CLEC9A','CADM1','XCR1','BATF3',
                                                         'CD1C', 'CD1E', 'CLEC10A',
                                                         'CPA3', 'TPSAB1', 'MS4A2',
                                                         'CD3E','CD8A','CD79A','MKI67','STMN1','CCL2','KRT19','EPCAM'),
                                   Fibroblast_cell_marker=c("ACTA2","COL1A1",'COL1A2', 'COL3A1', "FAP","DCN","POSTN","VIM","PDGFRB","LUM", 'IL6',
                                                            "ACTA2","COL1A2","PDGFRB",
                                                            "IL6","RGS5","GJA4","MCAM","MYH11","CCL8",
                                                            "POSTN","CTHRC1","COL6A3","MMP14","MMP11","COL5A1","COL5A2","LUM","DCN","VCAN",
                                                            "CST1","IGF1","FBLN1","C3","C7","CXCL1","IGFBP6","SLPI","SAA1",
                                                            "CD74","HLA-DRA","HLA-DRB1",
                                                            "KRT19","KRT8","SAA1"),
                                   Endothelial_cell_marker=c("PECAM1","VWF","CD34","PLVAP","HSPG2","RAMP2",
                                                             'SEMA3G','FBLN5','GJA5','JAG1','PPP1R14A',
                                                             'ACKR1','SELP','SELE','VWF',
                                                             'RGCC',"KIT",
                                                             'JAG1','HES1','ID1','ID2','ID3','ENG','PLVAP','HSPG2','APLNR',
                                                             'NRP1','ENPP2','THY1','FLT1','KDR'),
                                   Plasma_cell_marker = c('TOP2A',"SDC1","CD38","MZB1","IGHA1"),
                                   #Hepatocytes_cell_marker = c('APOE',"ALB","FGA","FGG","HPX","FN1","NKD1"),
                                   Granulocytes_cell_marker = c("ELF3","KRT19","MUC5B"),
                                   Macrophages_cell_marker = c('CD68','CD163',"MARCO","MRC1","MSR1"),
                                   
                                   Epithelial_cell_marker=c("EPCAM", "TG", "TPO","CEACAM6", 'KRT19','KRT7', 'KRT8','TFF3','KRT','ACTA2','CD10',
                                                            'CD13','CD326','CD34','CDH1','CDH4','L1CAM','MCAB51','OCLN','PAX2','SDC1','SNAIL')
)

#' @description
#' 计算细胞中基因富集在特定的细胞类型marker上的得分，从而确定细胞的细胞类型
#' @param FindAllMarkers_res 
#' @param ref_list description 
#' @param topn description
#' @param scale description
Predict_cellType_by_Enrichment <- function(FindAllMarkers_res, ref_list=NULL, topn=20, scale=TRUE,scaleCluster=TRUE){
  if(is.null(ref_list)){
    ref_list <- MainCellType_markers_enrich
  }
  toplist <- FindAllMarkers_res %>% group_by(cluster) %>% top_n(topn,wt = avg_log2FC)
  cell_types = names(ref_list)
  cluster = sort(unique(toplist$cluster))
  inter_mat = matrix(0, ncol = length(cluster), nrow = length(cell_types))
  
  for (i in 1:length(cell_types)) {
    for (j in 1:length(cluster)) {
      pm1 = ref_list[[i]]
      pm2 = toplist[toplist$cluster==cluster[j], "gene", drop=TRUE]
      inter_mat[i,j] = length(intersect(pm1, pm2))
    }
  }
  rownames(inter_mat) = c(cell_types)
  inter_mat <- inter_mat[!rowSums(inter_mat)==0,]
  
  if (is.null(rownames(inter_mat))) {
    print("没有富集到任何细胞类型")
    return(NULL)
  }else{
    norma = sapply(rownames(inter_mat), function(x)length(ref_list[[x]])) / rowSums(inter_mat)
    
    enri = apply(inter_mat,2,function(x)x/norma)
    enri = apply(enri, 2, function(x) x/sum(x))
    if(scale){
      if(scaleCluster){
        enri = scale(enri) #scale: column centering and scaling is performed
      }else{
        enri = t(scale(t(enri)))  #
      }
    }
    rownames(enri) = rownames(inter_mat)
    colnames(enri) = cluster
    enri[which(is.nan(enri))] = 0
    #print(enri)
    if(sum(enri)==0){
      print("没有富集到任何细胞类型")
      return(NULL)
    }else{
      p <- pheatmap::pheatmap(enri,border_color = 'white')
      return(list(enrich = enri,pheatmap = p))
    }
    
  }
}

Extract.predict.cell.type <- function(predict_data){
  predict_data[which(predict_data==0)] <- NA
  m <- list()
  n <- list()
  cluster <- colnames(predict_data)
  ref.celltype <- rownames(predict_data)
  for(ce in cluster) {
    #ce = cluster[7]
    m[[ce]] <- ref.celltype[which.max(x =predict_data[,ce] )]
    if (length(m[[ce]])== 0) {
      m[[ce]] = paste0('unknow_',ce)
    }else{
      m[[ce]] <- substr(m[[ce]],1,nchar(m[[ce]])-7)
    }
  }
  mm <- unlist(m)
  return(mm)
}

#' @description seurat clustering
RunClustering <- function(obj,harmony=FALSE,resolution=0.8,sctransform=FALSE,nfeatures = 3000){
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
  })
  if(sctransform){
    suppressPackageStartupMessages(library(sctransform))
    if(!any(grepl("percent.mt",colnames(obj[[]])))){
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
    }
    obj <- SCTransform(obj,vars.to.regress = "percent.mt", verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose = F)

  }else{
    obj <- NormalizeData(obj) %>% FindVariableFeatures(nfeatures = nfeatures) %>% ScaleData(features = row.names(obj)) 
    
  }
  VariableGenes = VariableFeatures(object = obj)
  IGgenes = VariableGenes[grep("^IG",VariableGenes)]
  TRgenes = VariableGenes[grep("^TR",VariableGenes)]
  ACgenes = VariableGenes[grep("^AC",VariableGenes)]
  ALgenes = VariableGenes[grep("^AL",VariableGenes)]
  HISTgenes = VariableGenes[grep("^HIST",VariableGenes)]
  LINgenes = VariableGenes[grep("^LIN",VariableGenes)]
  MTgenes = VariableGenes[grep("^MT",VariableGenes)]
  RPSgenes = VariableGenes[grep("^RPS",VariableGenes)]
  RPLgenes = VariableGenes[grep("^RPL",VariableGenes)]
  VariableGenes1 = VariableGenes[!VariableGenes %in% c(RPSgenes,TRgenes, IGgenes,ACgenes, ALgenes,HISTgenes,LINgenes,MTgenes,RPLgenes)]
  if(sctransform){
    assay="SCT"
  }else{
    assay="RNA"
  }
  obj[[assay]]@var.features = VariableGenes1
  
  obj <- RunPCA(obj,features = VariableGenes1,verbose = FALSE)
  if(harmony){
    obj <- RunHarmony(obj, group.by.vars = "orig.identNew")
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20)
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20,force.recalc=T) %>% FindClusters(resolution = resolution)
  }else{
    obj <- RunUMAP(obj, dims = 1:20)
    obj <- FindNeighbors(obj, dims = 1:20,force.recalc=T) %>% FindClusters(resolution = resolution)
  }
  return(obj)
}


#' @param color list of colors with name of "seurat_clusters"
subset_ana <- function(object,subset_clusters,assay="RNA",resolution = 0.8,marker.list = NULL,color=NULL,Findmarker=FALSE){
  obj_sub <- subset(object,seurat_clusters %in% subset_clusters)
  if(grepl("ATAC|peaks",assay)){
    obj_sub <- runClustering.scATAC(obj_sub,assay = assay,TopFeatures.cutoff = 'q5',resolution = resolution)
    
  }else if(grepl("RNA|SCT",assay)){
    obj_sub <- RunClustering(obj_sub,resolution = resolution)
  }
  p1=DimPlot(obj_sub, reduction = "umap", group.by = 'seurat_clusters', label = TRUE, label.size = 5, cols = color[["seurat_clusters"]],
             repel = TRUE) + ggtitle("Clusters")#+ NoLegend() 
  features_plot <- paste0(c("nCount_","nFeature_"),assay)
  pv <- VlnPlot_v2(obj_sub,group_by="seurat_clusters",features.plot= features_plot,
                   color_pattern=color[["seurat_clusters"]],ncol=1)
  
  
  if(is.null(marker.list)){
    marker.list <- list(immune_cell = 'PTPRC',
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
  
  marker_df <- do.call(rbind,lapply(1:length(marker.list), function(x,marker.list){
    celltype=names(marker.list)[x]
    mks <- as.character(marker.list[[x]])
    df <- data.frame(Gene=mks,CellType=rep(celltype,length(mks)))
    return(df)
  },MainCellType_markers_plot))
  p2 <- cellType.Plot(obj_sub,violinMarker=marker_df,dotMarker=marker_df,col_palette=color[["seurat_clusters"]])
  if(Findmarker){
    features.use <- rownames(obj_sub)
    features.use <- features.use[!grepl("^IG|^TR|^AC|^AL|^HIST|^LIN|^MT|^RPS",features.use)]
    markerGene <- FindAllMarkers(obj_sub,features=features.use, only.pos = TRUE,assay = assay)
  }else{markerGene=NULL}
  return(list(object=obj_sub,
              Figure=list(DimPlot=p1,VlnPlot=pv,Vln_Dot=p2),
              Findmarker_res=markerGene ))
}


suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(dplyr)
  library(Seurat)
  library(Signac)
  library(patchwork)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

#' @description 
#'      Extract.CellRanger.scMultiomic.10x, then analysi by seurat
#' @param File.dir path
#' @param CallingPeak calling peaks or not
#' 
Extract.CellRanger.scMultiomic.10x <- function(File.dir,CallPeaks=F,group.by = NULL,macs2.path = NULL){
  
  sampleName = names(File.dir)
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  #seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  seqlevelsStyle(annotation) <- 'UCSC'
  genome(annotation) <- 'hg38'
  
  #######
  results = lapply(sampleName, function(i){
    
    cat(format(Sys.time(), "[%b-%d,%H:%M:%S] "));
    cat("Reading the Sample of", i, "......");
    
    # load the RNA and ATAC data
    counts_path = file.path(File.dir[[i]],"filtered_feature_bc_matrix.h5")
    fragment_path = file.path(File.dir[[i]],"atac_fragments.tsv.gz");
    counts <- Read10X_h5(counts_path)
    
    # create a Seurat object
    temp.seurat <- CreateSeuratObject(counts = counts$`Gene Expression`,
                                      project = i, 
                                      assay = "RNA")
    temp.seurat$sampleID = i
    # create ATAC assay and add it to the object
    temp.seurat[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,
                                                  sep = c(":", "-"),
                                                  fragments = fragment_path,
                                                  annotation = annotation)
    # peak calling
    if(CallPeaks){
      cat('peaking calling')
      DefaultAssay(temp.seurat) <- 'ATAC'
      peaks = CallPeaks(temp.seurat,group.by=group.by,macs2.path=macs2.path)
      peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
      peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
      # quantify counts in each peak
      macs2_counts <- FeatureMatrix(fragments = Fragments(temp.seurat), features = peaks, cells = colnames(temp.seurat))
      temp.seurat[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts, fragments = fragment_path, annotation = annotation )
    }
    #产生线粒体基因count占比信息，方便未来过滤使用
    DefaultAssay(temp.seurat) <- 'RNA'
    temp.seurat[["percent.mt"]] = PercentageFeatureSet(temp.seurat, pattern = "^MT-")
    
    cat( ncol(temp.seurat),'cells,',nrow(temp.seurat),'features.\n')
    return(temp.seurat)
  })
  names(results) = sampleName;
  return(results)
}

#' @description 
#' @param obj seuarat 对象
#' 
Calculate.QC.index.scMulti <- function(obj,assay="peaks"){
  DefaultAssay(obj) <- assay
  obj <- NucleosomeSignal(object = obj)  
  obj <- TSSEnrichment(object = obj, fast = TRUE,assay=assay)
  obj$log10_nFrags_ATAC <- log10(obj$nFeature_ATAC)
  obj$log10_nCount_ATAC <- log10(obj$nCount_ATAC)
  obj$log10_nFeature_RNA <- log10(obj$nFeature_RNA)
  obj$log10_nCount_RNA <- log10(obj$nCount_RNA)
  #obj$log10_nCount_peaks <- log10(obj$nCount_peaks)
  #obj$log10_nFrags_peaks <- log10(obj$nFeature_peaks)
  return(obj)
}

#' @description
#' A short description...
#' 
runClustering.scMulti <- function(obj,algorithm = 3,dims.list = list(1:40, 2:30),
                                  assay.rna="RNA",assay.atac="peaks"){
  DefaultAssay(obj) <- assay.rna
  obj <- SCTransform(obj, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  DefaultAssay(obj) <- assay.atac
  obj <- FindTopFeatures(obj, min.cutoff = 5)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  obj <- FindMultiModalNeighbors(object = obj,reduction.list = list("pca","lsi"),dims.list = dims.list,
                                 modality.weight.name = "RNA.weight",verbose = T)
  obj <- RunUMAP(object = obj,nn.name = "weighted.nn",reduction.name = "wnn.umap", seed.use = 34,reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = algorithm, verbose = FALSE,random.seed = 343)
  
  return(obj)
  
}

#' @description 
#' @param obj seuarat object
#' 
Calculate.QC.index.scRNA <- function(obj,assay="RNA"){
  DefaultAssay(obj) <- assay
  if(!any(grepl("percent.mt",colnames(obj[[]])))){
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  }
  if(!any(grepl("percent.ribo",colnames(obj[[]])))){
    # Ribosomal
    obj <- PercentageFeatureSet(obj, "^RP[SL]|^rp[sl]", col.name = "percent.ribo")
     }
  if(!any(grepl("percent.hb",colnames(obj[[]])))){
    # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
    obj <- PercentageFeatureSet(obj, "^HB[^(P|E|S)]|^hb[^(p|e|s)]", col.name = "percent.hb")
  }
  if(!any(grepl("log10_nFeature_RNA",colnames(obj[[]])))){
    obj$log10_nFeature_RNA <- log10(obj$nFeature_RNA)
  }
  if(!any(grepl("log10_nCount_RNA",colnames(obj[[]])))){
    obj$log10_nCount_RNA <- log10(obj$nCount_RNA)
  }
  return(obj)
}




