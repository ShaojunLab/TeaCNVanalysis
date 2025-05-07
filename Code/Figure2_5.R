#Figure2,5

suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(ggsci)
  library(Signac)
  library(Seurat)
  library(TeaCNV)
 })

work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
setwd(work_dir)
source('./Code/Function/atac_visualizeFun.R')
#plotdir <- '~/Library/Mobile Documents/com~apple~CloudDocs/code/code_CNV/TeaCNVAnalysis'
plotdir <- './'
# setwd(plotdir)
color_celltype <- c(Endothelial="#CCFFCC",Epithelial="#C5A48A",Immune="#CCCCCC",Stromal="#CC99FF")
cols_Palette = c("#CCCCCC","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")

cnv_plots_comb <- function(ID,outres,
	cols_Palette = c("#CCCCCC","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780"),
	plotdir="./"
	){
	#(1)
	TeaCNV::plot_combine_seg(outres$clonalest,ylim=NULL,
	    outplot_name=paste0(ID,".clonalCNV_final.pdf"),show_dots=FALSE,
	    outdir=plotdir)
	#(2)
	CNmt <- do.call(cbind,lapply(sort(names(outres$clonalest)),function(cluster){
	  clonal_res <- outres$clonalest[[cluster]]$input_BinSegRatio
	  integerCNV <- outres$clonalest[[cluster]]$seg.dat
	  integerCNV <- integerCNV[,c("segName","relativeCN","integerCN")]
	  clonal_res <- clonal_res[,!grepl("relativeCN|integerCN",colnames(clonal_res))]
	  clonal_res <- left_join(clonal_res,integerCNV,by="segName")
	  rownames(clonal_res) <- clonal_res$binID
	  return(clonal_res$integerCN)
	}))
	rownames(CNmt) <- rownames(outres$clonalest[[1]]$input_BinSegRatio)
	colnames(CNmt) <- sort(names(outres$clonalest))
	clone_info <- data.frame(row.names=colnames(CNmt),clone=colnames(CNmt))
	color_r <- cols_Palette[1:length(unique(clone_info$clone))]
	    names(color_r) <- sort(unique(clone_info$clone))
	    left_anno_cols <- list()
	    left_anno_cols[["clone"]] <- color_r
	    height <- ifelse(ncol(CNmt)>2,0.35*(ncol(CNmt))+0.5,ifelse(ncol(CNmt)==1,1.2,1.5))
	    p_cloneCN <- heatmap4peakMt(mat=CNmt,
	                          meta_info=clone_info,
	                          sep_by="-",
	                          outdir= plotdir,
	                          value.type="CNV",
	                          clust_rows=F,
	                          show_legend_row = T,
	                          legend_titles="integer CN",
	                          fileout_name=paste0(ID,".heatmap_cloneCNA"),
	                          col_list=left_anno_cols,
	                          column.title = NULL,
	                          width=10,height=height,device="pdf") 

	#(3)
	plt_mt2 <- outres$cellbinRatio_raw
	plt_mt2 <- smooth_and_denoise(plt_mt2,window=60)
	left_anno_cols <- list()
	cellMeta2 <- outres$cellinfo[,c("clone"),drop=F]
	cellMeta2 <- na.omit(cellMeta2)
	color_r2 <- cols_Palette[1:length(unique(cellMeta2$clone))]
	names(color_r2) <- sort(unique(cellMeta2$clone))
	left_anno_cols[['clone']] <- color_r2
	p_scRatio <- heatmap4peakMt(mat=plt_mt2[,rownames(cellMeta2)],
	                         meta_info=cellMeta2,
	                         sep_by="-",
	                         outdir= plotdir,value.type="ratio",
	                         clust_rows=F,clustering_method_rows = "ward.D2",
	                        show_legend_row = T,
	                         fileout_name=paste0(ID,".heatmap_scRatioRaw"),
	                         col_list=left_anno_cols,
	                         width=10,height=5,device="pdf")
}
smooth_median_by_group <- function(df,k = 10,adjust=1) {
	df %>%
	  group_by(integerCN) %>%
	  mutate(
	    smoothed_binRatio = rollmean(binRatio, k = k, fill = NA, align = "right")  # 计算步长为k的均值
	  ) %>%
	  filter(!is.na(smoothed_binRatio)) %>% 
	  # mutate(
	  #   relativeCN_mean = median(smoothed_binRatio, na.rm = TRUE)  # 计算平滑后的中值
	  # )%>%   
	  mutate(
	     relativeCN_mean = {
	      d <- density(smoothed_binRatio,adjust=adjust)
	      max_density_idx <- which.max(d$y)
	      d$x[max_density_idx]
	    }
	  )%>%
	  as.data.frame()
}
DEseg <- function(outres){
	clonalest <- outres$clonalest
	res <- c()
	for(i in 1:length(clonalest)){
	  df_bin <- clonalest[[i]]$input_BinSegRatio
	  df_bin <- df_bin[,c("binID","Chromosome","Start","End","binRatio","length_bin","SegMean","segName")]
	  df_seg <- clonalest[[i]]$seg.dat
	  df_seg <- df_seg[,c("segName","relativeCN","integerCN")]
	  df_bin <- left_join(df_bin,df_seg,by=c("segName"))
	  df_bin$clone <- names(clonalest)[i]
	  df_bin$wei <-  df_bin$length_bin/sum(df_bin$length_bin)
	  res <- rbind(res,df_bin)
	}
	#select differential CNV region accross clones
	res_cnv <- unique(res[,c("Chromosome","segName","integerCN","clone")])
	res_cnv <- na.omit(res_cnv[res_cnv$integerCN != "2",])
	res_cnv$segstart <- as.numeric(sapply(strsplit(res_cnv$segName,"_|-|:"),"[",2))
	res_cnv$segend <- as.numeric(sapply(strsplit(res_cnv$segName,"_|-|:"),"[",3))
	filtered_res_cnv <- res_cnv %>%
	group_by(Chromosome, integerCN) %>%
	filter(n_distinct(clone) != length(unique(res_cnv$clone))) %>%
	ungroup()%>%
	as.data.frame()
	merged_df <- filtered_res_cnv %>%
	    dplyr::group_by(Chromosome)%>%
	    dplyr::mutate(newSeg = cumsum(integerCN != lag(integerCN, default = first(integerCN)) | segstart >= lag(segend + 2e6, default = first(segstart)))) %>%
	    ungroup()%>%
	    dplyr::group_by(Chromosome,newSeg, integerCN) %>%
	    summarise(
	      start = first(segstart),    
	      end = last(segend)          
	    ) %>%
	    ungroup()%>%
	    as.data.frame()
	merged_df <- merged_df[!duplicated(merged_df[, c("Chromosome", "integerCN","start","end")]), ]
return(list(res=res,merged_df=merged_df))
}
custom_theme <-
    theme_classic()+ 
    theme(plot.background=element_blank(),
          legend.position='right',
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(color="black",size=15),
          axis.text=element_text(color="black",size=12),
          axis.line.y.right = element_blank(),
          axis.text.y.right=element_blank(),
          axis.ticks.y.right=element_blank(),
          axis.ticks=element_line(color="black",linewidth=0.5),
          plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_blank()
        )


## ccRCC1
ID = "ccRCC1"
rdsFile <- "./Rdata/Kidney_scATAC_2Sample.rds"
obj <- readRDS(rdsFile)
#Fig.2b
color_list <- list(sampleName=c("#FFCC00","#009999"),
                   CellType = color_celltype)
group_by <- c("sampleName","CellType") 
pd <- DimPlot.multi(obj,group_by = group_by,color_list=color_list,legend_position = "right")
ggsave(paste0(plotdir,"/DimPlot_",group_by[1],".pdf"),pd,width =10,height = 5)

#Fig.2c
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)

##Supplementary Fig. 4
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,"_immune.rds")
outres2 <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres2,plotdir=plotdir)

CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,"_Endo.rds")
outres3 <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres3,plotdir=plotdir)

##Supplementary Fig. 3
DefaultAssay(obj) <- "peaks"
Idents(obj)<- "CellType"
library(Signac)
library(ggplot2)
markers <- c("PTPRC","PECAM1","CD34","COL3A1","COL1A2","EPCAM","KRT8","KRT18")
Coverageplt <- list()
for(gene in markers){
	#(1) Plotting aggregated signal
	cov_plot <- CoveragePlot(
	  object = obj,
	  region = gene,
	  annotation = TRUE,
	  peaks = F,
	  links = F
	)+scale_fill_manual(values=color_celltype)
	
	ggsave(paste0(plotdir,"/",ID,"_Coverageplt_",gene,".pdf"),cov_plot,width =5,height = 5)

}



##ccRCC2
ID = "ccRCC2"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
#Fig.2d,e
cnv_plots_comb(ID,outres,plotdir=plotdir)

##Supplementary Fig. 5
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,"_immune.rds")
outres2 <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres2,plotdir=plotdir)

CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,"_Endo.rds")
outres3 <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres3,plotdir=plotdir)


##.ccRCC3
ID = "ccRCC3"
rdsFile <- "./Rdata/Kidney_scMutiOmic_2Sample.rds"
obj <- readRDS(rdsFile)


##Supplementary Fig. 2

DefaultAssay(obj) <- "RNA"
FeaturePlot(obj,features = markers,
                   reduction = "umap",pt.size = 0.1,
                   max.cutoff = 'q95',ncol = 4) 
DefaultAssay(obj) <- "peaks"
Idents(obj)<- "CellType"
library(Signac)
library(ggplot2)
Coverageplt <- list()
for(gene in markers){
	#(1) Plotting aggregated signal
	cov_plot <- CoveragePlot(
	  object = obj,
	  region = gene,
	  annotation = TRUE,
	  peaks = F,
	  links = F
	)+scale_fill_manual(values=color_celltype)
	
	ggsave(paste0(plotdir,"/",ID,"_Coverageplt_",gene,".pdf"),cov_plot,width =5,height = 5)

}



#Fig.2a
group_by <- c("sampleName","CellType") 
pd <- DimPlot.multi(obj,group_by = group_by,color_list=color_list,legend_position = "right")
ggsave(paste0(plotdir,"/DimPlot_",group_by[1],".pdf"),pd,width =10,height = 5)
#Fig.2f,g
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
##
cnv_plots_comb(ID,outres,plotdir=plotdir)

##ccRCC4
ID = "ccRCC4"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
#Fig.2h,i
cnv_plots_comb(ID,outres,plotdir=plotdir)


#Fig.5
##BRCA
ID = "BRCA"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)

##PDAC
ID = "PDAC"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)


##HNSCC
ID = "HNSCC"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)

##CRC
ID = "CRC"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)

##Fig.OV
ID = "OV"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)




library(GenomicRanges)
library(Signac)
library(TeaCNV)
library(gridExtra)

outdir0 <- "./GSE240822"
if(!file.exists(outdir0)){dir.create(outdir0,recursive=T)}
IDs <- c("BRCA","PDAC","HNSCC","CESC","OV")
SampleInfo <- data.frame(cancer=c("BRCA","PDAC","HNSCC","CRC","OV"),
            sample=c("HT243B1-S1H4","PM581P1-T1","P5590-N1","HT250C1-Th1K1","VF027V1-S2"))

for(i in 1:5){
  ID <- SampleInfo[i,1]
  outdir <- paste0(outdir0,"/",ID);if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  sampleID <- SampleInfo[i,2]
  print(ID)
  print(sampleID)
  outdir_clt <- paste0(outdir,"/atac2cnv/",sampleID)

  datadir <- paste0("./GSE240822/",ID)
  inputdata <- paste0(datadir,"/",ID,"_AllCell_Final_SuerObj.RDS")
  obj <- readRDS(inputdata)

  CNVres <- readRDS(paste0(outdir_clt,"/final.CNVres.rds"))
  cellinfo <- CNVres$cellinfo
  cellinfo <- cellinfo[!is.na(cellinfo$clone),]
  dim(cellinfo)
  color_r <- cols_Palette[1:length(unique(cellinfo$clone))]
  names(color_r) <- sort(unique(cellinfo$clone))

  obj_sub <- obj[,rownames(cellinfo)]
  obj_sub<- AddMetaData(obj_sub,cellinfo)
  Idents(obj_sub)<- obj_sub$clone

  #UMAP RNA
  DefaultAssay(obj_sub) <- "RNA"
  obj_sub <- SCTransform(obj_sub, verbose = FALSE) %>% FindVariableFeatures(nfeatures = 3000)
  VariableGenes = VariableFeatures(object = obj_sub)
  removeG <- VariableGenes[grep("^IG|^TR|^AC|^AL|^HIST|^LIN|^MT|^RPS",VariableGenes)]
  VariableGenes  <- VariableGenes[!VariableGenes%in%removeG]
  obj_sub@assays$SCT@var.features = VariableGenes
  obj_sub <- obj_sub%>% RunPCA(features = VariableGenes,verbose = FALSE) %>% 
    RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  p_rna <- DimPlot(obj_sub, reduction = "umap.rna", group.by = "clone",pt.size =1.2,label.size = 5,label = F,repel = T,cols=color_r)+ ggtitle("RNA")#&NoLegend()


  #UMAP ATAC
  DefaultAssay(obj_sub) <- "peaks"
  obj_sub <- RunTFIDF(obj_sub)
  obj_sub <- FindTopFeatures(obj_sub, min.cutoff = 'q5')
  obj_sub <- RunSVD(obj_sub)
  obj_sub <- RunUMAP(object = obj_sub, reduction = 'lsi', dims = 2:30,reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  p_atac <- DimPlot(obj_sub, reduction = "umap.atac", group.by = "clone",pt.size =1.2,label.size = 5,label = F,repel = T,cols=color_r)+ ggtitle("ATAC")#&NoLegend()
  #ggsave(paste0(outdir_clt,"/umapATAC_Epi_clone.pdf"),p_atac,height=4.5,width =4.5)

  umap_coords_rna <- Embeddings(obj_sub, "umap.rna")
  umap_coords_atac <- Embeddings(obj_sub, "umap.atac")
  obj_sub@meta.data <- cbind(obj_sub@meta.data, umap_coords_rna,umap_coords_atac)

  mat_r <- CNVres$cellbinRatio_raw
  mat_r <- mat_r[,rownames(cellinfo)]

  if(ID=="PDAC" & sampleID=="PM581P1-T1"){
    mat_r <- mat_r[grepl("chr18|chr19",rownames(mat_r)),]
    CNV_segs <- data.frame(Region=c("chr18_21599651_80247759","chr19_29212136_40571061"))#from "CNVvalue_ano"

  }
  if(ID=="BRCA" & sampleID=="HT243B1-S1H4"){
    mat_r <- mat_r[grepl("chr1|chr10|chr11|chr16|chr20",rownames(mat_r)),]

    CNV_segs <- data.frame(Region=c("chr11_108663863_134332689","chr20_49045582_56468036"))#from "CNVvalue_ano"

  }
  if(ID=="HNSCC" & sampleID=="P5590-N1"){
    # mat_r <- mat_r[grepl("chr3|chr4|chr5|chr13|chr22",rownames(mat_r)),]
    CNV_segs <- data.frame(Region=c("chr3_143970877_198081878","chr5_52799370_135142136","chr13_19633124_114314884"))#from "CNVvalue_ano"

  }

  if(ID=="CRC" & sampleID=="HT250C1-Th1K1"){
    mat_r <- mat_r[grepl("chr2",rownames(mat_r)),]

    CNV_segs <- data.frame(Region=c("chr2_44660_242089223"))#from "CNVvalue_ano"

  }

  if(ID=="OV" & sampleID=="VF027V1-S2"){
    # mat_r <- mat_r[grepl("chr1|chr7|chr19",rownames(mat_r)),]
    # CNV_segs <- data.frame(Region=c("chr1_36322681_54801704","chr19_29211648_39032790"))#from "CNVvalue_ano"
    # CNV_segs <- data.frame(Region=c("chr1_999532_121185055","chr7_148953_56106829","chr8_730578_43141640","chr19_39236582_47484450"))#from "CNVvalue_ano"
    CNV_segs <- data.frame(Region=c("chr1_999532_121185055","chr1_165598925_165630123","chr5_160008300_181160654","chr7_148953_56106829","chr8_730578_43141640","chr19_39236582_47484450"))#from "CNVvalue_ano"

  }

  new_assay <- CreateChromatinAssay(
    counts  = mat_r,
    sep = c("-", "-"),
    genome = "GRCh38",
    min.cells = 1
  )
  new_obj <- CreateSeuratObject(counts = new_assay, assay = "subPeaks")
  new_obj <- FindTopFeatures(new_obj, min.cutoff = 'q5')
  new_obj <- ScaleData(new_obj,features = rownames(new_obj) )
  Ncells <- ncol(new_obj)
  if(Ncells<50 & Ncells>1){
    npcs.pca=floor(Ncells/2)
    dims.max=min(npcs.pca,15)
    n.neighbors = min(30,Ncells)
  }else{npcs.pca=50;dims.max=15;n.neighbors=30}
  new_obj <- RunPCA(new_obj, features = VariableFeatures(object = new_obj),npcs=npcs.pca)
  new_obj <- RunUMAP(new_obj, reduction = 'pca',dims = 1:dims.max,n.neighbors=n.neighbors)
  new_obj<- AddMetaData(new_obj,cellinfo)

  pd <- DimPlot(new_obj, reduction = "umap", group.by = "clone",pt.size =1.2,label.size = 5,label = F,repel = T,cols=color_r)+ ggtitle("CNV")#&NoLegend()
  pd+p_atac+p_rna
  ggsave(paste0(outdir_clt,"/umap_Epi_clone.pdf"),pd+p_atac+p_rna,height=4.5,width =13.5)

  markers <- c("EPCAM","KRT19","KRT18","KRT8")
  DefaultAssay(obj_sub) <- "RNA"
  existing_genes <- markers[markers %in% rownames(obj_sub)]
  exp_mt<- GetAssayData(obj_sub, assay = "RNA", slot = "data")
  avg_exp <- colMeans(exp_mt[existing_genes, ])
  obj_sub$TumorMarkers_avg <- avg_exp

  new_obj<- AddMetaData(new_obj,obj_sub@meta.data[,"TumorMarkers_avg",drop=F])
  pm <- FeaturePlot(object = new_obj, features = "TumorMarkers_avg",ncol=1,pt.size =1.2,
     min.cutoff = "q10", max.cutoff = "q90")& scale_colour_gradientn(colours = c("#F1F2F2","yellow","red")) 
  pm <- lapply(pm, function(x) x  +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),legend.position = "bottom"
    )) 

  umap_coords <- Embeddings(new_obj, "umap")
  colnames(umap_coords) <- c("cnvUMAP_1","cnvUMAP_2")
  obj_sub<- AddMetaData(obj_sub,umap_coords)

  #计算CNV区间的ATAC nCount
  DefaultAssay(obj_sub) <- "peaks"
  counts_mat <- GetAssayData(obj_sub[["peaks"]], slot = "counts")
  peak.gr <- StringToGRanges(rownames(obj_sub))



  CNV_segs <- CNV_segs %>%
  tidyr::separate(
    col = Region,
    into = c("chr", "start", "end"),
    sep = ":|-|_",   
    remove = FALSE
  ) %>%
  dplyr::mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  dplyr::select(chr, start, end, everything())

  tb_nCount <- c()
  for(s in 1:nrow(CNV_segs) ){
    CNV_segs_gr <- GRanges(
      seqnames = CNV_segs$chr[s],
      ranges = IRanges(start = CNV_segs$start[s], end = CNV_segs$end[s]),
      CNV = CNV_segs$Region[s]
    )
    overlapping_peaks <- subsetByOverlaps(peak.gr, CNV_segs_gr)
    peaks_key <- as.character(overlapping_peaks)
    peaks_key <- gsub(":","-",peaks_key)
    counts_mat_seg <- counts_mat[peaks_key,,drop=F]
    counts_mat_seg <- as.matrix(counts_mat_seg)
    nCount_seg <- colSums(counts_mat_seg)

    tb_nCount <- rbind(tb_nCount,nCount_seg)
  }
  tb_nCount <- as.data.frame(t(tb_nCount))
  colnames(tb_nCount) <-CNV_segs$Region

  saveRDS(obj_sub,paste0(outdir_clt,"/SeuratObj_final_clone.rds"))
  new_obj<- AddMetaData(new_obj,tb_nCount)
  saveRDS(new_obj,paste0(outdir_clt,"/SeuratObj_CNVumap.rds"))

  pp <- FeaturePlot(object = new_obj, features = CNV_segs$Region,ncol=4,pt.size =1.2,
     min.cutoff = "q5", max.cutoff = "q95")& scale_colour_gradientn(colours = c("#F1F2F2","yellow","red")) 
  plot.list <- lapply(pp, function(x) x  +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),legend.position = "bottom"
    ))

  ncol <- ifelse(length(plot.list)>4,4,length(plot.list))
  combined_plot <- grid.arrange(
    grobs = c(pm,plot.list),
    ncol = ncol+1
  )
  width <- ifelse(length(plot.list)>4,16,4*(length(plot.list)+1))
  height <- ifelse(length(plot.list)>4,4.5*(ceiling(length(plot.list)/4)),4.5)
  ggsave(paste0(outdir_clt, "/umapcnv_CNVregion_atac.pdf"), combined_plot, height =height, width = width)
}


