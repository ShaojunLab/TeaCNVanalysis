##Fig.2ab UMAP


library(GenomicRanges)
library(Signac)
library(TeaCNV)
library(gridExtra)
cols_Palette <- c("#CCCCCC","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
setwd(work_dir)
source("./Code/Function/funs_dataProcess.R")
source("./Code/Function/atac_visualizeFun.R")
source("./Code/Function/colorPalette.R")
source("./Code/Function/fun_Grange.R")
outdir <- work_dir




for(ID in c("ccRCC2","ccRCC3","ccRCC4")){
	if(ID %in% c("ccRCC3","ccRCC4")){
		rdsFile <- "./Rdata/Kidney_scMutiOmic_2Sample.rds"
		obj <- readRDS(rdsFile)
	}else if(ID %in% c("ccRCC2")){
		rdsFile <- "./Rdata/Kidney_scATAC_2Sample.rds"
		obj <- readRDS(rdsFile)
	}


	#ID = "ccRCC3"
	outdir_clt <- paste0(outdir,'/',ID)
	CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
	CNVres <- readRDS(CNVresFile_sc)
	cellinfo <- CNVres$cellinfo
	cellinfo <- cellinfo[!is.na(cellinfo$clone),]
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

	umap_coords_rna <- Embeddings(obj_sub, "umap.rna")
	umap_coords_atac <- Embeddings(obj_sub, "umap.atac")
	obj_sub@meta.data <- cbind(obj_sub@meta.data, umap_coords_rna,umap_coords_atac)

	mat_r <- CNVres$cellbinRatio_raw
	mat_r <- mat_r[,rownames(cellinfo)]

	if(ID=="ccRCC3"){
	  mat_r <- mat_r[grepl("chr8|chr16",rownames(mat_r)),]
	  CNV_segs <- data.frame(Region=c("chr8_231970_43141576","chr16_46688301_90082747","chr8_47260174_145052768","chr16_9930_31874205"))#from "density_histogram_segs.r"

	}
	if(ID=="ccRCC4"){
	  mat_r <- mat_r[grepl("chr8|chr20|chr22",rownames(mat_r)),]
	  CNV_segs <- data.frame(Region=c("chr20_289957_64287454","chr22_17158443_50784147"))#from "density_histogram_segs.r"

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

	#Fig.
	pd <- DimPlot(new_obj, reduction = "umap", group.by = "clone",pt.size =1.2,label.size = 5,label = F,repel = T,cols=color_r)+ ggtitle("CNV")#&NoLegend()
	pd+p_atac+p_rna
	ggsave(paste0(outdir_clt,"/umap_Epi_clone.pdf"),pd+p_atac+p_rna,height=4.5,width =13.5)

	#fig.
	markers <- c("EPCAM","KRT19","KRT18","KRT8")
	DefaultAssay(obj_sub) <- "SCT"
	existing_genes <- markers[markers %in% rownames(obj_sub)]
	exp_mt<- GetAssayData(obj_sub, assay = "RNA", slot = "data")
	avg_exp <- colMeans(exp_mt[existing_genes, ])
	obj_sub$TumorMarkers_avg <- avg_exp

	new_obj<- AddMetaData(new_obj,obj_sub@meta.data[,"TumorMarkers_avg",drop=F])

	if(ID%in% c(c("ccRCC3","ccRCC4"))){min_cutoff="q20";max_cutoff = "q80"}else{min_cutoff="q10";max_cutoff = "q90"}
	pm <- FeaturePlot(object = new_obj, features = "TumorMarkers_avg",ncol=1,pt.size =1.2,
	   min.cutoff = min_cutoff, max.cutoff = max_cutoff)& scale_colour_gradientn(colours = c("#F1F2F2","yellow","red")) 
	pm <- lapply(pm, function(x) x  +
	  theme(
	    axis.ticks = element_blank(),
	    axis.text = element_blank(),
	    axis.title = element_blank(),
	    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
	    axis.line = element_blank(),legend.position = "bottom"
	  )) 
	pm <- grid.arrange(grobs = pm)
	ggsave(paste0(outdir_clt,"/umap_TumorMarkers_avg.pdf"),pm,height=4.5,width =4.5)


	umap_coords <- Embeddings(new_obj, "umap")
colnames(umap_coords) <- c("cnvUMAP_1","cnvUMAP_2")

obj_sub<- AddMetaData(obj_sub,umap_coords)

#
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
    start = as.integer(start),
    end = as.integer(end)
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
  counts_mat_seg <- counts_mat[peaks_key,]
  nCount_seg <- colSums(counts_mat_seg)

  tb_nCount <- rbind(tb_nCount,nCount_seg)
}
tb_nCount <- as.data.frame(t(tb_nCount))
colnames(tb_nCount) <-CNV_segs$Region

# saveRDS(obj_sub,paste0(outdir_clt,"/SeuratObj_final_clone.rds"))
new_obj<- AddMetaData(new_obj,tb_nCount)
# saveRDS(new_obj,paste0(outdir_clt,"/SeuratObj_CNVumap.rds"))

##Fig.
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
width <- ifelse(length(plot.list)>4,16,4*length(plot.list))
combined_plot <- grid.arrange(
  grobs = plot.list, 
  ncol = ncol
)
ggsave(paste0(outdir_clt, "/umapcnv_CNVregion_atac.pdf"), combined_plot, height =4.5, width = width)



}
