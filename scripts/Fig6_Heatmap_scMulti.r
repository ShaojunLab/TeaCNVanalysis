#This is for visualization of scMultiomic of HCC
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(tidyverse)
  library(hrbrthemes)
  library(viridis)
  library(gridExtra)
  library(ggpubr)
  library(RColorBrewer)
  library(EnsDb.Hsapiens.v86)
  library(ggsci)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(ComplexHeatmap)
  library(patchwork)
  library(circlize)
  library(parallelDist)
  library(paletteer)
  library(TeaCNV)
})
script_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/code/code_CNV/github/TeaCNV/"
setwd(script_path)
# source("colorPalette.R")
# source("fun_Grange.R")
# source("fun_pseudoBulk_mat.R")
# source("atac_visualizeFun.R")
# source("mydataProcess.R")
#source("~/Library/Mobile Documents/com~apple~CloudDocs/code/HCC/fun_gene_atac_cor.R")
# source("~/Library/Mobile Documents/com~apple~CloudDocs/code/R/fun_Heatmap_PeakGeneLinks.R")
cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
cols13 <- c("#A9011B","#E4A826","#D2DCAD","#DCD6B2","#94BEA7","#4E7989","#75A0AE","#80944E","#547DB1","#9055A2","#D43B36","#F4C28F","#BAAFD1")
cols_Palette_clone <- c("#CCCCCC",
    "#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
names(cols_Palette_clone) <- as.character(seq(1:9))
cols_heat <- c("#8c510a", "#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e")



data_path = paste0(work_dir, '/Rdata/')
inputdata <- paste0(data_path,"/HCC_MultiOmic_12sam_clean.rds")
obj <- readRDS(inputdata)
gtf_file <- "~/Library/Mobile Documents/com~apple~CloudDocs/reference/genes.gtf.gz"
annotations_all <- import(gtf_file)
annotations <- annotations_all[annotations_all$type=="gene",]
work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
cnvdir <- paste0(work_dir,'/TeaCNV_HCC_MultiOmic')
sampleIDs <- list.files(cnvdir)
outdir <- paste0(work_dir,"/Figures")
setwd(outdir)

#Figure legend: Heatmap showing chromatin accessibility and gene expression of N significantly linked peak-gene pairs (columns: left, CRE accessibility; right, linked gene expression). Rows were clustered using k-means clustering (k = 5).
tb <- as.data.frame(table(obj$sampleID,obj$tissue))
tb <- tb %>%
dplyr::filter(Var2=="HCC_tumor")%>%
dplyr::arrange(desc(Freq))
sampleNm <- data.frame(sampleID=unique(tb$Var1),sampleID_new= c(paste0("HCC",seq(1:10)),paste0("HCCn",seq(1:2))))
names(cols13) <- sampleNm$sampleID_new

color_list <- list(dataset=unname(cols13),
                   sampleID_new = cols13[1:12])



VariableGenes <- obj[["SCT"]]@var.features
if(length(Links(obj[["peaks"]]))==0){
    obj <- RegionStats(object = obj,assay = 'peaks',genome = BSgenome.Hsapiens.UCSC.hg38)
    obj<- LinkPeaks(
      object = obj,
      peak.assay = "peaks",
      expression.assay = "SCT",
      genes.use = VariableGenes,
      pvalue_cutoff = 0.05, #default
      score_cutoff = 0.05, #coef.result>score_cutoff, default 0.05
      distance = 5e+05 #default  5e+05
    )

}
gplinks <- Links(obj[["peaks"]])
gplinks <- as.data.frame(gplinks)
write.table(gplinks,paste0(outdir,"/Links_500kb_VarGene.txt"),col.names = T,row.names = F,quote=F,sep="\t")
# gplinks <- read.table(paste0(outdir,"/Links_500kb_VarGene.txt"),header=T,sep="\t")
gplinks[gplinks$gene %in% c("VEGFR","CCND1"),]
which(VariableGenes%in% c("VEGFR","CCND1"))

sigLink <- gplinks[gplinks$pvalue<0.05&gplinks$score>0.05,,drop=F]
sigLink <- sigLink[order(sigLink$pvalue,sigLink$score,decreasing = c(FALSE,TRUE)),]
sigLink <- sigLink[!(sigLink$seqnames%in% c("chrX","chrY")),]
write.table(sigLink,paste0(outdir,"/sigLinks_500kb_VarGene.txt"),col.names = T,row.names = F,quote=F,sep="\t")
# sigLink <- read.table(paste0(outdir,"/sigLinks_500kb_VarGene.txt"),header=T,sep="\t")
# sigLink <- sigLink[order(sigLink$pvalue,sigLink$score,decreasing = c(FALSE,TRUE)),]
# sigLink <- sigLink[!(sigLink$seqnames%in% c("chrX","chrY")),]

sigLink[sigLink$gene %in% c("VEGFR","CCND1"),]

peaks <- rownames(obj[["peaks"]])
peaks_gr <- StringToGRanges(regions = peaks)
#Visualization of scMultiomic of HCC
meta_all <- c()
CNmt_all <- data.frame(matrix(ncol = 0, nrow = length(peaks)))
for (i in 1:length(sampleIDs)){
    clt <- as.character(sampleIDs[i])
    print(clt)
    outdir_clt <- paste0(cnvdir,'/',clt,"/filtered_ana_binSize1")

    input_obj <- readRDS(paste0(outdir_clt,"/EiCNVs.obj"))
    cnvres <- readRDS(paste0(outdir_clt,"/final.CNVres.rds")) 

    cell_anno_new <- input_obj@cell_anno.filt
    cellMeta <- cell_anno_new[!cell_anno_new$subCluster %in%"reference",,drop=F]
    cellMeta  <- cellMeta[!cellMeta$clone_merged %in%"removed",,drop=F]
    cellMeta <- cellMeta[order(cellMeta$clone_merged,cellMeta$subCluster),,drop=F]
    cellinfo <- cnvres$cellinfo
    cellMeta <- left_join(cellMeta,cellinfo,by=c("row"="cellname"))

    cellMeta <- cellMeta %>%
      add_count(clone) %>%
      as.data.frame()
    colnames(cellMeta)[colnames(cellMeta)=="n"] <- "Ncells"
    rownames(cellMeta)<- cellMeta$row
    cellMeta$CellProportion <- cellMeta$Ncells/nrow(cellMeta)
    
    meta_all <- rbind(meta_all,cellMeta)


    cloneInfo_table <- unique(cellMeta[,c("clone","aneuploidy_score","Ncells","CellProportion")])
    #CNV
    CNmt <- do.call(cbind,lapply(sort(names(cnvres$clonalest)),function(cluster,peaks_gr){
      # clonal_res <- cnvres$clonalest[[cluster]]$input_BinSegRatio
      integerCNV <- cnvres$clonalest[[cluster]]$seg.dat
      # integerCNV <- integerCNV[,c("segName","relativeCN","integerCN")]
      # clonal_res <- clonal_res[,!grepl("relativeCN|integerCN",colnames(clonal_res))]
      # clonal_res <- left_join(clonal_res,integerCNV,by="segName")
      # rownames(clonal_res) <- clonal_res$binID
      integerCNV <- integerCNV[,c("chr","start","end","segName","integerCN")]
      CNmt_peak <- align_Grange2bin(as.data.frame(peaks_gr),integerCNV)

      return(CNmt_peak$integerCN)
      },peaks_gr))
    # rownames(CNmt) <- rownames(cnvres$clonalest[[1]]$input_BinSegRatio)
    colnames(CNmt) <- sort(names(cnvres$clonalest))
    rownames(CNmt) <- peaks

    if(as.character(colnames(CNmt)[1])=="0"){colnames(CNmt) <- as.numeric(colnames(CNmt))+1}


    # ##scale clonal size
    Ncol <- round(cloneInfo_table$CellProportion*nrow(cellMeta))
    new_CNmt <- do.call(cbind, lapply(1:ncol(CNmt), function(i) {
      matrix(rep(CNmt[, i], Ncol[i]), ncol = Ncol[i])
    }))
    rownames(new_CNmt) <- rownames(CNmt)
    colnames(new_CNmt) <- rownames(cellMeta)

    # new_CNmt_bed <- data.frame(chr = sapply(strsplit(rownames(new_CNmt), "_|:|-"), "[", 1),
    # start = as.numeric(sapply(strsplit(rownames(new_CNmt), "_|:|-"), "[", 2)), 
    # end = as.numeric(sapply(strsplit(rownames(new_CNmt), "_|:|-"), "[", 3)),
    # new_CNmt)
    # CNmt_peak <- align_Grange2bin(as.data.frame(peaks_gr),new_CNmt_bed)
    # rownames(CNmt_peak) <- CNmt_peak$binID
    # CNmt_peak <- CNmt_peak[,10:ncol(CNmt_peak)]

    CNmt_all <- cbind(CNmt_all,new_CNmt)

}

meta_ID <- unique(obj@meta.data[,c("sampleID","sampleID_new")])
meta_all <- left_join(meta_all,meta_ID,by="sampleID")

#
meta_all<- meta_all %>%
     group_by(sampleID_new) %>%
     mutate(clone = if (min(as.numeric(as.character(clone))) == 0) as.character(as.numeric(as.character(clone))+1) else clone)%>%
     ungroup() %>%
     as.data.frame()
 table(meta_all$clone)
meta_all <- meta_all[order(meta_all$sampleID_new,meta_all$clone),]
meta_all$clone_final <- paste(meta_all$sampleID_new,meta_all$clone,sep="_") 

meta_all$color_sample <- color_list$sampleID_new[match(meta_all$sampleID_new,names(color_list$sampleID_new))]
names(cols_Palette_clone) <- names(table(meta_all$clone))
meta_all$color_clone <- cols_Palette_clone[match(meta_all$clone,names(cols_Palette_clone))]
rownames(meta_all) <- meta_all$row

write.csv(meta_all,paste0(outdir,"/CellMeta_all.csv"),row.names = T)
meta_all$sampleID_new <- factor(meta_all$sampleID_new,levels= unique(meta_all$sampleID_new))
meta_all$clone_final <- factor(meta_all$clone_final,levels= unique(meta_all$clone_final))
meta_all <- droplevels(meta_all)

CNmt_all <- CNmt_all[,rownames(meta_all)]
saveRDS(CNmt_all,paste0(outdir,"/clonalCNV_peakMT_all.rds"))

table(meta_all$sampleID_new,meta_all$clone_final)
# meta_all$clone[meta_all$sampleID_new %in%"HCC1" &meta_all$clone=="3"] <-"1"



# meta_all <- read.csv(paste0(outdir,"/CellMeta_all.csv"),row.names = 1)
# meta_all$sampleID_new <- factor(meta_all$sampleID_new,levels= unique(meta_all$sampleID_new))
# meta_all$clone_final <- factor(meta_all$clone_final,levels= unique(meta_all$clone_final))

colOrder = seeds <- rownames(meta_all)
obj_sub <- subset(obj,cells =seeds)
if(length(obj_sub@neighbors)<1){
    DefaultAssay(obj_sub) <- "SCT"
    obj_sub <- obj_sub %>% 
        SCTransform(verbose = FALSE)%>%
        RunPCA(features = obj_sub@assays$SCT@var.features,verbose = FALSE)%>%
        FindNeighbors(dims = 1:30,force.recalc=T,return.neighbor=TRUE,k.param=100,assay = "SCT")
    DefaultAssay(obj_sub) <- "peaks"
    obj_sub <- FindNeighbors(object = obj_sub, force.recalc=T,reduction = 'lsi', return.neighbor=TRUE,dims = 2:30,assay = "peaks")

}

pseudoExp <- gen_pseudo_mat(obj_sub,seeds,Neighbors=50,assay="SCT",expMatrix=obj_sub[["SCT"]]@data)
pseudoPeak <- gen_pseudo_mat(obj_sub,seeds,Neighbors=50,assay="SCT",expMatrix=obj_sub[["peaks"]]@data)

saveRDS(pseudoExp,paste0(outdir,"/pseudoExp_50MetaCell.rds"))
saveRDS(pseudoPeak,paste0(outdir,"/pseudoPeak_50MetaCell.rds"))
write.csv(seeds,paste0(outdir,"/seeds.csv"))
#pseudoExp <- readRDS(paste0(outdir,"/pseudoExp_50MetaCell.rds"))
#pseudoPeak <- readRDS(paste0(outdir,"/pseudoPeak_50MetaCell.rds"))
#CNmt_all <- readRDS(paste0(outdir,"/clonalCNV_peakMT_all-1015.rds"))

#VariableGenes <- VariableFeatures(FindVariableFeatures(obj_sub, selection.method = "vst", nfeatures = round(0.75*dim(obj_sub)[1])))
# VariableGenes <- obj_sub[["SCT"]]@var.features
# peaks_in <- find_local_peaks(obj_sub,VariableGenes,annotations=annotations,peak.assa="peaks",distance=500000,CollapseTranscript=F) #Default distance is 500 kb
# peaks_in <- peaks_in[!(peaks_in$seqnames%in% c("chrX","chrY")),]
# #write.table(peaks_in,paste0(outdir_clt,"/totalLinks_500kb_VarGene.txt"),col.names = T,row.names = F,quote=F,sep="\t")

length(intersect(rownames(pseudoExp),unique(sigLink$gene)))
subexp=pseudoExp[rownames(pseudoExp) %in%unique(sigLink$gene),colOrder]
subpeak=pseudoPeak[unique(sigLink$peak),colOrder]

# system.time(
# sigLink <- gene_atac_cor.my(gene.mat=subexp,peak.mat=subpeak,links.df=peaks_in,pvalue_cutoff = 0.05,
#                             score_cutoff = 0.05)
# )
# colnames(sigLink)[colnames(sigLink)=="gene_name"] <- "gene"   

peak_exp <- as.matrix(subpeak)
colnames(peak_exp) <- gsub("\\.","-",colnames(peak_exp))
peak_exp <- peak_exp[rownames(peak_exp)%in%sigLink$peak,]
index1=match(sigLink$peak,row.names(peak_exp))
cell_indx1 <- match(colOrder,colnames(peak_exp))
peak_exp=as.matrix(peak_exp[index1,cell_indx1])
dim(peak_exp)
saveRDS(peak_exp,paste0(outdir,"/AllSample_50MetaCell_sigLink_ATAC.rds"))
# peak_exp <- readRDS(paste0(outdir,"/AllSample_50MetaCell_sigLink_ATAC.rds"))


rna_exp <- as.matrix(subexp)
colnames(rna_exp) <- gsub("\\.","-",colnames(rna_exp))
rna_exp <- rna_exp[rownames(rna_exp)%in%sigLink$gene,]
index2=match(sigLink$gene,row.names(rna_exp))
cell_indx2 <- match(colOrder,colnames(rna_exp))
rna_exp=as.matrix(rna_exp[index2,cell_indx2])
dim(rna_exp)
length(unique(rownames(rna_exp)))
saveRDS(rna_exp,paste0(outdir,"/AllSample_50MetaCell_sigLink_RNA.rds"))
# rna_exp <- readRDS(paste0(outdir,"/AllSample_50MetaCell_sigLink_RNA.rds"))
    

cnv_exp <- as.matrix(CNmt_all)  ##using CNV 
rownames(cnv_exp) <- gsub("\\_","-",rownames(cnv_exp))
colnames(cnv_exp) <- gsub("\\.","-",colnames(cnv_exp))
index3=match(sigLink$peak,rownames(cnv_exp))
cell_indx3 <- match(colOrder,colnames(cnv_exp))
cnv_exp=as.matrix(cnv_exp[index3,cell_indx3])
any(is.na(cnv_exp))
length(which(rowSums(cnv_exp)>0)) #number of non-zero peaks
dim(cnv_exp)
cnv_exp[1:10,1:5]
saveRDS(cnv_exp,paste0(outdir,"/AllSample_50MetaCell_sigLink_CNV.rds"))
# cnv_exp<- readRDS(paste0(outdir,"/AllSample_50MetaCell_sigLink_CNV.rds"))


mATAC <- .rowZscores(peak_exp,min = -4, max = 4,limit=T)
mRNA <- .rowZscores(rna_exp,min = -4, max = 4,limit=T)
mCNV <- cnv_exp

#row order
# library(parallelDist)
set.seed(123)
# d = dist(mRNA, method = "euclidean")
d = parallelDist::parDist(mRNA, method = "euclidean")
hclust_dist<- as.dist(d)
hclust_dist[is.na(hclust_dist)] <- 0
hclust_dist[is.nan(hclust_dist)] <- 0
sum(is.infinite(hclust_dist)) 

clust=hclust(hclust_dist,method = "ward.D2")
rowOrder <- clust$order
saveRDS(clust,paste0(outdir,"/hclust_rna.rds"))
# clust <- readRDS(paste0(outdir,"/hclust_rna.rds"))
# rowOrder <- clust$order


colOrder <- rownames(meta_all)
colData <- meta_all[,c("sampleID_new","clone_final"),drop=F]
colnames(colData) <- c("sample","clone")



meta_all_clone <- unique(meta_all[,c("sampleID_new","color_sample", "clone_final", "color_clone")])
left_anno_cols <- list()
left_anno_cols[['clone']] <- meta_all_clone$color_clone
names(left_anno_cols[['clone']]) <- meta_all_clone$clone_final
left_anno_cols[['sample']] <- unique(meta_all_clone$color_sample)
names(left_anno_cols[['sample']]) <- unique(meta_all_clone$sampleID_new)
################add Correlation annotation for RNA
####2. Correlation between clone-level gene expression and CN
cloneID <- unique(meta_all$clone_final)

metaCNV <- do.call(cbind,lapply(cloneID, function(clone,meta_all,cnv_exp){
  cellnames <- meta_all$row[meta_all$clone_final==clone]
  subdata <- cnv_exp[,cellnames]
  return(apply(subdata, 1, mean))
},meta_all,cnv_exp))

metaRNA <- do.call(cbind,lapply(cloneID, function(clone,meta_all,rna_exp){
  cellnames <- meta_all$row[meta_all$clone_final==clone]
  subdata <- rna_exp[,cellnames]
  return(apply(subdata, 1, mean))
},meta_all,rna_exp))
colnames(metaCNV) <- cloneID
colnames(metaRNA) <- cloneID

vargene <- unique(sigLink$gene)
geneindex <- match(vargene,row.names(metaRNA))
metaRNA1 <- metaRNA[geneindex,]
metaCNV1 <- metaCNV[geneindex,]

clone.cor <- do.call(rbind,lapply(1:dim(metaRNA1)[1], function(r,metaCNV1,metaRNA1,vargene){
  if (sum(!is.na(metaCNV1[r,]))>15){
    res2 <- cor.test(metaCNV1[r,],metaRNA1[r,])
    return(c(res2$estimate,res2$p.value,vargene[r]))
  }
},metaCNV1,metaRNA1,vargene))
clone.cor <- data.frame(gene = clone.cor[,3],RNA.r = clone.cor[,1],RNA.p = clone.cor[,2])

nRNA <- sum(!is.na(clone.cor$RNA.r))
index_CNVexpl <- which(clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05&!is.na(clone.cor$RNA.r))
genes_CNVexpl <- clone.cor$gene[index_CNVexpl]
n2 <- length(index_CNVexpl) 
RNA.CNV <- c(n2,nRNA-n2)
names(RNA.CNV) <- c(paste0("CNV explained: ",round(n2/nRNA*100,2),"%"),"Others")
index_NA <- which(is.na(clone.cor$RNA.r))
genes_NA <- clone.cor$gene[index_NA]

row_ano_right <- data.frame(row.name=rownames(rna_exp),Group="cCREs")
row_ano_right$Group[row_ano_right$row.name %in%genes_CNVexpl] <- "CNV explained"
row_ano_right$Group[row_ano_right$row.name %in%genes_NA] <- NA
left_anno_cols[['Group']] <- c("#fdcfa2", "#9ecae1")
names(left_anno_cols[['Group']])<- c("CNV explained","cCREs")

###combined heatmap    
pl <- c()
for(pp in c("ATAC","RNA","CNV")){
   heatmap_legend_param_ls = list(
                        direction = "horizontal",
                        #title_position = "leftcenter-rot",
                        legend_height = unit(3, "cm"))
   row_ha=NULL

    if(pp=="ATAC"){
        pal_col = paletteContinuous("solarExtra")
        mat <-mATAC 
        fileName <- paste0(pp," Z-Scores\n", nrow(mat), " Links, ",length(unique(rownames(mat)))," peaks")
    }
    if(pp=="RNA"){
        # pal_col = paletteContinuous("blueYellow")
        mat <-mRNA
        fileName <- paste0(pp," Z-Scores\n",nrow(mat), " Links, ", length(unique(rownames(mat))), " genes")
        
        brks = seq(-4, 4, length.out = 11)
        colorss <- colorRampPalette(c("#6ABDB2",'#94D5CB',"#BCE5DF","#DFF3F0",'#E4F4F2', "#FFFFFF","#F2E2B9","#E5CC91","#D5AF66","#bf812d","#8c510a"))(11) 
        pal_col <-colorRamp2(brks,colorss)

       row_ha = rowAnnotation(df = row_ano_right[rowOrder,'Group',drop=F],
        show_annotation_name = TRUE,
        col = left_anno_cols, 
        gp = gpar(col = "NA"))


    }
    if(pp=="CNV"){
        # pal_col = colorRampPalette(colors = c("darkblue", "white", "darkred"))(16)#length(bk)
        mat <-mCNV

        fileName <- paste0("Integer CN\n",nrow(mat), " Links, ", length(unique(rownames(mat))), " bins")
        # colData <- clone_info
        # colOrder <- rownames(clone_info)
        colorss <- colorRampPalette(c("#3d8bff","#C4D8F5", "grey95","#f9dcc4","#f79d65","#f27059","#85182a"))(7)
        pal_col <-colorRamp2(c(1,2,3,4,5,6,7),colorss)
        at_brk <- c(1:7)
        label_brk <- c(as.character(c(1:6)),"7+")
        heatmap_legend_param_ls <- list(
                  direction = "horizontal",
                  at = at_brk,
                  labels = label_brk,
                  legend_height = unit(3, "cm"), #图例长度
                  color_bar = "discrete")


    }
    ht1Anno <- HeatmapAnnotation(
        df = colData,
        col = left_anno_cols, 
        show_annotation_name = TRUE,
        gp = gpar(col = "NA"),
        annotation_legend_param =
        list(nrow = min(5, max(round(nrow(colData)/5), 1)))
    )

    ht_tumor <- Heatmap(mat[rowOrder,colOrder],
                        name = fileName,
                        show_column_names = F, 
                        show_row_names = F,
                        border = F,
                        row_title = NULL, #"%s",
                        row_title_gp = gpar(fontsize = 10),
                        heatmap_legend_param = heatmap_legend_param_ls,
                        col=pal_col,
                        cluster_columns = F,
                        cluster_rows = F,
                        #Annotation
                        top_annotation = ht1Anno,
                        right_annotation = row_ha)
     pl[[pp]] <- ht_tumor
}
pdf(paste0(outdir, "/Heatmap_P2GLinks_VarGene.scale-241031.pdf"),width = 12,height = 8)
ht_list = pl[['CNV']]+pl[['ATAC']]+pl[['RNA']]
draw(ht_list,merge_legend = TRUE,heatmap_legend_side="bottom",annotation_legend_side = "bottom")
dev.off()



###---------------------------------------###
###(2)plot heatmap basd on genomic regions
###---------------------------------------###
#（2-1） all samples
genes <- rownames(pseudoExp)
genes_gr <- annotations[annotations$gene_name %in% genes,]
# shrink to TSS position
tss.positions <- resize(genes_gr, width = 1, fix = 'start')
tss.positions <- tss.positions[!grepl("^chrM|^Mt|^MT|^chrX|^chrY|^X|^Y|^GL|^KI",as.character(seqnames(tss.positions))) ]
tss.positions <- Extend(
    x = tss.positions,
    upstream = 1000,
    downstream = 1000,
    from.midpoint = TRUE
  )

genes_gr <- as.data.frame(tss.positions) %>%
  group_by(gene_name) %>%
  dplyr::filter(gene_version == max(gene_version)) %>%
  ungroup()%>%as.data.frame()
genes_gr <- genes_gr[,c("seqnames","start","end","gene_name")]

peaks <- rownames(pseudoPeak)
peaks_gr <- StringToGRanges(regions = peaks)
peaks_gr$peaks <- peaks

peaks_gr_al <- align_Grange2bin(genes_gr,as.data.frame(peaks_gr))
peaks_gr_al <- na.omit(peaks_gr_al) ###Heatmap rows
head(genes_gr)
index2 <- peaks_gr_al$gene_name
colOrder <- rownames(meta_all)
cell_indx2 <- match(colOrder,colnames(pseudoExp))
rna_exp=as.matrix(pseudoExp[index2,cell_indx2])
# rna_exp <- readRDS(paste0(outdir,"/AllSample_50MetaCell_RNA.rds"))
saveRDS(rna_exp,paste0(outdir,"/AllSample_50MetaCell_RNA.rds"))
mRNA <- .rowZscores(rna_exp,min = -2, max = 2,limit=T)
dim(mRNA)
mRNA[1:3,1:3]


peak_exp <- pseudoPeak[peaks_gr_al$peaks,colOrder]
saveRDS(peak_exp,paste0(outdir,"/AllSample_50MetaCell_ATAC.rds"))
# peak_exp <- readRDS(paste0(outdir,"/AllSample_50MetaCell_ATAC.rds"))
mATAC <- .rowZscores(peak_exp,min = -3, max = 3,limit=T)
dim(mATAC)



cnv_exp <- as.matrix(CNmt_all)  ##using CNV 
rownames(cnv_exp) <- gsub("\\_","-",rownames(cnv_exp))
colnames(cnv_exp) <- gsub("\\.","-",colnames(cnv_exp))
cnv_exp <- cnv_exp[!grepl("chrX|chrY",rownames(cnv_exp)),,drop=F]
cell_indx3 <- match(colOrder,colnames(cnv_exp))
mCNV <- cnv_exp[peaks_gr_al$peaks,cell_indx3]
saveRDS(mCNV,paste0(outdir,"/AllSample_50MetaCell_CNV.rds"))
dim(mCNV)
mCNV[1:3,1:3]
#mCNV <- readRDS(paste0(outdir,"/AllSample_50MetaCell_CNV.rds"))


#heatmap

colData <- meta_all[,c("sampleID_new","clone_final"),drop=F]
colnames(colData) <- c("sample","clone")

meta_all_clone <- unique(meta_all[,c("sampleID_new", "clone_final", "color_clone","color_sample")])
left_anno_cols <- list()
left_anno_cols[['clone']] <- meta_all_clone$color_clone
names(left_anno_cols[['clone']]) <- meta_all_clone$clone_final
left_anno_cols[['sample']] <- unique(meta_all_clone$color_sample)
names(left_anno_cols[['sample']]) <- unique(meta_all_clone$sampleID_new)

pl2 <- c()
for(pp in c("ATAC","RNA","CNV")){
   heatmap_legend_param_ls = list(
                        direction = "horizontal",
                        #title_position = "leftcenter-rot",
                        legend_height = unit(3, "cm"))

   row_split <- peaks_gr_al$chr
    if(pp=="ATAC"){
        colorss = colorRampPalette(colorPalettes[["solarExtra"]])(5)     
        pal_col <-colorRamp2(c(-3,-1,0,1,3),colorss)
        mat <-mATAC[,colOrder] 
        fileName <- paste0(pp," Z-Scores\n")

    }
    if(pp=="RNA"){
        mat <-mRNA[,colOrder] 
        fileName <- paste0(pp," Z-Scores\n",length(unique(peaks_gr_al$gene_name)), " genes")
        brks = seq(-2, 2, length.out = 11)
        colorss <- colorRampPalette(c("#6ABDB2",'#94D5CB',"#BCE5DF","#DFF3F0",'#E4F4F2', "#FFFFFF","#F2E2B9","#E5CC91","#D5AF66","#bf812d","#8c510a"))(11) 
        pal_col <-colorRamp2(brks,colorss)     
    }
    if(pp=="CNV"){
        mat <-mCNV[,colOrder]

        fileName <- paste0("Integer CN\n")
        colorss <- colorRampPalette(c("#3d8bff","#C4D8F5", "grey95","#f9dcc4","#f79d65","#f27059","#85182a"))(7)
        pal_col <-colorRamp2(c(1,2,3,4,5,6,7),colorss)
        at_brk <- c(1:7)
        label_brk <- c(as.character(c(1:6)),"7+")
        heatmap_legend_param_ls <- list(
                  direction = "horizontal",
                  at = at_brk,
                  labels = label_brk,
                  legend_height = unit(3, "cm"), #图例长度
                  color_bar = "discrete")
    }
    ht1Anno <- HeatmapAnnotation(
        df = colData,
        col = left_anno_cols, 
        show_annotation_name = TRUE,
        gp = gpar(col = "NA"),
        annotation_legend_param =
        list(nrow = min(5, max(round(nrow(colData)/5), 1)))
    )

    

    ht_tumor <- Heatmap(mat,
                        name = fileName,
                        show_column_names = F, 
                        show_row_names = F,
                        border = F,
                        row_title = NULL, #"%s",
                        row_title_gp = gpar(fontsize = 10),
                        heatmap_legend_param = heatmap_legend_param_ls,
                        col=pal_col,
                        cluster_columns = F,
                        cluster_rows = F,
                        #Annotation
                        row_split=row_split,
                        row_gap = unit(0.5, "mm"),
                        top_annotation = ht1Anno)
     pl2[[pp]] <- ht_tumor
}
pdf(paste0(outdir, "/Heatmap_GenomicPositionTSS.pdf"),width = 12,height = 8)
ht_list = pl2[["CNV"]]+pl2[["ATAC"]]+pl2[["RNA"]]
draw(ht_list,merge_legend = TRUE,heatmap_legend_side="bottom",annotation_legend_side = "bottom")
dev.off()

#save clone-level data
cellMeta <- data.frame(cells=rownames(meta_all),clone=meta_all[,c("clone_final")])
rownames(cellMeta) <- rownames(meta_all)

cnv_exp_b <- pseudo_bulk_v2(mCNV[,rownames(cellMeta)],group_ano = cellMeta,method ="median",adjust_zero=TRUE)
peak_exp_b <- pseudo_bulk_v2(peak_exp[,rownames(cellMeta)],group_ano = cellMeta,method ="mean",adjust_zero=TRUE)
rna_exp_b <- pseudo_bulk_v2(rna_exp[,rownames(cellMeta)],group_ano = cellMeta,method ="mean",adjust_zero=TRUE)
saveRDS(cnv_exp_b,paste0(outdir,"/AllSample_50MetaCell_CNA_clonalMean.rds"))
saveRDS(peak_exp_b,paste0(outdir,"/AllSample_50MetaCell_ATAC_clonalMean.rds"))
saveRDS(rna_exp_b,paste0(outdir,"/AllSample_50MetaCell_RNA_clonalMean.rds"))







#(2-2) for each sample
gplinks <- read.table(paste0(outdir,"/Links_500kb_VarGene.txt"),header=T,sep="\t")
sigLink <- gplinks[gplinks$pvalue<0.05&gplinks$score>0.05,,drop=F]
sigLink <- sigLink[order(sigLink$pvalue,sigLink$score,decreasing = c(FALSE,TRUE)),]
sigLink <- sigLink[!(sigLink$seqnames%in% c("chrX","chrY")),]

plot4tss=FALSE
plot4VarGene = FALSE
wholeGenome = TRUE

for (i in 1:length(sampleIDs)){
    clt <- as.character(sampleIDs[i])
    print(clt)
    outdir_clt <- paste0(cnvdir,'/',clt,"/filtered_ana_binSize1")
    input_obj <- readRDS(paste0(outdir_clt,"/EiCNVs.obj"))
    cnvres <- readRDS(paste0(outdir_clt,"/final.CNVres.rds")) 

    cell_anno_new <- input_obj@cell_anno.filt
    cellMeta <- cell_anno_new[!cell_anno_new$subCluster %in%"reference",,drop=F]
    cellMeta  <- cellMeta[!cellMeta$clone_merged %in%"removed",,drop=F]
    cellMeta <- cellMeta[order(cellMeta$clone_merged,cellMeta$subCluster),,drop=F]
    colnames(cellMeta) <- gsub("clone_merged","clone",colnames(cellMeta))
    cellMeta <- cellMeta %>%
      add_count(clone) %>%
      as.data.frame()
    colnames(cellMeta)[colnames(cellMeta)=="n"] <- "Ncells"
    rownames(cellMeta)<- cellMeta$row
    cellMeta$CellProportion <- cellMeta$Ncells/nrow(cellMeta)
    if("0" %in% as.character(cellMeta$clone)){cellMeta$clone <- as.numeric(as.character(cellMeta$clone))+1}


    obj_clt <- subset(obj,cells = rownames(cellMeta))

    peaks_clt <- rownames(obj_clt[["peaks"]])
    peaks_clt_gr <- StringToGRanges(regions = peaks_clt)

    cloneInfo_table <- unique(cellMeta[,c("clone","aneuploidy_score","Ncells","CellProportion")])
    #CNV
    CNmt <- do.call(cbind,lapply(sort(names(cnvres$clonalest)),function(cluster){
      clonal_res <- cnvres$clonalest[[cluster]]$input_BinSegRatio
      integerCNV <- cnvres$clonalest[[cluster]]$seg.dat
      integerCNV <- integerCNV[,c("segName","relativeCN","integerCN")]
      clonal_res <- clonal_res[,!grepl("relativeCN|integerCN",colnames(clonal_res))]
      clonal_res <- left_join(clonal_res,integerCNV,by="segName")
      rownames(clonal_res) <- clonal_res$binID
      return(clonal_res$integerCN)
      }))
    rownames(CNmt) <- rownames(cnvres$clonalest[[1]]$input_BinSegRatio)
    colnames(CNmt) <- sort(names(cnvres$clonalest))
    if(as.character(colnames(CNmt)[1])=="0"){colnames(CNmt) <- as.numeric(colnames(CNmt))+1}


    # ##scale clonal size
    Ncol <- round(cloneInfo_table$CellProportion*nrow(cellMeta))
    new_CNmt <- do.call(cbind, lapply(1:ncol(CNmt), function(i) {
      matrix(rep(CNmt[, i], Ncol[i]), ncol = Ncol[i])
    }))
    rownames(new_CNmt) <- rownames(CNmt)
    colnames(new_CNmt) <- rownames(cellMeta)
    if(plot4tss){
        new_CNmt_bed <- data.frame(chr = sapply(strsplit(rownames(new_CNmt), "_|:|-"), "[", 1),
            start = as.numeric(sapply(strsplit(rownames(new_CNmt), "_|:|-"), "[", 2)), 
            end = as.numeric(sapply(strsplit(rownames(new_CNmt), "_|:|-"), "[", 3)),
            new_CNmt)
        CNmt_peak <- align_Grange2bin(as.data.frame(peaks_clt_gr),new_CNmt_bed)
        rownames(CNmt_peak) <- CNmt_peak$binID
        CNmt_peak <- CNmt_peak[,10:ncol(CNmt_peak)] 
    }else if(wholeGenome){
        CNmt_peak <- new_CNmt
    }



    cnv_exp <- as.matrix(CNmt_peak)  ##using CNV 
    rownames(cnv_exp) <- gsub("\\_","-",rownames(cnv_exp))
    colnames(cnv_exp) <- gsub("\\.","-",colnames(cnv_exp))
    cnv_exp <- cnv_exp[!grepl("chrX|chrY",rownames(cnv_exp)),,drop=F]
    dim(cnv_exp)
    cnv_exp[1:3,1:3]

    colOrder <- rownames(cellMeta)
    bins_gr <- StringToGRanges(regions = rownames(cnv_exp))



    #RNA/ATAC profile
    if(length(obj_clt@neighbors)<1){
        DefaultAssay(obj_clt) <- "SCT"
        obj_clt <- obj_clt %>% 
            SCTransform(verbose = FALSE)%>%
            RunPCA(features = obj_clt@assays$SCT@var.features,verbose = FALSE)%>%
            FindNeighbors(dims = 1:30,force.recalc=T,return.neighbor=TRUE,k.param=100,assay = "SCT")
        DefaultAssay(obj_clt) <- "peaks"
        obj_clt <- FindNeighbors(object = obj_clt, force.recalc=T,reduction = 'lsi', return.neighbor=TRUE,dims = 2:30,assay = "peaks")

    }

    pseudoExp <- gen_pseudo_mat(obj_clt,rownames(cellMeta),Neighbors=50,assay="SCT",expMatrix=obj_clt[["SCT"]]@data)
    pseudoPeak <- gen_pseudo_mat(obj_clt,rownames(cellMeta),Neighbors=50,assay="SCT",expMatrix=obj_clt[["peaks"]]@data)


    genes <- rownames(pseudoExp)
    genes_gr <- annotations[annotations$gene_name %in% genes,]

    if(plot4tss){
        # shrink to TSS position
        tss.positions <- resize(genes_gr, width = 1, fix = 'start')
        tss.positions <- tss.positions[!grepl("^chrM|^Mt|^MT|^chrX|^chrY|^X|^Y|^GL|^KI",as.character(seqnames(tss.positions))) ]
        
        tss.positions <- Extend(
            x = tss.positions,
            upstream = 1000,
            downstream = 1000,
            from.midpoint = TRUE
          )

        genes_gr <- as.data.frame(tss.positions) %>%
          group_by(gene_name) %>%
          dplyr::filter(gene_version == max(gene_version)) %>%
          ungroup()%>%as.data.frame()
    }





    genes_gr <- as.data.frame(genes_gr)
    genes_gr <- genes_gr[,c("seqnames","start","end","gene_name")]



    if(plot4tss){
        peaks <- rownames(pseudoPeak)
        peaks_gr <- StringToGRanges(regions = peaks)
        peaks_gr$peaks <- peaks
        peaks_gr_al <- align_Grange2bin(genes_gr,as.data.frame(peaks_gr))
        peaks_gr_al <- na.omit(peaks_gr_al)
        index2 <- peaks_gr_al$gene_name

        peak_exp <- pseudoPeak[peaks_gr_al$peaks,]

        cell_indx2 <- match(colOrder,colnames(pseudoExp))
        rna_exp=as.matrix(pseudoExp[index2,cell_indx2])
        count_lim <- round(quantile(rna_exp,0.95,na.rm=T));
        rna_exp[rna_exp>count_lim] <- count_lim
    }else if(plot4VarGene){
        length(intersect(rownames(pseudoExp),unique(sigLink$gene)))
        subexp=pseudoExp[rownames(pseudoExp) %in%unique(sigLink$gene),]
        subpeak=pseudoPeak[unique(sigLink$peak),]

        peak_exp <- as.matrix(subpeak)
        colnames(peak_exp) <- gsub("\\.","-",colnames(peak_exp))
        peak_exp <- peak_exp[rownames(peak_exp)%in%sigLink$peak,]
        index1=match(sigLink$peak,row.names(peak_exp))
        cell_indx1 <- match(colOrder,colnames(peak_exp))
        peak_exp=as.matrix(peak_exp[index1,cell_indx1])
        dim(peak_exp)

        rna_exp <- as.matrix(subexp)
        colnames(rna_exp) <- gsub("\\.","-",colnames(rna_exp))
        rna_exp <- rna_exp[rownames(rna_exp)%in%sigLink$gene,]
        index2=match(sigLink$gene,row.names(rna_exp))
        cell_indx2 <- match(colOrder,colnames(rna_exp))
        rna_exp=as.matrix(rna_exp[index2,cell_indx2])
        dim(rna_exp)
        length(unique(rownames(rna_exp)))

    }else if(wholeGenome){
        genes_peaks <- align_Grange2bin(as.data.frame(bins_gr),genes_gr)
        

        peak_exp <- pseudoPeak[rownames(cnv_exp),colOrder]
        dim(peak_exp)

        index2=match(genes_peaks$gene_name,row.names(pseudoExp))  
        cell_indx2 <- match(colOrder,colnames(pseudoExp))
        rna_exp=as.matrix(pseudoExp[index2,cell_indx2])
        # count_lim <- round(quantile(rna_exp,0.95,na.rm=T));
        # rna_exp[rna_exp>count_lim] <- count_lim
    }

    # count_lim_ATAC <- round(quantile(peak_exp,0.95,na.rm=T));
    # peak_exp[peak_exp>count_lim_ATAC] <- count_lim_ATAC
    if(plot4tss){
        mATAC <- smooth_and_denoise(peak_exp,window=10)
        mCNV <- cnv_exp[peaks_gr_al$peaks,]
    }else{
        mATAC <- peak_exp
        mCNV <- cnv_exp
    }

    mATAC <- .rowZscores(mATAC,min = -2, max = 2,limit=T)
    mRNA <- .rowZscores(rna_exp,min = -2, max = 2,limit=T)
    dim(mRNA)
    mRNA[1:3,1:3]
    dim(mATAC)
    dim(mCNV)
    

    #heatmap
    colData <- cellMeta[,c("clone"),drop=F]
    colnames(colData) <- c("clone")

    left_anno_cols <- list()
    color_r2 <- cols_Palette[1:length(unique(cellMeta$clone))]
    names(color_r2) <- sort(unique(cellMeta$clone))
    left_anno_cols[['clone']] <- color_r2


    pl2 <- c()
    for(pp in c("ATAC","RNA","CNV")){
       heatmap_legend_param_ls = list(
                            direction = "horizontal",
                            #title_position = "leftcenter-rot",
                            legend_height = unit(3, "cm"))
       row_split <- genes_peaks$chr

        if(pp=="ATAC"){
            colorss = colorRampPalette(colorPalettes[["solarExtra"]])(5)  
            pal_col <-colorRamp2(c(-4,-2,0,2,4),colorss)
            mat <-mATAC[,colOrder] 
            fileName <- paste0(pp," Z-Scores\n")
        }
        if(pp=="RNA"){
            mat <-mRNA[,colOrder] 
            fileName <- paste0(pp," Z-Scores\n",length(na.omit(unique(rownames(mat)))), " genes")
            # brks <- quantile(mat,c(0.05,seq(0.125,0.875,by=0.125),0.95),na.rm=TRUE)
            # colorss <- colorRampPalette(c("#0067C3", "#48D3E7","#AFFBF9", "#FFFFFF", "#FFFF00","#B37B02", "#6B1D03"))(9) #6: van Gogh
            # colorss = colorRampPalette(colorPalettes[["blueYellow"]])(9) 
            # pal_col <-colorRamp2(brks,colorss)
            # pal_col <-colorRamp2(c(-2, 0, 2), c("#23A5B4", "white", "#F8FA0D"))
            brks = seq(-2, 2, length.out = 11)
            colorss <- colorRampPalette(c("#6ABDB2",'#94D5CB',"#BCE5DF","#DFF3F0",'#E4F4F2', "#FFFFFF","#F2E2B9","#E5CC91","#D5AF66","#bf812d","#8c510a"))(11) 
            pal_col <-colorRamp2(brks,colorss)

        }
        if(pp=="CNV"){
            # pal_col = colorRampPalette(colors = c("darkblue", "white", "darkred"))(16)#length(bk)
            mat <-mCNV[,colOrder]

            fileName <- paste0("Integer CN")
            # colData <- clone_info
            # colOrder <- rownames(clone_info)
            colorss <- colorRampPalette(c("#3d8bff","#C4D8F5", "grey95","#f9dcc4","#f79d65","#f27059","#85182a"))(7)
            pal_col <-colorRamp2(c(1,2,3,4,5,6,7),colorss)
            at_brk <- c(1:7)
            label_brk <- c(as.character(c(1:6)),"7+")
            heatmap_legend_param_ls <- list(
                      direction = "horizontal",
                      at = at_brk,
                      labels = label_brk,
                      legend_height = unit(3, "cm"), #图例长度
                      color_bar = "discrete")
        }
        ht1Anno <- HeatmapAnnotation(
            df = colData,
            col = left_anno_cols, 
            show_annotation_name = TRUE,
            gp = gpar(col = "NA"),
            annotation_legend_param =
            list(nrow = min(5, max(round(nrow(colData)/5), 1)))
        )

        ht_tumor <- Heatmap(mat,
                            name = fileName,
                            show_column_names = F, 
                            show_row_names = F,
                            border = F,
                            row_title = NULL, #"%s",
                            row_title_gp = gpar(fontsize = 10),
                            heatmap_legend_param = heatmap_legend_param_ls,
                            col=pal_col,
                            cluster_columns = F,
                            cluster_rows = F,
                            #Annotation
                            row_split=row_split,
                            row_gap = unit(0.5, "mm"),
                            top_annotation = ht1Anno)

         pl2[[pp]] <- ht_tumor
    }
    if(plot4tss){
        outfile = paste0(outdir_clt, "/Heatmap_GenomicPositionTSS_",clt,".pdf")
    }else if(wholeGenome){
        outfile = paste0(outdir_clt, "/Heatmap_GenomicPosition_",clt,".pdf")
    }
    pdf(outfile,width = 12,height = 8)
    ht_list = pl2[[3]]+pl2[[1]]+pl2[[2]]
    draw(ht_list,merge_legend = TRUE,heatmap_legend_side="bottom",annotation_legend_side = "bottom")
    dev.off()

}


