setwd('~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/')

library(Seurat)
library(ggplot2)
library(ggsci)
source('./Code/Function/VlnPlotSelf.R')
cb_pattern = c('#a9011b', '#e4a826', '#d2dcad', '#dcd6b2','#94bea7', '#4e7989', '#75a0ae', '#80944e', '#547db1', '#9055a2', '#d43b36', '#f4c28f')
 
obj <- readRDS('·/HCC_MultiOmic_12sam_clean.rds')
 

## Figure4A
DimPlot(obj,group.by =   'sampleID_new',cols = cb_pattern,label = T)

## 

 
## Figure4-S1A
Idents(obj) = obj$sampleID_new
P1 = VlnPlotSelf(obj,features ='atac_fragments',lim = c(0,100000))
P2 = VlnPlotSelf(obj,features = c('atac_peak_region_fragments'))
P3 = VlnPlotSelf(obj,features = c('TSS.enrichment'))
P4 = VlnPlotSelf(obj,features = c('nFeature_RNA'))
P5 = VlnPlotSelf(obj,features = c('nCount_RNA'))
P = P1+P2+P3+P4+P5 
P
ggsave('../Figures/raw/Figure4//HCC_qc.pdf',plot = P,width = 22,height =8)


## Figure4-SB
DimPlot(obj,group.by = c('tissue','CellType.final'))


## Figure4-S1C
FeaturePlot(obj,features = c('PTPRC','CD3D','CD3E','CD68','PECAM1','VWF','ACTA2','COL1A1','ALB','HNF4A'))



 
