
#library(tsne)
#library(plotly)
library(umap) 
library(HelloRanges)
library(reshape)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggseqlogo)
library(chromVAR)

library(EnsDb.Hsapiens.v86)
#library(ggplot2)
#library(ggsci)
#library(ggpubr)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
library(Signac)
library(dplyr)
library(phangorn)
library("ggtree")
library(ggplot2)
library(tidytree)
library("TreeTools")
library(Seurat)
library(msigdbr)
library(data.table)
library(RColorBrewer)
library(Signac)

WNNcolor <-brewer.pal(12, "Set3") 
work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
datapath <- paste0(work_dir,'/TeaCNV_HCC_MultiOmic')
setwd(datapath)
samplelist <- list.files()
samplelist <- samplelist[grep("JS",samplelist)]
respath <- paste0(work_dir,"/clonal_complexity_res/links")
ifelse(!dir.exists(file.path(respath)), dir.create(file.path(respath)), FALSE)


#################2. clone水平基因表达与CN的相关性
datapath2 <- paste0(work_dir,"/Figures")
setwd(datapath2)
CellMeta_all <- read.csv("CellMeta_all.txt",sep="\t")
CNVdata <- readRDS("AllSample_50MetaCell_CNV.rds")
RNAdata <- readRDS("AllSample_50MetaCell_RNA.rds")
ATACdata <- readRDS("AllSample_50MetaCell_ATAC.rds")
cloneID <- unique(CellMeta_all$clone_final)
metaCNV <- do.call(cbind,lapply(cloneID, function(clone,CellMeta_all,CNVdata){
  cellnames <- CellMeta_all$row[CellMeta_all$clone_final==clone]
  subdata <- CNVdata[,cellnames]
  return(apply(subdata, 1, median))
},CellMeta_all,CNVdata))

metaATAC <- do.call(cbind,lapply(cloneID, function(clone,CellMeta_all,ATACdata){
  cellnames <- CellMeta_all$row[CellMeta_all$clone_final==clone]
  subdata <- ATACdata[,cellnames]
  return(apply(subdata, 1, median))
},CellMeta_all,ATACdata))

metaRNA <- do.call(cbind,lapply(cloneID, function(clone,CellMeta_all,RNAdata){
  cellnames <- CellMeta_all$row[CellMeta_all$clone_final==clone]
  subdata <- RNAdata[,cellnames]
  return(apply(subdata, 1, median))
},CellMeta_all,RNAdata))

colnames(metaCNV) <- cloneID
colnames(metaATAC) <- cloneID
colnames(metaRNA) <- cloneID

clone.cor <- do.call(rbind,lapply(1:dim(metaCNV)[1], function(r,metaCNV,metaRNA,metaATAC){
  if (sum(!is.na(metaCNV[r,]))>15){
    res1 <- cor.test(metaCNV[r,],metaATAC[r,])
    res2 <- cor.test(metaCNV[r,],metaRNA[r,])
    return(c(res1$estimate,res1$p.value,res2$estimate,res2$p.value))
  }
},metaCNV,metaRNA,metaATAC))
clone.cor <- data.frame(ATAC.r = clone.cor[,1],ATAC.p = clone.cor[,2],RNA.r = clone.cor[,3],RNA.p = clone.cor[,4])
nATAC <- sum(!is.na(clone.cor$ATAC.r))
n1 <- length(clone.cor$ATAC.r[clone.cor$ATAC.r > 0.3 & clone.cor$ATAC.p < 0.05&!is.na(clone.cor$ATAC.r)]) 
nRNA <- sum(!is.na(clone.cor$RNA.r))
n2 <- length(clone.cor$RNA.r[clone.cor$RNA.r > 0.3 & clone.cor$RNA.p < 0.05&!is.na(clone.cor$RNA.r)]) 
ATAC.CNV <- c(n1,nATAC-n1)
names(ATAC.CNV) <- c(paste0("CNV explained: ",round(n1/nATAC*100,2),"%"),"Others")
RNA.CNV <- c(n2,nRNA-n2)
names(RNA.CNV) <- c(paste0("CNV explained: ",round(n2/nRNA*100,2),"%"),"Others")
# "#3399CC"
# "#FF9933"
pdf(file=paste0(respath,"/clone.peak.RNA.CNV.density.pdf"),width = 4,height = 4)
plot(density(clone.cor$ATAC.r[!is.na(clone.cor$ATAC.r)],width=0.1),xlim=c(-1,1),ylim=c(0,1.4),col="#FF9933",lwd=2,main="",ylab="Density",xlab="Correlation")
par(new=T)
plot(density(clone.cor$RNA.r[!is.na(clone.cor$RNA.r)],width=0.1),xlim=c(-1,1),ylim=c(0,1.4),axes=F,col="#3399CC",lwd=2,main="",ylab="",xlab="")
abline(v=0.3,lty=2,col="grey")
legend("topleft",legend = c("ATAC","RNA"),col=c("#FF9933","#3399CC"),lwd=2,bty="n")
dev.off()
pdf(file=paste0(respath,"/clone.peak.RNA.CNV.explained.pie.pdf"),width = 8,height = 4)
par(mfrow=c(1,2),mar=c(1,6,2,6))
pie(ATAC.CNV,col=WNNcolor,main="ATAC")
pie(RNA.CNV,col=WNNcolor,main="RNA")
dev.off()

############4. 克隆水平variable gene表达与CN的相关性
setwd(datapath2)
vargenelist <- read.csv("sigLinks_500kb_VarGene.txt",sep="\t")
CellMeta_all <- read.csv("CellMeta_all.txt",sep="\t")
CNVdata <- readRDS("AllSample_50MetaCell_sigLink_CNV.rds")
RNAdata <- readRDS("AllSample_50MetaCell_sigLink_RNA.rds")
ATACdata <- readRDS("AllSample_50MetaCell_sigLink_ATAC.rds")
cloneID <- unique(CellMeta_all$clone_final)
metaCNV <- do.call(cbind,lapply(cloneID, function(clone,CellMeta_all,CNVdata){
  cellnames <- CellMeta_all$row[CellMeta_all$clone_final==clone]
  subdata <- CNVdata[,cellnames]
  return(apply(subdata, 1, mean))
},CellMeta_all,CNVdata))

metaATAC <- do.call(cbind,lapply(cloneID, function(clone,CellMeta_all,ATACdata){
  cellnames <- CellMeta_all$row[CellMeta_all$clone_final==clone]
  subdata <- ATACdata[,cellnames]
  return(apply(subdata, 1, mean))
},CellMeta_all,ATACdata))

metaRNA <- do.call(cbind,lapply(cloneID, function(clone,CellMeta_all,RNAdata){
  cellnames <- CellMeta_all$row[CellMeta_all$clone_final==clone]
  subdata <- RNAdata[,cellnames]
  return(apply(subdata, 1, mean))
},CellMeta_all,RNAdata))

colnames(metaCNV) <- cloneID
colnames(metaATAC) <- cloneID
colnames(metaRNA) <- cloneID


vargene <- unique(vargenelist$gene)
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
write.table(clone.cor,paste0(respath,"/vargene.CNV.cor.csv"),sep=",",col.names = T,row.names = F,quote = F)


vargene_all <- vargenelist$gene
clone.cor_allVar <- do.call(rbind,lapply(1:dim(metaCNV)[1], function(r,metaCNV,metaRNA,metaATAC,vargene_all){
  if (sum(!is.na(metaCNV[r,]))>15 &sum(!is.na(metaATAC[r,]))>15 & sum(!is.na(metaRNA[r,]))>15){
    res1 <- cor.test(metaCNV[r,],metaATAC[r,])
    res2 <- cor.test(metaCNV[r,],metaRNA[r,])
    return(c(res1$estimate,res1$p.value,res2$estimate,res2$p.value,vargene_all[r]))
  }
},metaCNV,metaRNA,metaATAC,vargene_all))


#plot pie chart
clone.cor_allVar <- data.frame(ATAC.r = clone.cor_allVar[,1],ATAC.p = clone.cor_allVar[,2],RNA.r = clone.cor_allVar[,3],RNA.p = clone.cor_allVar[,4])
nATAC <- sum(!is.na(clone.cor_allVar$ATAC.r))
n1 <- length(clone.cor_allVar$ATAC.r[clone.cor_allVar$ATAC.r > 0.3 & clone.cor_allVar$ATAC.p < 0.05&!is.na(clone.cor_allVar$ATAC.r)]) 
nRNA <- sum(!is.na(clone.cor_allVar$RNA.r))
n2 <- length(clone.cor_allVar$RNA.r[clone.cor_allVar$RNA.r > 0.3 & clone.cor_allVar$RNA.p < 0.05&!is.na(clone.cor_allVar$RNA.r)]) 
ATAC.CNV_allVar <- c(n1,nATAC-n1)
names(ATAC.CNV_allVar) <- c(paste0("CNV explained: ",round(n1/nATAC*100,2),"%"),"Others")
RNA.CNV_allVar <- c(n2,nRNA-n2)
names(RNA.CNV_allVar) <- c(paste0("CNV explained: ",round(n2/nRNA*100,2),"%"),"Others")

pdf(file=paste0(respath,"/clone.peak.RNA.CNV.explained.pie_VarGenes.pdf"),width = 8,height = 4)
par(mfrow=c(1,2),mar=c(1,6,2,6))
pie(ATAC.CNV_allVar,col=WNNcolor,main="ATAC")
pie(RNA.CNV_allVar,col=WNNcolor,main="RNA")
dev.off()




#####5. CNV和peak解释的variable gene富集的pathway
library(fgsea)
library(tidyr)
set.seed(42)
inpGS <- data.table(msigdbr(species = "Homo sapiens", category = "H"))
inpGS$gs_name <- gsub("HALLMARK_", "", inpGS$gs_name)
inpGS$gs_name <- gsub("_SIGNALING", "", inpGS$gs_name)
inpGS <- split(inpGS$gene_symbol, inpGS$gs_name)


MP <- read.csv(paste0(respath,"/TableS2_meta.program.csv")) #Table S2 from https://doi.org/10.1038/s41586-023-06130-4
MPname <- colnames(MP)
MPname <- unlist(lapply(MPname, function(ele){
  ele <- unlist(strsplit(ele,split="[.]"))
  ele <- ele[ele!=""]
  ele <- ele[-grep("MP",ele)]
  ele <- paste(ele,collapse = " ")
  return(ele)
}))
colnames(MP) <- MPname

MP <- as.list(MP)

library(hypeR)
hyp_obj.hall <- hypeR(inpGS, list(gene = clone.cor$gene))
res <- hyp_obj.hall$data
hallmarkhyp <- do.call(rbind,lapply(1:length(res), function(i,res){
  pathway <- names(res)[i]
  output <- c(pathway,res[[i]]$data$pval,res[[i]]$data$fdr,res[[i]]$data$overlap,res[[i]]$data$signature)
  return(output)
},res))
hallmarkhyp <- data.frame(pathway=hallmarkhyp[,1],pval = as.numeric(hallmarkhyp[,2]),fdr = as.numeric(hallmarkhyp[,3]),overlap=as.numeric(hallmarkhyp[,4]),signature =as.numeric(hallmarkhyp[,5]))

hyp_obj <- hypeR(inpGS, list(gene = clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05]))
res <- hyp_obj$data
hallmarkhyp1 <- do.call(rbind,lapply(1:length(res), function(i,res){
  pathway <- names(res)[i]
  output <- c(pathway,res[[i]]$data$pval,res[[i]]$data$fdr,res[[i]]$data$overlap,res[[i]]$data$signature)
  return(output)
},res))
hallmarkhyp.cnv <- data.frame(pathway=hallmarkhyp1[,1],pval = as.numeric(hallmarkhyp1[,2]),fdr = as.numeric(hallmarkhyp1[,3]),overlap=as.numeric(hallmarkhyp1[,4]),signature =as.numeric(hallmarkhyp1[,5]))

hyp_obj <- hypeR(inpGS, list(gene = setdiff(clone.cor$gene,clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05])))
res <- hyp_obj$data
hallmarkhyp1 <- do.call(rbind,lapply(1:length(res), function(i,res){
  pathway <- names(res)[i]
  output <- c(pathway,res[[i]]$data$pval,res[[i]]$data$fdr,res[[i]]$data$overlap,res[[i]]$data$signature)
  return(output)
},res))
hallmarkhyp.peak <- data.frame(pathway=hallmarkhyp1[,1],pval = as.numeric(hallmarkhyp1[,2]),fdr = as.numeric(hallmarkhyp1[,3]),overlap=as.numeric(hallmarkhyp1[,4]),signature =as.numeric(hallmarkhyp1[,5]))

hallmarkhyp$group <- "All"
hallmarkhyp$size <- dim(clone.cor)[1]
hallmarkhyp.cnv$group <- "CNV"
hallmarkhyp.cnv$size <- length(clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05])
hallmarkhyp.peak$group <- "peak"
hallmarkhyp.peak$size <- dim(clone.cor)[1]-length(clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05])
hallmark.gsea <- rbind(hallmarkhyp,hallmarkhyp.cnv,hallmarkhyp.peak)
hallmark.gsea$ratio <- hallmark.gsea$overlap/hallmark.gsea$size
hallmark.gsea1 <- hallmark.gsea[hallmark.gsea$group!="All",]
hallmark.gsea1 <- hallmark.gsea1[hallmark.gsea1$pval < 0.01,]
hallmark.gsea1 <- hallmark.gsea1[order(hallmark.gsea1$group,hallmark.gsea1$ratio),]
pathway <- unique(hallmark.gsea1$pathway)
index <- match(pathway,hallmark.gsea1$pathway)
hallmark.gsea1 <- hallmark.gsea1[index,]
hallmark.gsea1$pathway <- factor(hallmark.gsea1$pathway,levels = pathway)

hyp_obj.mp <- hypeR(MP, list(gene = clone.cor$gene))
res <- hyp_obj.mp$data
hallmarkhyp <- do.call(rbind,lapply(1:length(res), function(i,res){
  pathway <- names(res)[i]
  output <- c(pathway,res[[i]]$data$pval,res[[i]]$data$fdr,res[[i]]$data$overlap,res[[i]]$data$signature)
  return(output)
},res))
hallmarkhyp <- data.frame(pathway=hallmarkhyp[,1],pval = as.numeric(hallmarkhyp[,2]),fdr = as.numeric(hallmarkhyp[,3]),overlap=as.numeric(hallmarkhyp[,4]),signature =as.numeric(hallmarkhyp[,5]))

hyp_obj <- hypeR(MP, list(gene = clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05]))
res <- hyp_obj$data
hallmarkhyp1 <- do.call(rbind,lapply(1:length(res), function(i,res){
  pathway <- names(res)[i]
  output <- c(pathway,res[[i]]$data$pval,res[[i]]$data$fdr,res[[i]]$data$overlap,res[[i]]$data$signature)
  return(output)
},res))
hallmarkhyp.cnv <- data.frame(pathway=hallmarkhyp1[,1],pval = as.numeric(hallmarkhyp1[,2]),fdr = as.numeric(hallmarkhyp1[,3]),overlap=as.numeric(hallmarkhyp1[,4]),signature =as.numeric(hallmarkhyp1[,5]))

hyp_obj <- hypeR(MP, list(gene = setdiff(clone.cor$gene,clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05])))
res <- hyp_obj$data
hallmarkhyp1 <- do.call(rbind,lapply(1:length(res), function(i,res){
  pathway <- names(res)[i]
  output <- c(pathway,res[[i]]$data$pval,res[[i]]$data$fdr,res[[i]]$data$overlap,res[[i]]$data$signature)
  return(output)
},res))
hallmarkhyp.peak <- data.frame(pathway=hallmarkhyp1[,1],pval = as.numeric(hallmarkhyp1[,2]),fdr = as.numeric(hallmarkhyp1[,3]),overlap=as.numeric(hallmarkhyp1[,4]),signature =as.numeric(hallmarkhyp1[,5]))

hallmarkhyp$group <- "All"
hallmarkhyp$size <- dim(clone.cor)[1]
hallmarkhyp.cnv$group <- "CNV"
hallmarkhyp.cnv$size <- length(clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05])
hallmarkhyp.peak$group <- "peak"
hallmarkhyp.peak$size <- dim(clone.cor)[1]-length(clone.cor$gene[clone.cor$RNA.r > 0 & clone.cor$RNA.p < 0.05])
MP.gsea <- rbind(hallmarkhyp,hallmarkhyp.cnv,hallmarkhyp.peak)
MP.gsea$ratio <- MP.gsea$overlap/MP.gsea$size
MP.gsea1 <- MP.gsea[MP.gsea$group!="All",]
MP.gsea1 <- MP.gsea1[MP.gsea1$pval < 0.01,]
MP.gsea1 <- MP.gsea1[order(MP.gsea1$group,MP.gsea1$ratio),]
pathway <- unique(MP.gsea1$pathway)
index <- match(pathway,MP.gsea1$pathway)
MP.gsea1 <- MP.gsea1[index,]
MP.gsea1$pathway <- factor(MP.gsea1$pathway,levels = pathway)

p1 <- ggplot(data=hallmark.gsea1, aes(x=pathway, y=ratio,fill=group)) +
  geom_bar(stat="identity")+ coord_flip()+
  scale_fill_manual(values = c("#FDCFA2","#9ECAE1"))

p2 <- ggplot(data=MP.gsea1, aes(x=pathway, y=ratio,fill=group)) +
  geom_bar(stat="identity")+ coord_flip()+
  scale_fill_manual(values = c("#FDCFA2","#9ECAE1"))
pdf(file=paste0(respath,"/variable.gene.pathway.enriched.pdf"),width = 12,height = 6)
p1+p2
dev.off()


