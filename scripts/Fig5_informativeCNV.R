
####identifying clonal specific CNA region
clonalCN.Out <- function(ratiodata,binID){
  CN <- do.call(rbind,lapply(ratiodata, function(ele,binID){
    bindata <- ele$input_BinSegRatio
    bindata <- bindata[match(binID,bindata$binID),]
    segdata <- ele$seg.dat
    binindex <- match(bindata$segName,segdata$segName)
    bindata$CN <- segdata$integerCN[binindex]
    return(bindata$CN)
  },binID))
  colnames(CN) <- binID
  row.names(CN) <- names(ratiodata)
  return(CN)
}

sigCNA.clone <- function(celldata,cellinfo,cnh,breaks){
  celldata <- celldata[,colnames(celldata)%in%colnames(cnh)]
  breakindex <- which(breaks == 1)
  stat.test <- do.call(rbind,lapply(1:length(breakindex), function(i,breakindex,celldata,cellinfo,cnh){
    if (i != length(breakindex)){
      subdata <- celldata[,breakindex[i]:(breakindex[i+1]-1)]
      subCN <- cnh[,breakindex[i]:(breakindex[i+1]-1)]
      start <- colnames(celldata)[breakindex[i]]
      end <- colnames(celldata)[breakindex[i+1]-1]
    }else{
      subdata <- celldata[,breakindex[i]:dim(celldata)[2]]
      subCN <- cnh[,breakindex[i]:dim(celldata)[2]]
      start <- colnames(celldata)[breakindex[i]]
      end <- colnames(celldata)[dim(celldata)[2]]
    }
    if (start==end){
      subdata.mean <- subdata
      subCN.mean <- subCN
    }else{
      subdata.mean <- apply(subdata,1,mean)
      subCN.mean <- apply(subCN,1,mean)
    }
    index <- match(cellinfo$cellname,names(subdata.mean))
    cellinfo$value <- subdata.mean[index]
    stat.test <- cellinfo %>%
      t_test(value ~ clone) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
    stat.mean <- cellinfo %>%
      group_by(clone) %>%
      summarise_at(vars(value), list(name = mean))
    index1 <- match(stat.test$group1,stat.mean$clone)
    stat.test$group1.ratio <- stat.mean$name[index1]
    index2 <- match(stat.test$group2,stat.mean$clone)
    stat.test$group2.ratio <- stat.mean$name[index2]
    index1 <- match(stat.test$group1,names(subCN.mean))
    stat.test$group1.CN <- subCN.mean[index1]
    index2 <- match(stat.test$group2,names(subCN.mean))
    stat.test$group2.CN <- subCN.mean[index2]
    stat.test$FC.ratio <- stat.test$group1.ratio-stat.test$group2.ratio
    stat.test$FC.CN <- stat.test$group1.CN-stat.test$group2.CN
    stat.test$diff <- stat.test$FC.ratio*stat.test$FC.CN
    stat.test$start <- start
    stat.test$end <- end
    return(stat.test)
  },breakindex,celldata,cellinfo,cnh))
  return(stat.test)
}



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

work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
datapath <- paste0(work_dir,'/TeaCNV_HCC_MultiOmic')
setwd(datapath)
samplelist <- list.files()
samplelist <- samplelist[grep("JS",samplelist)]
respath <- paste0(work_dir,"/clonal_complexity_res/informativeCNV")
ifelse(!dir.exists(file.path(respath)), dir.create(file.path(respath)), FALSE)


informativeCNV <- list()

sampleID <- "JS230712J1"
####Extract EiCNA estimation
CNVres <- readRDS(paste0("./",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))

cellinfo <- CNVres$cellinfo
cellinfo <- cellinfo[!is.na(cellinfo$clone),]
celldata <- CNVres$cellbinRatio_deNoise
celldata <- celldata[,cellinfo$cellname]
ratiodata <- CNVres$clonalest
binID <- row.names(celldata)
CNdata <- clonalCN.Out(ratiodata,binID)
cloneID <- as.character(unique(cellinfo$clone))
metacloneCN <- CNdata
row.names(metacloneCN) <- cloneID
colnames(metacloneCN) <- binID
cnh <- do.call(cbind,apply(metacloneCN, 2, function(x){
  if (!(NA %in% x)){
    return(x)
  }
}))
breaks <- c(1)
breaks <- c(breaks,unlist(lapply(1:(dim(cnh)[2]-1), function(j,cnh){
  if (sum(cnh[,j]==cnh[,j+1])== dim(cnh)[1]){
    return(0)
  }else{
    return(1)
  }
},cnh)))
celldata <- t(celldata)
seg.res <- sigCNA.clone(celldata,cellinfo,cnh,breaks)
sig.CNA <- as.data.frame(seg.res[seg.res$p.adj < 0.01&seg.res$diff > 0,])
sig.CNA$startpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$start,split="-"))[,2])
sig.CNA$endpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$end,split="-"))[,3])
sig.CNA$seglength <- (sig.CNA$endpos-sig.CNA$startpos)/1000000
sig.CNA <- sig.CNA[sig.CNA$seglength > 10,]
sig.seg <- unique(sig.CNA[,c("start","end")])
info.cnh <- apply(sig.seg, 1, function(x,cnh){
  index1 <- which(colnames(cnh)==x[1])
  index2 <- which(colnames(cnh)==x[2])
  subdata <- cnh[,index1:index2]
  return(apply(subdata, 1, mean))
},cnh)
informativeCNV[[length(informativeCNV)+1]] <- sig.CNA
names(informativeCNV)[length(informativeCNV)] <- "JS230712J1"
colnames(info.cnh) <- paste(do.call(rbind,strsplit(sig.seg$start,split="-"))[,1],do.call(rbind,strsplit(sig.seg$start,split="-"))[,2],do.call(rbind,strsplit(sig.seg$end,split="-"))[,3],sep="-")

dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
#treeNJ  <- NJ(dm2)
#cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
#PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
#PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
#plot(treeNJ,main="NJ",cex=2)
#plot(PARtree1,main="MP",cex=2)
dev.off()
pdf(file=paste0(respath,"/",sampleID,".cloneTree.heatmap.pdf"),width = 10,height = 5)
pheatmap::pheatmap(info.cnh,cluster_rows = F,cluster_cols = F)
dev.off()


####

sampleID <- "JS230725J1"
####Extract EiCNA estimation
CNVres <- readRDS(paste0(datapath,"/",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
cellinfo <- CNVres$cellinfo
cellinfo <- cellinfo[!is.na(cellinfo$clone),]
celldata <- CNVres$cellbinRatio_deNoise
celldata <- celldata[,cellinfo$cellname]
ratiodata <- CNVres$clonalest
binID <- row.names(celldata)
CNdata <- clonalCN.Out(ratiodata,binID)
cloneID <- as.character(unique(cellinfo$clone))
metacloneCN <- CNdata
row.names(metacloneCN) <- cloneID
colnames(metacloneCN) <- binID
cnh <- do.call(cbind,apply(metacloneCN, 2, function(x){
  if (!(NA %in% x)){
    return(x)
  }
}))
breaks <- c(1)
breaks <- c(breaks,unlist(lapply(1:(dim(cnh)[2]-1), function(j,cnh){
  if (sum(cnh[,j]==cnh[,j+1])== dim(cnh)[1]){
    return(0)
  }else{
    return(1)
  }
},cnh)))
celldata <- t(celldata)
seg.res <- sigCNA.clone(celldata,cellinfo,cnh,breaks)
sig.CNA <- as.data.frame(seg.res[seg.res$p.adj < 0.01&seg.res$diff > 0,])
sig.CNA$startpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$start,split="-"))[,2])
sig.CNA$endpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$end,split="-"))[,3])
sig.CNA$seglength <- (sig.CNA$endpos-sig.CNA$startpos)/1000000
sig.CNA <- sig.CNA[sig.CNA$seglength > 5,]
informativeCNV[[length(informativeCNV)+1]] <- sig.CNA
names(informativeCNV)[length(informativeCNV)] <- "JS230725J1"
sig.seg <- unique(sig.CNA[,c("start","end")])
info.cnh <- apply(sig.seg, 1, function(x,cnh){
  index1 <- which(colnames(cnh)==x[1])
  index2 <- which(colnames(cnh)==x[2])
  subdata <- cnh[,index1:index2]
  return(apply(subdata, 1, mean))
},cnh)

colnames(info.cnh) <- paste(do.call(rbind,strsplit(sig.seg$start,split="-"))[,1],do.call(rbind,strsplit(sig.seg$start,split="-"))[,2],do.call(rbind,strsplit(sig.seg$end,split="-"))[,3],sep="-")

dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
treeNJ  <- NJ(dm2)
cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
plot(treeNJ,main="NJ",cex=2)
plot(PARtree1,main="MP",cex=2)
dev.off()
pdf(file=paste0(respath,"/",sampleID,".cloneTree.heatmap.pdf"),width = 10,height = 5)
pheatmap::pheatmap(info.cnh[c("5","4","2","1","3","6"),],cluster_rows = F,cluster_cols = F)
dev.off()

#select Informative segments for buliding tree
infoSegs <- c("chr1-233613204-248907521","chr3-15606371-87228012","chr4-123842-13485207","chr11-206917-5685346","chr11-69891281-104969829","chr22-31970093-38424595")
info.cnh <- info.cnh[,infoSegs]
dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
treeNJ  <- NJ(dm2)
cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree-InfoSegs.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
plot(treeNJ,main="NJ",cex=2)
plot(PARtree1,main="MP",cex=2)
dev.off()



sampleID <- "JS230830J2"
####Extract EiCNA estimation
CNVres <- readRDS(paste0(datapath,"/",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
cellinfo <- CNVres$cellinfo
cellinfo <- cellinfo[!is.na(cellinfo$clone),]
celldata <- CNVres$cellbinRatio_deNoise
celldata <- celldata[,cellinfo$cellname]
ratiodata <- CNVres$clonalest
binID <- row.names(celldata)
CNdata <- clonalCN.Out(ratiodata,binID)
cloneID <- as.character(unique(cellinfo$clone))
metacloneCN <- CNdata
row.names(metacloneCN) <- cloneID
colnames(metacloneCN) <- binID
cnh <- do.call(cbind,apply(metacloneCN, 2, function(x){
  if (!(NA %in% x)){
    return(x)
  }
}))
breaks <- c(1)
breaks <- c(breaks,unlist(lapply(1:(dim(cnh)[2]-1), function(j,cnh){
  if (sum(cnh[,j]==cnh[,j+1])== dim(cnh)[1]){
    return(0)
  }else{
    return(1)
  }
},cnh)))
celldata <- t(celldata)
seg.res <- sigCNA.clone(celldata,cellinfo,cnh,breaks)
sig.CNA <- as.data.frame(seg.res[seg.res$p.adj < 0.01&seg.res$diff > 0,])
sig.CNA$startpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$start,split="-"))[,2])
sig.CNA$endpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$end,split="-"))[,3])
sig.CNA$seglength <- (sig.CNA$endpos-sig.CNA$startpos)/1000000
sig.CNA <- sig.CNA[sig.CNA$seglength > 10,]
informativeCNV[[length(informativeCNV)+1]] <- sig.CNA
names(informativeCNV)[length(informativeCNV)] <- "JS230830J2"

sig.seg <- unique(sig.CNA[,c("start","end")])
info.cnh <- apply(sig.seg, 1, function(x,cnh){
  index1 <- which(colnames(cnh)==x[1])
  index2 <- which(colnames(cnh)==x[2])
  subdata <- cnh[,index1:index2]
  return(apply(subdata, 1, mean))
},cnh)

colnames(info.cnh) <- paste(do.call(rbind,strsplit(sig.seg$start,split="-"))[,1],do.call(rbind,strsplit(sig.seg$start,split="-"))[,2],do.call(rbind,strsplit(sig.seg$end,split="-"))[,3],sep="-")

dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
treeNJ  <- NJ(dm2)
cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
plot(treeNJ,main="NJ",cex=2)
plot(PARtree1,main="MP",cex=2)
dev.off()
pdf(file=paste0(respath,"/",sampleID,".cloneTree.heatmap.pdf"),width = 10,height = 5)
pheatmap::pheatmap(info.cnh[c("5","3","1","2","4","0"),],cluster_rows = F,cluster_cols = F)
dev.off()




sampleID <- "JS230907J3"
####Extract EiCNA estimation
CNVres <- readRDS(paste0(datapath,"/",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
cellinfo <- CNVres$cellinfo
cellinfo <- cellinfo[!is.na(cellinfo$clone),]
celldata <- CNVres$cellbinRatio_deNoise
celldata <- celldata[,cellinfo$cellname]
ratiodata <- CNVres$clonalest
binID <- row.names(celldata)
CNdata <- clonalCN.Out(ratiodata,binID)
cloneID <- as.character(unique(cellinfo$clone))
metacloneCN <- CNdata
row.names(metacloneCN) <- cloneID
colnames(metacloneCN) <- binID
cnh <- do.call(cbind,apply(metacloneCN, 2, function(x){
  if (!(NA %in% x)){
    return(x)
  }
}))
breaks <- c(1)
breaks <- c(breaks,unlist(lapply(1:(dim(cnh)[2]-1), function(j,cnh){
  if (sum(cnh[,j]==cnh[,j+1])== dim(cnh)[1]){
    return(0)
  }else{
    return(1)
  }
},cnh)))
celldata <- t(celldata)
seg.res <- sigCNA.clone(celldata,cellinfo,cnh,breaks)
sig.CNA <- as.data.frame(seg.res[seg.res$p.adj < 0.01&seg.res$diff > 0,])
sig.CNA$startpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$start,split="-"))[,2])
sig.CNA$endpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$end,split="-"))[,3])
sig.CNA$seglength <- (sig.CNA$endpos-sig.CNA$startpos)/1000000
sig.CNA <- sig.CNA[sig.CNA$seglength > 10,]
informativeCNV[[length(informativeCNV)+1]] <- sig.CNA
names(informativeCNV)[length(informativeCNV)] <- "JS230907J3"

sig.seg <- unique(sig.CNA[,c("start","end")])
info.cnh <- apply(sig.seg, 1, function(x,cnh){
  index1 <- which(colnames(cnh)==x[1])
  index2 <- which(colnames(cnh)==x[2])
  subdata <- cnh[,index1:index2]
  return(apply(subdata, 1, mean))
},cnh)

colnames(info.cnh) <- paste(do.call(rbind,strsplit(sig.seg$start,split="-"))[,1],do.call(rbind,strsplit(sig.seg$start,split="-"))[,2],do.call(rbind,strsplit(sig.seg$end,split="-"))[,3],sep="-")

dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
treeNJ  <- NJ(dm2)
cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
plot(treeNJ,main="NJ",cex=2)
plot(PARtree1,main="MP",cex=2)
dev.off()
pdf(file=paste0(respath,"/",sampleID,".cloneTree.heatmap.pdf"),width = 10,height = 5)
# pheatmap::pheatmap(info.cnh[c("2","1","7","9","5","3","8","6","4"),],cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(info.cnh[c("5","3","8","6","4","2","1","9","7"),],cluster_rows = F,cluster_cols = F)
dev.off()


sampleID <- "JS230915J1"
####Extract EiCNA estimation
CNVres <- readRDS(paste0(datapath,"/",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
cellinfo <- CNVres$cellinfo
cellinfo <- cellinfo[!is.na(cellinfo$clone),]
celldata <- CNVres$cellbinRatio_deNoise
celldata <- celldata[,cellinfo$cellname]
ratiodata <- CNVres$clonalest
binID <- row.names(celldata)
CNdata <- clonalCN.Out(ratiodata,binID)
cloneID <- as.character(unique(cellinfo$clone))
metacloneCN <- CNdata
row.names(metacloneCN) <- cloneID
colnames(metacloneCN) <- binID
cnh <- do.call(cbind,apply(metacloneCN, 2, function(x){
  if (!(NA %in% x)){
    return(x)
  }
}))
breaks <- c(1)
breaks <- c(breaks,unlist(lapply(1:(dim(cnh)[2]-1), function(j,cnh){
  if (sum(cnh[,j]==cnh[,j+1])== dim(cnh)[1]){
    return(0)
  }else{
    return(1)
  }
},cnh)))
celldata <- t(celldata)
seg.res <- sigCNA.clone(celldata,cellinfo,cnh,breaks)
sig.CNA <- as.data.frame(seg.res[seg.res$p.adj < 0.01&seg.res$diff > 0,])
sig.CNA$startpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$start,split="-"))[,2])
sig.CNA$endpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$end,split="-"))[,3])
sig.CNA$seglength <- (sig.CNA$endpos-sig.CNA$startpos)/1000000
sig.CNA <- sig.CNA[sig.CNA$seglength > 5,]
informativeCNV[[length(informativeCNV)+1]] <- sig.CNA
names(informativeCNV)[length(informativeCNV)] <- "JS230915J1"

sig.seg <- unique(sig.CNA[,c("start","end")])
info.cnh <- apply(sig.seg, 1, function(x,cnh){
  index1 <- which(colnames(cnh)==x[1])
  index2 <- which(colnames(cnh)==x[2])
  subdata <- cnh[,index1:index2]
  return(apply(subdata, 1, mean))
},cnh)

colnames(info.cnh) <- paste(do.call(rbind,strsplit(sig.seg$start,split="-"))[,1],do.call(rbind,strsplit(sig.seg$start,split="-"))[,2],do.call(rbind,strsplit(sig.seg$end,split="-"))[,3],sep="-")

dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
treeNJ  <- NJ(dm2)
cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
plot(treeNJ,main="NJ",cex=2)
plot(PARtree1,main="MP",cex=2)
dev.off()
pdf(file=paste0(respath,"/",sampleID,".cloneTree.heatmap.pdf"),width = 10,height = 5)
pheatmap::pheatmap(info.cnh[c("2","1","4","3"),],cluster_rows = F,cluster_cols = F)
dev.off()


sampleID <- "JS230915J2"
####Extract EiCNA estimation
CNVres <- readRDS(paste0(datapath,"/",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))

cellinfo <- CNVres$cellinfo
cellinfo <- cellinfo[!is.na(cellinfo$clone),]
celldata <- CNVres$cellbinRatio_deNoise
celldata <- celldata[,cellinfo$cellname]
ratiodata <- CNVres$clonalest
binID <- row.names(celldata)
CNdata <- clonalCN.Out(ratiodata,binID)
cloneID <- as.character(unique(cellinfo$clone))
metacloneCN <- CNdata
row.names(metacloneCN) <- cloneID
colnames(metacloneCN) <- binID
cnh <- do.call(cbind,apply(metacloneCN, 2, function(x){
  if (!(NA %in% x)){
    return(x)
  }
}))
breaks <- c(1)
breaks <- c(breaks,unlist(lapply(1:(dim(cnh)[2]-1), function(j,cnh){
  if (sum(cnh[,j]==cnh[,j+1])== dim(cnh)[1]){
    return(0)
  }else{
    return(1)
  }
},cnh)))
celldata <- t(celldata)
seg.res <- sigCNA.clone(celldata,cellinfo,cnh,breaks)
sig.CNA <- as.data.frame(seg.res[seg.res$p.adj < 0.01&seg.res$diff > 0,])
sig.CNA$startpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$start,split="-"))[,2])
sig.CNA$endpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$end,split="-"))[,3])
sig.CNA$seglength <- (sig.CNA$endpos-sig.CNA$startpos)/1000000
sig.CNA <- sig.CNA[sig.CNA$seglength > 5,]
informativeCNV[[length(informativeCNV)+1]] <- sig.CNA
names(informativeCNV)[length(informativeCNV)] <- "JS230915J2"

sig.seg <- unique(sig.CNA[,c("start","end")])
info.cnh <- apply(sig.seg, 1, function(x,cnh){
  index1 <- which(colnames(cnh)==x[1])
  index2 <- which(colnames(cnh)==x[2])
  subdata <- cnh[,index1:index2]
  return(apply(subdata, 1, mean))
},cnh)

colnames(info.cnh) <- paste(do.call(rbind,strsplit(sig.seg$start,split="-"))[,1],do.call(rbind,strsplit(sig.seg$start,split="-"))[,2],do.call(rbind,strsplit(sig.seg$end,split="-"))[,3],sep="-")

dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
treeNJ  <- NJ(dm2)
cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
plot(treeNJ,main="NJ",cex=2)
plot(PARtree1,main="MP",cex=2)
dev.off()
pdf(file=paste0(respath,"/",sampleID,".cloneTree.heatmap.pdf"),width = 10,height = 5)
#pheatmap::pheatmap(info.cnh[c("2","1","3","4"),],cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(info.cnh[c("2","3","1","4"),],cluster_rows = F,cluster_cols = F)

dev.off()


sampleID <- "JS230919J1"
####Extract EiCNA estimation
CNVres <- readRDS(paste0(datapath,"/",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
cellinfo <- CNVres$cellinfo
cellinfo <- cellinfo[!is.na(cellinfo$clone),]
celldata <- CNVres$cellbinRatio_deNoise
celldata <- celldata[,cellinfo$cellname]
ratiodata <- CNVres$clonalest
binID <- row.names(celldata)
CNdata <- clonalCN.Out(ratiodata,binID)
cloneID <- as.character(unique(cellinfo$clone))
metacloneCN <- CNdata
row.names(metacloneCN) <- cloneID
colnames(metacloneCN) <- binID
cnh <- do.call(cbind,apply(metacloneCN, 2, function(x){
  if (!(NA %in% x)){
    return(x)
  }
}))
breaks <- c(1)
breaks <- c(breaks,unlist(lapply(1:(dim(cnh)[2]-1), function(j,cnh){
  if (sum(cnh[,j]==cnh[,j+1])== dim(cnh)[1]){
    return(0)
  }else{
    return(1)
  }
},cnh)))
celldata <- t(celldata)
seg.res <- sigCNA.clone(celldata,cellinfo,cnh,breaks)
sig.CNA <- as.data.frame(seg.res[seg.res$p.adj < 0.01&seg.res$diff > 0,])
sig.CNA$startpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$start,split="-"))[,2])
sig.CNA$endpos <- as.numeric(do.call(rbind,strsplit(sig.CNA$end,split="-"))[,3])
sig.CNA$seglength <- (sig.CNA$endpos-sig.CNA$startpos)/1000000
sig.CNA <- sig.CNA[sig.CNA$seglength > 13,]
informativeCNV[[length(informativeCNV)+1]] <- sig.CNA
names(informativeCNV)[length(informativeCNV)] <- "JS230919J1"

sig.seg <- unique(sig.CNA[,c("start","end")])
info.cnh <- apply(sig.seg, 1, function(x,cnh){
  index1 <- which(colnames(cnh)==x[1])
  index2 <- which(colnames(cnh)==x[2])
  subdata <- cnh[,index1:index2]
  return(apply(subdata, 1, mean))
},cnh)

colnames(info.cnh) <- paste(do.call(rbind,strsplit(sig.seg$start,split="-"))[,1],do.call(rbind,strsplit(sig.seg$start,split="-"))[,2],do.call(rbind,strsplit(sig.seg$end,split="-"))[,3],sep="-")

#select Informative segments for buliding tree
infoSegs <-  colnames(info.cnh)[c(6,7,9,10,12,19,21,23)]
info.cnh <- info.cnh[,infoSegs]


dm2  <- dist(info.cnh)
treeUPGMA  <- upgma(dm2)
treeNJ  <- NJ(dm2)
cnh3=as.phyDat(info.cnh,type="USER",levels=c(0:max(info.cnh)))
PARtree=pratchet(cnh3,minit = 10,maxit = 1000,k=5)
PARtree1 <- acctran(PARtree, cnh3)
pdf(file=paste0(respath,"/",sampleID,".cloneTree.pdf"),width = 6,height = 4)
par(mfrow=c(1,3))
plot(treeUPGMA,main="UPGMA",cex=2)
plot(treeNJ,main="NJ",cex=2)
plot(PARtree1,main="MP",cex=2)
dev.off()
pdf(file=paste0(respath,"/",sampleID,".cloneTree.heatmap.pdf"),width = 10,height = 5)
pheatmap::pheatmap(info.cnh[c("4","5","6","2","3","1"),],cluster_rows = F,cluster_cols = F) #MPtree
dev.off()


saveRDS(informativeCNV,paste0(respath,"/informativeCNV.rds"))


###-----------------------------------------###
###plot density for all candidate Infomative CNVs
###-----------------------------------------###
script_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/code/code_CNV/github/TeaCNV/"
setwd(script_path)
library(TeaCNV)
# source("funs_InformationSegs.r")
# source("fun_Grange.R")
# source("figurePlot.R")
library(ComplexHeatmap)
cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")

informativeCNV <- readRDS(paste0(respath,"/informativeCNV.rds"))
outdir <- datapath
sampleIDs <- names(informativeCNV)
sampleIDs
for(i in 1:length(sampleIDs)){
  clt <- as.character(sampleIDs[i])
  print(clt)
  outdir_clt <- paste0(outdir,'/',clt,"/filtered_ana_binSize1");ifelse(!dir.exists(file.path(outdir_clt)), dir.create(file.path(outdir_clt)), FALSE)

  CNVresFile_sc <- paste0(outdir_clt,"/final.CNVres.rds")
  outres <- readRDS(CNVresFile_sc)  

  sig.CNA <- informativeCNV[[clt]]

  sig.seg <- unique(sig.CNA[,c("start","end")])
  chr <- sapply(strsplit(sig.seg$start,"-|_|:"),"[",1)
  start <- as.numeric(sapply(strsplit(sig.seg$start,"-|_|:"),"[",2))
  end <- as.numeric(sapply(strsplit(sig.seg$end,"-|_|:"),"[",3))
  sig_seg <- paste(chr,start,end,sep="-")
  seg_bed <- data.frame(chr,start,end,sig_seg)


  pl_clt <- densityPlot_segs(outres$clonalest,seg_bed)

  ncol <- ifelse(length(pl_clt)>4,4,length(pl_clt))
  width <- ifelse(length(pl_clt)>4,14,3.5*length(pl_clt))
  pcom <- cowplot::plot_grid(plotlist=pl_clt,ncol = ncol,align="v")
  ggsave(paste0(respath,"/",clt,".DensityPlot_InfoSegs.pdf"),pcom, width=width, height=2.5*ceiling(length(pl_clt)/4),device = pdf,bg="white",limitsize = FALSE)  




  if(clt == "JS230725J1"){
    # seg <- data.frame(chr="chr1",strat=143735634,end = 229099475)
    # pl_clt.seg <- densityPlot_segs(outres$clonalest,seg)
    # ggsave(paste0(respath,"/",clt,".DensityPlot_InfoSegs_chr1q.pdf"),pl_clt.seg[[1]], width=3.5, height=2.5,device = pdf,bg="white",limitsize = FALSE)  
    seg <- seg_bed[c(7,12,14,20,21,29),]
    colOdr <- c("5","4","2","1","3","6")
  }
  if(clt == "JS230830J2"){
    seg <- seg_bed[c(6,22,58,59,66,72,85),]
    colOdr <- c("5","3","1","2","4","0")
  }

  if(clt == "JS230907J3"){
    seg <- seg_bed[c(5,24,46,50,51,52),]
    colOdr <- c("5","3","8","6","4","2","1","9","7")
  }

  if(clt == "JS230915J1"){
      ##select Informative segments
    seg <- seg_bed[c(9,10,11,13,14),]
    colOdr <- c("2","1","4","3")
  }


  if(clt == "JS230915J2"){
    seg_bed2 <- data.frame(chr="chr9",start=9954,end=101640330)
    segnm <- paste(seg2,collapse ="-")
    pl_clt.seg <- densityPlot_segs(outres$clonalest,seg_bed2)
    ggsave(paste0(respath,"/",clt,".DensityPlot_InfoSegs_",segnm,".pdf"),pl_clt.seg[[1]], width=3.5, height=2.5,device = pdf,bg="white",limitsize = FALSE)  
    #select Informative segments
    seg <- seg_bed[1:4,]
    colOdr <- c("2","3","1","4")
  }
  if(clt == "JS230919J1"){
    seg <- seg_bed[c(6,7,9,10,12,19,21,23),]
    colOdr <- c("4","5","6","2","3","1")
  }


  #heatmap for Informative segmets with true length
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
    CNmt <- na.omit(CNmt)
       
    row_CNmt <- rownames(CNmt)
    chr_CNmt <- sapply(strsplit(row_CNmt,"_|-|:"),"[",1)
    start_CNmt <- sapply(strsplit(row_CNmt,"_|-|:"),"[",2)
    end_CNmt <- sapply(strsplit(row_CNmt,"_|-|:"),"[",3)
    row_CNmt_bed <- data.frame(chr=chr_CNmt,start=start_CNmt,end=end_CNmt)

    row_align <- align_Grange2bin(row_CNmt_bed,seg)
    row_align <- row_align[!is.na(row_align$sig_seg),,drop=FALSE]
    rownames(row_align) <- gsub("_","-",row_align$binID)

    CNmt_filt <- CNmt[rownames(CNmt)%in% rownames(row_align),,drop=F]

    clone_info <- data.frame(row.names=colnames(CNmt_filt),clone=colnames(CNmt_filt))
    color_r <- cols_Palette[1:length(unique(clone_info$clone))]
    names(color_r) <- sort(unique(clone_info$clone))
    clone_info$clone <- factor(clone_info$clone,levels=colOdr)
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r
    height <- ifelse(ncol(CNmt_filt)>2,0.35*(ncol(CNmt_filt))+0.5,ifelse(ncol(CNmt_filt)==1,1.2,1.5))
    p_cloneCN <- heatmap4peakMt(mat=CNmt_filt[,colOdr],
                      meta_info=clone_info[colOdr,,drop=F],
                      sep_by="-",
                      outdir= respath,
                      value.type="CNV",
                      clust_rows=F,
                      show_legend_row = T,
                      legend_titles="integer CN",
                      fileout_name=paste0(clt,".heatmap_cloneCNA_InfoSegs"),
                      col_list=left_anno_cols,
                      column.title = NULL,
                      width=10,height=height,device="pdf") 

    seg$sig_seg

}



