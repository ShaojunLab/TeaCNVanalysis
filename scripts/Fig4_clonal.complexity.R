library(data.table)
library(Seurat)
library(patchwork)
library(Signac)
library(dplyr)
library(msigdbr)
library(ggplot2)
library(fgsea)
library(tidyr)
library("tidyverse")
library("purrr")
library(umap) 
library(HelloRanges)
library(reshape)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggseqlogo)
library(chromVAR)
library(cowplot) 

cols_clone <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'

datapath <- paste0(work_dir,'/TeaCNV_HCC_MultiOmic')
setwd(datapath)
samplelist <- list.files()
samplelist <- samplelist[grep("JS",samplelist)]

respath <- paste0(work_dir,"/clonal_complexity_res")
ifelse(!dir.exists(file.path(respath)), dir.create(file.path(respath)), FALSE)
subclone.summ <- c()
for (sampleID in samplelist){
  CNVres <- readRDS(paste0("./",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
  cellinfo <- CNVres$cellinfo
  cellinfo <- cellinfo[!is.na(cellinfo$clone),]
  subres <- as.data.frame(table(cellinfo$clone))
  if(min(as.numeric(as.character(subres$Var1)))==0){
    subres$Var1 <- as.numeric(as.character(subres$Var1))+1
  }
  subres$sampleID <- sampleID
  ploidy <- c()
  for (ele in CNVres$clonalest){
    ploidy <- c(ploidy,ele$ploidy)
  }
  subres$ploidy <- ploidy
  subclone.summ <- rbind(subclone.summ,subres)

}
n.subclone <- subclone.summ %>%
  group_by(sampleID) %>%
  summarise(count=n())
n.subclone <- as.data.frame(n.subclone)
n.subclone <- n.subclone[order(n.subclone$count),]
n.subclone$sampleID <- factor(n.subclone$sampleID,levels = n.subclone$sampleID)
n.subclone$clonality <- n.subclone$count

#add sample label
SampleInfo <- read.csv(paste0(work_dir,"/AnalysisData/SampleInfo.csv"))
n.subclone <- left_join(n.subclone,SampleInfo[,1:2],by="sampleID")
n.subclone$sampleID_new <- factor(n.subclone$sampleID_new,levels = c(paste0("HCC",seq(1:10))))
subclone.summ <- left_join(subclone.summ,SampleInfo[,1:2],by="sampleID")

Nampli <- nrow(subclone.summ[subclone.summ$ploidy>2.5,,drop=F])
print(paste0("Total number of clones: ",length(unique(subclone.summ$Var1))))
#40
print(paste0(Nampli," clones had ploidy greater than 2.5"))
#28


plot_theme <- theme(plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 7),
    axis.ticks=element_line(color="black",linewidth=0.5),
    axis.text=element_text(color="black",size=7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.background=element_blank(),
    legend.key=element_blank(),
    legend.text=element_text(color="black",size=7),
    legend.title=element_blank())


p1 <- ggplot(data=n.subclone, aes(x=sampleID_new, y=clonality)) +
  geom_bar(stat="identity",fill = "grey",  width = 0.5)+
  theme_classic() +
  plot_theme+
  labs(y = "Clonality")

subclone.summ$sampleID <- factor(subclone.summ$sampleID,levels = n.subclone$sampleID)
subclone.summ$sampleID_new <- factor(subclone.summ$sampleID_new,levels = n.subclone$sampleID_new)

size.subclone <- subclone.summ %>%
  group_by(sampleID_new) %>%
  summarise(size=sum(Freq))
size.subclone <- as.data.frame(size.subclone)
subclone.summ$size <- size.subclone$size[match(subclone.summ$sampleID_new,size.subclone$sampleID_new)]
subclone.summ$pro <- subclone.summ$Freq/subclone.summ$size
max.subclone <- subclone.summ %>%
  group_by(sampleID_new) %>%
  # filter(Var1==1)%>%
  summarise(size=max(pro))

#order clone group_by sample
subclone.summ$sampleID_new <- factor(subclone.summ$sampleID_new,levels = as.character(max.subclone$sampleID_new[order(max.subclone$size,decreasing = T)]))
# subclone.summ$Var1 <- as.numeric(as.character(subclone.summ$Var1))
# subclone.summ$Var1 <- factor(subclone.summ$Var1)
subclone.summ <- subclone.summ %>%
  group_by(sampleID_new) %>%
  arrange(desc(pro), .by_group = TRUE)%>%
  mutate(Var2=paste0(sampleID_new,"_",Var1))%>%
  mutate(fill = factor(Var2, levels = Var2))%>%
  as.data.frame()
color_fill <- c("#D1D2D3","#939597","#D5E8F7","#8BC4EB","#9B98C6","#FEDBAE","#FBAF3F","#CEE5AD","#A1CE61","#C85248")
names(color_fill) <- unique(subclone.summ$Var1)
subclone.summ$color <- color_fill[match(subclone.summ$Var1,names(color_fill))]

diversity <- unlist(lapply(samplelist, function(ID){
  -sum(subclone.summ$pro[subclone.summ$sampleID==ID]*log(subclone.summ$pro[subclone.summ$sampleID==ID],base = 2))
}))
names(diversity) <- samplelist
n.subclone$diversity <- diversity[as.character(n.subclone$sampleID)]
n.subclone$sampleID_new <- factor(n.subclone$sampleID_new,levels = as.character(n.subclone$sampleID_new[order(n.subclone$diversity)]))
ploidyCV <- unlist(lapply(samplelist, function(ID){
  sd(subclone.summ$ploidy[subclone.summ$sampleID==ID])/mean(subclone.summ$ploidy[subclone.summ$sampleID==ID])
}))
names(ploidyCV) <- samplelist
n.subclone$ploidyCV <- ploidyCV[as.character(n.subclone$sampleID)]
n.subclone$ploidyCV[is.na(n.subclone$ploidyCV)] <- 0

write.csv(subclone.summ,paste0(respath,"/subclone.summ.res.csv"))
# subclone.summ <- read.csv(paste0(respath,"/subclone.summ.res.csv"),row.names=1)


p2 <- ggplot(data=n.subclone, aes(x=sampleID_new, y=diversity)) +
  geom_bar(stat="identity",  width = 0.5)+
  theme_classic() +
  plot_theme+
  labs(y = "Diversity")

p3 <- ggplot(data=subclone.summ, aes(x=sampleID_new, y=pro, fill=fill)) +
  geom_bar(stat="identity", position = "stack",width = 0.5)+
  theme_classic() +
  theme(legend.position = "none")+
  plot_theme+
  labs(y = "Proportion")+
  scale_fill_manual(values = subclone.summ$color)

p4 <- ggplot(data=subclone.summ, aes(x=sampleID_new, y=ploidy,color=Var1)) +
  geom_point(size = 2)+theme_classic() +
  plot_theme+
  scale_color_manual(values = color_fill)
legend <- get_legend(p4)
p3_com <- (p3+legend)+plot_layout(ncol = 2, widths = c(3,1))


color_sample <- c("#A9011B","#E4A826","#D2DCAD","#DCD6B2","#94BEA7","#4E7989","#75A0AE","#80944E","#547DB1","#9055A2","#D43B36")

names(color_sample) <- level(n.subclone$sampleID_new)
p5 <- ggplot(data=n.subclone, aes(x=clonality, y=ploidyCV,color=sampleID_new)) +
  geom_point(size=2)+theme_classic() +
  theme(plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid = element_blank(), 
    # axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 7),
    axis.ticks=element_line(color="black",linewidth=0.5),
    axis.text=element_text(color="black",size=7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.background=element_blank(),
    legend.key=element_blank(),
    legend.text=element_text(color="black",size=7),
    legend.title=element_blank())+
  scale_color_manual(values = color_sample)

pdf(file=paste0(respath,"/subclonecomplexity.pdf"),width = 16,height = 4)
(p1+p2+p3_com+p4+p5)+plot_layout(ncol = 5, widths = c(1,1,1.5,1,1))
dev.off()


#####
CNV.stat <- c()
for (sampleID in samplelist){
  CNVres <- readRDS(paste0("./",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
  for (i in 1:length(CNVres$clonalest)){
    seg.dat <- CNVres$clonalest[[i]]$seg.dat
    seg.dat$start <- as.numeric(seg.dat$start)
    seg.dat$end <- as.numeric(seg.dat$end)
    seg.dat$integerCN <- as.numeric(seg.dat$integerCN)
    CNVseg <- seg.dat[seg.dat$integerCN !=2,]
    genoFrac <- sum(CNVseg$end-CNVseg$start)/sum(seg.dat$end-seg.dat$start)
    gainFrac <- sum(CNVseg$end[CNVseg$integerCN>2]-CNVseg$start[CNVseg$integerCN>2])/sum(seg.dat$end-seg.dat$start)
    lossFrac <- sum(CNVseg$end[CNVseg$integerCN<2]-CNVseg$start[CNVseg$integerCN<2])/sum(seg.dat$end-seg.dat$start)
    a <- c(i,sampleID,genoFrac,gainFrac,lossFrac)
    CNV.stat <- rbind(CNV.stat,a)
  }
}

CNV.stat <- data.frame(clone=CNV.stat[,1],sampleID=CNV.stat[,2],CNVFrac = as.numeric(CNV.stat[,3]),gainFrac=as.numeric(CNV.stat[,4]),lossFrac=as.numeric(CNV.stat[,5]))
CNV.stat <- left_join(CNV.stat,SampleInfo[,1:2],by="sampleID")

CNV.stat$sampleID_new <- factor(CNV.stat$sampleID_new,levels=unique(CNV.stat$sampleID_new[order(CNV.stat$CNVFrac)]))

p <- ggplot(data=CNV.stat, aes(x=sampleID_new, y=CNVFrac, fill=clone)) +
  geom_bar(stat="identity", position=position_dodge())+
  # scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=cols_clone)+
  theme_classic() +
  plot_theme

pdf(file=paste0(respath,"/CNV.farction.pdf"),width = 5,height = 5)
p
dev.off()

CNV.stat1 <- CNV.stat[,c(1,4,5,6)]
CNV.stat1$lossFrac <- -CNV.stat1$lossFrac
CNV.statlong <- melt(setDT(CNV.stat1), id.vars = c("clone","sampleID_new"), variable.name = "Frac")
CNV.statlong$ID <- paste0(CNV.statlong$sampleID_new,":",CNV.statlong$clone)
p <- ggplot(CNV.statlong) + 
  aes(x = sampleID_new, y = value, fill = clone) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed')+
  # scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=cols_clone)+
  theme_classic() +
  plot_theme

pdf(file=paste0(respath,"/CNV.gain.loss.farction.pdf"),width = 5,height = 5)
p
dev.off()


CNV.length <- c()
allCNV <- c()
for (sampleID in samplelist){
  CNVres <- readRDS(paste0("./",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
  for (i in 1:length(CNVres$clonalest)){
    seg.dat <- CNVres$clonalest[[i]]$seg.dat
    seg.dat$start <- as.numeric(seg.dat$start)
    seg.dat$end <- as.numeric(seg.dat$end)
    seg.dat$integerCN <- as.numeric(seg.dat$integerCN)
    CNVseg <- seg.dat[seg.dat$integerCN !=2,]
    allCNV <- rbind(allCNV,CNVseg[,c("chr","start","end","abspos","absend","integerCN")])
    CNV.length <- rbind(CNV.length,CNVseg[,c(1:3,ncol(CNVseg))])
  }
}
CNV.length$length <- (CNV.length$end-CNV.length$start)/1000000
write.csv(CNV.length,paste0(respath,"/CNV_seg.length.csv"))

### 
CNV.length_stat <- data.frame(
  min=min(CNV.length$length),
  median=median(CNV.length$length),
  max=max(CNV.length$length))
CNV.length_stat
      # min   median      max
# 2.14383 20.51314 146.9686


p <- ggplot(CNV.length, aes(x=length)) + 
  geom_histogram(aes(y=..density..), binwidth=0.2, fill="#56B4E9")+
  theme_classic()

pdf(file=paste0(respath,"/CNV.length.pdf"),width = 5,height = 5)
p
dev.off()



chrom <- unique(allCNV$chr)
gain <- c()
loss <- c()
for (chr in chrom){
  subseg <-allCNV[allCNV$chr==chr,]
  subseg <- subseg[order(subseg$start),]
  start <- unique(subseg$start)
  end <- unique(subseg$end)
  breakpoint <- unique(c(start,end))
  breakpoint <- breakpoint[order(breakpoint)]
  for (i in 2:length(breakpoint)){
    subdata <-  subseg[subseg$start<=breakpoint[i-1]&subseg$end>=breakpoint[i],]
    gaindata <- subdata[subdata$integerCN > 2,]
    lossdata <- subdata[subdata$integerCN < 2,]
    if (dim(gaindata)[1]>0){
      gvalue <- sum(gaindata$integerCN-2) 
    }else{
      gvalue <- 0
    }
    if (dim(lossdata)[1] > 0){
      lvalue <- sum(lossdata$integerCN-2) 
    }else{
      lvalue <- 0
    }
    gain <- rbind(gain,c(chr,breakpoint[i-1],breakpoint[i],gvalue))
    loss <- rbind(loss,c(chr,breakpoint[i-1],breakpoint[i],lvalue))
  }
}
gain <- data.frame(chr=gain[,1],start=as.numeric(gain[,2]),end=as.numeric(gain[,3]),value=as.numeric(gain[,4]))
loss <- data.frame(chr=loss[,1],start=as.numeric(loss[,2]),end=as.numeric(loss[,3]),value=as.numeric(loss[,4]))

allseg <- c()
for (sampleID in samplelist){
  CNVres <- readRDS(paste0("./",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
  for (i in 1:length(CNVres$clonalest)){
    seg.dat <- CNVres$clonalest[[i]]$seg.dat
    seg.dat$start <- as.numeric(seg.dat$start)
    seg.dat$end <- as.numeric(seg.dat$end)
    allseg <- rbind(allseg,seg.dat[,1:3])
  }
}

chrom <- unique(allseg$chr)
chr.length <- unlist(lapply(chrom, function(chr,allseg){
  max(allseg$end[allseg$chr==chr])
},allseg))

chr.length <- c(0,chr.length)

gainabs <- do.call(rbind,lapply(unique(gain$chr), function(chr,chrom,chr.length,gain){
  subdata <- gain[gain$chr==chr,]
  index <- which(chrom==chr)
  subdata$absstart <- subdata$start+sum(chr.length[1:index])
  subdata$absend <- subdata$end+sum(chr.length[1:index])
  return(subdata)
},chrom,chr.length,gain))


lossabs <- do.call(rbind,lapply(unique(gain$chr), function(chr,chrom,chr.length,loss){
  subdata <- loss[loss$chr==chr,]
  index <- which(chrom==chr)
  subdata$absstart <- subdata$start+sum(chr.length[1:index])
  subdata$absend <- subdata$end+sum(chr.length[1:index])
  return(subdata)
},chrom,chr.length,loss))

gainabs$CNV <- "gain"
lossabs$CNV <- "loss"
CNVvalue <- rbind(gainabs,lossabs)

CNVvalue$y <- 0
chr.length <- cumsum(chr.length[2:length(chr.length)])
p <- ggplot() + 
  scale_x_continuous(name="Genome") + 
  scale_y_continuous(name="Accumulate CNV") +
  geom_rect(data=CNVvalue, mapping=aes(xmin=absstart, xmax=absend, ymin=y, ymax=value, fill=CNV), alpha=0.5)+
  geom_vline(xintercept = chr.length, color = 'grey', linetype = 'dashed')+
  theme_classic()

pdf(file=paste0(respath,"/CNV.distribution.pdf"),width = 12,height = 3.5)
p
dev.off()

Genes = c("VEGFA","CCND1")
gtf_file <- "~/Library/Mobile Documents/com~apple~CloudDocs/reference/genes.gtf.gz"
annotations_all <- rtracklayer::import(gtf_file)
annotations <- annotations_all[annotations_all$type=="gene",]
Genes_gr <- annotations[annotations$gene_name %in% Genes,]
Genes_co <- data.frame(chr=seqnames(Genes_gr),start=start(Genes_gr),end=end(Genes_gr),gene_name=Genes)
Genes_ano <- align_Grange2bin(CNVvalue,Genes_co)
Genes_ano2 <- na.omit(Genes_ano)
Genes_ano2$seg.length <- Genes_ano2$end-Genes_ano2$start



allCNV.clone <- c()
for (sampleID in samplelist){
  CNVres <- readRDS(paste0("./",sampleID,"/filtered_ana_binSize1/final.CNVres.rds"))
  for (i in 1:length(CNVres$clonalest)){
    seg.dat <- CNVres$clonalest[[i]]$seg.dat
    seg.dat$start <- as.numeric(seg.dat$start)
    seg.dat$end <- as.numeric(seg.dat$end)
    seg.dat$integerCN <- as.numeric(seg.dat$integerCN)
    CNVseg <- seg.dat[seg.dat$integerCN !=2,]
    CNVseg$count <- 1
    allCNV.clone <- rbind(allCNV.clone,CNVseg[,c("chr","start","end","count","integerCN")])
  }
}

gain.count <- c()
loss.count <- c()
for (chr in chrom){
  subseg <-allCNV.clone[allCNV.clone$chr==chr,]
  subseg <- subseg[order(subseg$start),]
  start <- unique(subseg$start)
  end <- unique(subseg$end)
  breakpoint <- unique(c(start,end))
  breakpoint <- breakpoint[order(breakpoint)]
  for (i in 2:length(breakpoint)){
    subdata <-  subseg[subseg$start<=breakpoint[i-1]&subseg$end>=breakpoint[i],]
    gaindata <- subdata[subdata$integerCN > 2,]
    lossdata <- subdata[subdata$integerCN < 2,]
    gvalue <- dim(gaindata)[1]
    lvalue <- dim(lossdata)[1]
    gain.count <- rbind(gain.count,c(chr,breakpoint[i-1],breakpoint[i],gvalue))
    loss.count <- rbind(loss.count,c(chr,breakpoint[i-1],breakpoint[i],lvalue))
  }
}
gain.count <- data.frame(chr=gain.count[,1],start=as.numeric(gain.count[,2]),end=as.numeric(gain.count[,3]),value=as.numeric(gain.count[,4]))
loss.count <- data.frame(chr=loss.count[,1],start=as.numeric(loss.count[,2]),end=as.numeric(loss.count[,3]),value=as.numeric(loss.count[,4]))
loss.count$value <- -loss.count$value
chr.length <- unlist(lapply(chrom, function(chr,allseg){
  max(allseg$end[allseg$chr==chr])
},allseg))

chr.length <- c(0,chr.length)

gainabs.count <- do.call(rbind,lapply(unique(gain.count$chr), function(chr,chrom,chr.length,gain.count){
  subdata <- gain.count[gain.count$chr==chr,]
  index <- which(chrom==chr)
  subdata$absstart <- subdata$start+sum(chr.length[1:index])
  subdata$absend <- subdata$end+sum(chr.length[1:index])
  return(subdata)
},chrom,chr.length,gain.count))


lossabs.count <- do.call(rbind,lapply(unique(loss.count$chr), function(chr,chrom,chr.length,loss.count){
  subdata <- loss.count[loss.count$chr==chr,]
  index <- which(chrom==chr)
  subdata$absstart <- subdata$start+sum(chr.length[1:index])
  subdata$absend <- subdata$end+sum(chr.length[1:index])
  return(subdata)
},chrom,chr.length,loss.count))

gainabs.count$CNV <- "gain"
lossabs.count$CNV <- "loss"
CNV.count <- rbind(gainabs.count,lossabs.count)

CNV.count$y <- 0
chr.length <- cumsum(chr.length[2:length(chr.length)])
CNV.count$Frac <- CNV.count$value/dim(subclone.summ)[1]
mean(gain.count$value)/dim(subclone.summ)[1]
# [1] 0.6095825
mean(loss.count$value)/dim(subclone.summ)[1]
# [1] -0.08494271


p <- ggplot() + 
  scale_x_continuous(name="Genome") + 
  scale_y_continuous(name="Clone frequency") +
  geom_rect(data=CNV.count, mapping=aes(xmin=absstart, xmax=absend, ymin=y, ymax=Frac, fill=CNV), alpha=0.5)+
  geom_vline(xintercept = chr.length, color = 'grey', linetype = 'dashed')+
  theme_classic()

pdf(file=paste0(respath,"/CNV.clone.frequency.pdf"),width = 14,height = 3.5)
p
dev.off()






