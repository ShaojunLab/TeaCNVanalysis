###Simulate scATAC with different genomic CNV fraction
suppressMessages({
	library(simATAC)
	library(SingleCellExperiment)
	library(GenomicRanges)
	library(BSgenome.Hsapiens.UCSC.hg38)
	library(SummarizedExperiment)
	library(dplyr)
	library(data.table)
	library(plyranges)
	library(Matrix)
})

options(expressions=10000)
wkdir <- "./TeaCNV"
outdir <- paste0(wkdir,"/simulation");if(!file.exists(outdir)){dir.create(outdir,recursive=T)}

outdir_simData <- paste0(outdir,"/simData_diffCNVpct");if(!file.exists(outdir_simData)){dir.create(outdir_simData,recursive=T)}
setwd(outdir_simData)

###------------------------------------###
### Step 1: Caculate CNV fraction from HCC
###------------------------------------###
datapath <- "./TeaCNV"
setwd(datapath)
samplelist <- list.files()
samplelist <- samplelist[grep("JS",samplelist)]

SampleInfo <- read.csv("./SampleInfo.csv")

CNV.stat <- c()
for (sampleID in samplelist){
  CNVres <- readRDS(paste0("./",sampleID,"/final.CNVres.rds"))
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
write.csv(CNV.stat,paste0("CNV.stat_10sam.csv"),row.names=FALSE)

CNV.stat <- CNV.stat[order(CNV.stat$CNVFrac),]

###------------------------------------###
### Step 2: Seed bin-cell matrix selection 
###------------------------------------###
library(TeaCNV)
rawmtx_c <- readRDS(paste0("ATACrawmtx_corBatch_12sam.rds"))
inputdata <- paste0("HCC_MultiOmic_12sam_clean.rds")
obj <- readRDS(inputdata)
obj$group <- obj$CellType.final
obj$group[obj$group=="Hepatocytes"] <- obj$dataset[obj$group=="Hepatocytes"]
obj$group[obj$tissue=="HCC_adjacent"&grepl("JS2307|Cholangiocytes",obj$group)] <- "Adjacent.normal"

cell_anno_raw<- data.frame(row.names=rownames(obj@meta.data),
                           Cluster=obj$group,
                           sampleID = obj$dataset,
                           tissue = obj$tissue,
                           celltype = obj$CellType.final,
                           log10_nFrags_peaks=obj$log10_nFrags_peaks)
DefaultAssay(obj) <- "peaks"
ref_group_names <- "Adjacent.normal"

obj_sub <- subset(obj, CellType.final %in% c("Cholangiocytes","Hepatocytes"))
filtcell_res <- FiltCell.obj(obj_sub,assay = "peaks",filt_by=c("nCount_peaks","nCount_RNA"),
                             filt_per_group ="sampleID",by_zscore.x = F,log_ = TRUE)
Hepato_cells <- filtcell_res$cell_keep
nonHepato_cells <- rownames(cell_anno_raw)[!cell_anno_raw$celltype %in%c("Cholangiocytes","Hepatocytes")]
cells_select <- c(nonHepato_cells,Hepato_cells)
cell_anno <- cell_anno_raw[cells_select,,drop=F]
write.csv(cell_anno,paste0("cell_anno_10sam_input.csv"),row.names=FALSE)
rm(obj)


cells_ref <-  rownames(cell_anno)[cell_anno$Cluster%in%ref_group_names]
counts_ref <- rawmtx_c[,cells_ref]
saveRDS(counts_ref,"seed_counts_ref.rds")



###------------------------------------###
### Step 3: Save seed bin-cell count matrix
###------------------------------------###
CNV.stat <- read.csv(paste0("CNV.stat_10sam.csv"))
CNV.stat <- CNV.stat[order(CNV.stat$CNVFrac),]

cnvpath <- "./TeaCNV"
cnvProp_max <- c(0.3,0.4,0.6,0.8,0.9,1)

##(1) CNVFrac <0.3
for(p in 1:length(cnvProp_max)){
  cnvFrac=cnvProp_max[p] # 
  if(p==1){
    cnvFrac_min <- 0
  }else{cnvFrac_min <- cnvProp_max[p-1] }

  CNV.stat_sub <- CNV.stat[CNV.stat$CNVFrac<=cnvFrac & CNV.stat$CNVFrac > cnvFrac_min,,drop=FALSE]
  sam_tb <- as.data.frame(table(CNV.stat_sub$sampleID))
  sampleIDs <- sam_tb$Var1[sam_tb$Freq>=4] ##select 4 clones
  if(length(sampleIDs)>0){
    sampleID <- sampleIDs[1]
  }else{
    sampleID <- sam_tb$Var1[cumsum(sam_tb$Freq)<=4]
  }


  CNV.stat_sub <- CNV.stat_sub[CNV.stat_sub$sampleID%in% sampleID,,drop=F]

  if(length(sampleID)==1){
    cnvresFile <- paste0(cnvpath,"/",sampleID,"/final.CNVres.rds")
    CNVres <- readRDS(cnvresFile)
    clone_info <- CNVres$cellinfo

    gtCNV <- c()
    cells_in_clone <- list()
    for (i in CNV.stat_sub$clone[1:4]){
      seg.dat <- CNVres$clonalest[[i]]$seg.dat
      bin.dat <- CNVres$clonalest[[i]]$input_BinSegRatio
      bin.dat <- bin.dat[,!grepl("integerCN",colnames(bin.dat))]
      bin.dat <- left_join(bin.dat,seg.dat[,c("segName","integerCN")],by="segName")
      bin.dat <- bin.dat %>%
        dplyr::select(binID,Chromosome,Start,End,segName,integerCN)%>%
        dplyr::mutate(clone=names(CNVres$clonalest)[i],
                      sampleID=sampleID)%>%
        as.data.frame()
      gtCNV <- rbind(gtCNV,bin.dat)  

      cells_clone1 <- clone_info$cellname[clone_info$clone==as.character(i)]
      cells_clone1 <- intersect(colnames(rawmtx_c),cells_clone1)
      cells_in_clone[[length(cells_in_clone)+1]] <- cells_clone1
    }
   
  }else{
    gtCNV <- c()
    cells_in_clone <- list()
    for(s in sampleID){
      clone_use <- CNV.stat_sub$clone[CNV.stat_sub$sampleID==s]
      cnvresFile <- paste0(cnvpath,"/",s,"/final.CNVres.rds")
      CNVres <- readRDS(cnvresFile)
      clone_info <- CNVres$cellinfo

      for (i in clone_use){
        seg.dat <- CNVres$clonalest[[i]]$seg.dat
        bin.dat <- CNVres$clonalest[[i]]$input_BinSegRatio
        bin.dat <- bin.dat[,!grepl("integerCN",colnames(bin.dat))]
        bin.dat <- left_join(bin.dat,seg.dat[,c("segName","integerCN")],by="segName")
        bin.dat <- bin.dat %>%
          dplyr::select(binID,Chromosome,Start,End,segName,integerCN)%>%
          dplyr::mutate(clone=names(CNVres$clonalest)[i],
                        sampleID=s)%>%
          as.data.frame()
        gtCNV <- rbind(gtCNV,bin.dat)  

        cells_clone1 <- clone_info$cellname[clone_info$clone==as.character(i)]
        cells_clone1 <- intersect(colnames(rawmtx_c),cells_clone1)
        cells_in_clone[[length(cells_in_clone)+1]] <- cells_clone1

      }
    }
    bins_common <- gtCNV %>%
      group_by(sampleID) %>%
      summarise(item_set = list(unique(binID))) %>%
      pull(item_set) %>%
      Reduce(intersect, .)
    gtCNV <- gtCNV[gtCNV$binID %in% bins_common,,drop=F]

  }
  write.csv(gtCNV,paste0(outdir_simData,"/GroundTruth_CNVFrac",cnvFrac,".csv"),row.names=FALSE)



  gtCNV$Chromosome <- factor(gtCNV$Chromosome,levels = unique(gtCNV$Chromosome))
  chrom_sizes <- gtCNV %>%
    dplyr::filter(clone==unique(gtCNV$clone)[1])%>%
    group_by(Chromosome) %>%
    summarize(chr_len = max(End)) %>%
    arrange(Chromosome)
  chrom_sizes$offset <- cumsum(c(0, head(chrom_sizes$chr_len, -1)))
  chrom_sizes <- chrom_sizes %>%mutate(line_pos = offset+chr_len/2)
  gtCNV <- gtCNV[!is.na(gtCNV$integerCN),]

  gtCNV <- gtCNV %>%
    left_join(chrom_sizes[, c("Chromosome", "offset")], by = "Chromosome") %>%
    mutate(global_start = Start  + offset)


  pls <- c()
  for(cl in 1:4 ){
    clonei <-  unique(gtCNV$clone)[cl]
    pldata <- gtCNV %>%dplyr::filter(clone==clonei)

    p1 <- ggplot(pldata,aes(x = global_start, y = integerCN)) +
        geom_vline(xintercept = chrom_sizes$offset[2:nrow(chrom_sizes)], color = 'grey', linetype = 'dashed')+
        geom_step(direction = "hv",color = "darkgreen") +
        labs(title = paste0("Genome-wide CNV (clone ",cl,")"), x = "Genome Position", y = "Copy Number") +
        theme_minimal()+
        scale_x_continuous(
          breaks = chrom_sizes$line_pos,
          labels = gsub("chr|chr19|chr21","",unique(pldata$Chromosome))
        ) 
    pls[[clonei]] <- p1 

  }
  pcom <- cowplot::plot_grid(plotlist=pls,ncol = 1,align="v")
  pcom
  ggsave(paste0(outdir_simData,"/GroundTruth_Frac",cnvFrac,".pdf"),pcom,width = 10,height = 10)


  for(cl in 1:length(cells_in_clone)){
    cells_clonei <- cells_in_clone[[cl]]
    counts_clonei <- rawmtx_c[,cells_clonei]
    saveRDS(counts_clonei,paste0(outdir_simData,"/seed_counts_clone",cl,"_CNVFrac",cnvFrac,".rds"))
  }
   

}





###------------------------------------###
### Step 4: Simulation
###------------------------------------###
wkdir <- "./TeaCNV"
outdir <- paste0(wkdir,"/simulation");if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
outdir_simData <- paste0(outdir,"/simData_diffCNVpct");if(!file.exists(outdir_simData)){dir.create(outdir_simData,recursive=T)}
setwd(outdir_simData)
clonal_groups=c("Monoclonal","Biclonal","Triclonal","Tetraclonal")


##4.1 simulate reference data
counts_ref <- readRDS('seed_counts_ref.rds')
seeds <- sample(1:10000, 10)
outdir_ref <- paste0(outdir_simData,"/refCells");if(!file.exists(outdir_ref)){dir.create(outdir_ref,recursive=T)}
object <- simATACEstimate(counts_ref)
# params <- getParameters(object, c("nBins", "species", "nCells"))
# params
###simulate 10 times
for(nth in 1:10){
  object <- setParameters(object, seed = seeds[nth])
  sim_ref <- simATACSimulate(object,default = FALSE,nCells = 2000)
  count_sim <- counts(sim_ref)
  rownames(count_sim) <- rownames(counts_ref)
  # simATACCompare(assay(sim_ref), counts_ref, outdir_ref, 'RefCells')
  saveRDS(count_sim,paste0(outdir_ref,"/counts_ref",nth,".rds"))
}
rm(object)


##Starting Simulation...
size = 2000
cnvProp_max <- c(0.3,0.4,0.6,0.8,0.9,1)

##(1) CNVFrac <0.3
for(p in 1:length(cnvProp_max)){
  cnvFrac=cnvProp_max[p] 
  #cnvFrac=0.3 # 
  outdir_obs <- paste0(outdir_simData,"/CNVpct",cnvFrac);if(!file.exists(outdir_obs)){dir.create(outdir_obs,recursive=T)}

  counts_clone1 <- readRDS(paste0("seed_counts_clone1_CNVFrac",cnvFrac,".rds"))
  counts_clone2 <- readRDS(paste0("seed_counts_clone2_CNVFrac",cnvFrac,".rds"))
  counts_clone3 <- readRDS(paste0("seed_counts_clone3_CNVFrac",cnvFrac,".rds"))
  counts_clone4 <- readRDS(paste0("seed_counts_clone4_CNVFrac",cnvFrac,".rds"))

  object1 <- simATACEstimate(counts_clone1)
  object2 <- simATACEstimate(counts_clone2)
  object3 <- simATACEstimate(counts_clone3)
  object4 <- simATACEstimate(counts_clone4)
dim(counts_clone1)

  ##4.2 Simulation monoclone 
  group = "Monoclonal"
  outdir_obs2 <- paste0(outdir_obs,"/",group);if(!file.exists(outdir_obs2)){dir.create(outdir_obs2,recursive=T)}
  seeds <- sample(1:10000, 10)
  for(nth in 1:10){
    object1 <- setParameters(object1, seed = seeds[nth])
    sim1 <- simATACSimulate(object1,default = FALSE,nCells = size)
    count_sim1 <- counts(sim1)
    rownames(count_sim1) <- rownames(counts_clone1)
    colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))

    saveRDS(count_sim1,paste0(outdir_obs2,"/counts_",nth,".rds"))
  }

  ##4.3 Simulation Biclonal 
  group = "Biclonal"
  outdir_obs2 <- paste0(outdir_obs,"/",group);if(!file.exists(outdir_obs2)){dir.create(outdir_obs2,recursive=T)}
  seeds <- sample(1:10000, 10)
  NcRatio <- c(0.9,0.1)

  for(nth in 1:10){
    object1 <- setParameters(object1, seed = seeds[nth])
    sim1 <- simATACSimulate(object1,default = FALSE,nCells = size*NcRatio[1])
    count_sim1 <- counts(sim1)
    rownames(count_sim1) <- rownames(counts_clone1)
    colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))

    object2 <- setParameters(object2, seed = seeds[nth])
    sim2 <- simATACSimulate(object2,default = FALSE,nCells = size*NcRatio[2])
    count_sim2 <- counts(sim2)
    rownames(count_sim2) <- rownames(counts_clone2)
    colnames(count_sim2) <- paste0("C2_",colnames(count_sim2))

    counts_total <- cbind(count_sim1,count_sim2)
    saveRDS(counts_total,paste0(outdir_obs2,"/counts_",nth,".rds"))
  }

  ##4.4 Simulation Triclonal 
  group = "Triclonal"
  outdir_obs2 <- paste0(outdir_obs,"/",group);if(!file.exists(outdir_obs2)){dir.create(outdir_obs2,recursive=T)}
  seeds <- sample(1:10000, 10)
  NcRatio <- c(0.6,0.3,0.1)

  for(nth in 1:10){
    object1 <- setParameters(object1, seed = seeds[nth])
    sim1 <- simATACSimulate(object1,default = FALSE,nCells = size*NcRatio[1])
    count_sim1 <- counts(sim1)
    rownames(count_sim1) <- rownames(counts_clone1)
    colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))

    object2 <- setParameters(object2, seed = seeds[nth])
    sim2 <- simATACSimulate(object2,default = FALSE,nCells = size*NcRatio[2])
    count_sim2 <- counts(sim2)
    rownames(count_sim2) <- rownames(counts_clone2)
    colnames(count_sim2) <- paste0("C2_",colnames(count_sim2))

    #sample from clone 3
    object3 <- setParameters(object3, seed = seeds[nth])
    sim3 <- simATACSimulate(object3,default = FALSE,nCells = size*NcRatio[3])
    count_sim3 <- counts(sim3)
    rownames(count_sim3) <- rownames(counts_clone3)
    colnames(count_sim3) <- paste0("C3_",colnames(count_sim3))

    counts_total <- cbind(count_sim1,count_sim2,count_sim3)

    saveRDS(counts_total,paste0(outdir_obs2,"/counts_",nth,".rds"))
  }

  ##4.5 Simulation Tetraclonal 
  group = "Tetraclonal"
  outdir_obs2 <- paste0(outdir_obs,"/",group);if(!file.exists(outdir_obs2)){dir.create(outdir_obs2,recursive=T)}
  seeds <- sample(1:10000, 10)
  NcRatio <- c(0.4,0.3,0.2,0.1)

  for(nth in 1:10){
    object1 <- setParameters(object1, seed = seeds[nth])
    sim1 <- simATACSimulate(object1,default = FALSE,nCells = size*NcRatio[1])
    count_sim1 <- counts(sim1)
    rownames(count_sim1) <- rownames(counts_clone1)
    colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))

    object2 <- setParameters(object2, seed = seeds[nth])
    sim2 <- simATACSimulate(object2,default = FALSE,nCells = size*NcRatio[2])
    count_sim2 <- counts(sim2)
    rownames(count_sim2) <- rownames(counts_clone2)
    colnames(count_sim2) <- paste0("C2_",colnames(count_sim2))

    #sample from clone 3
    object3 <- setParameters(object3, seed = seeds[nth])
    sim3 <- simATACSimulate(object3,default = FALSE,nCells = size*NcRatio[3])
    count_sim3 <- counts(sim3)
    rownames(count_sim3) <- rownames(counts_clone3)
    colnames(count_sim3) <- paste0("C3_",colnames(count_sim3))

    #sample from clone 4
    object4 <- setParameters(object4, seed = seeds[nth])
    sim4 <- simATACSimulate(object4,default = FALSE,nCells = size*NcRatio[4])
    count_sim4 <- counts(sim4)
    rownames(count_sim4) <- rownames(counts_clone4)
    colnames(count_sim4) <- paste0("C4_",colnames(count_sim4))

    counts_total <- cbind(count_sim1,count_sim2,count_sim3,count_sim4)
      
    saveRDS(counts_total,paste0(outdir_obs2,"/counts_",nth,".rds"))
  }

}



