# This program will access the teaCNV's performance by compairing the output of different methods.
##Figure 3,4
options(expressions=10000)
suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(Matrix)
    library(data.table)
    library(magrittr)
    #library(IRdisplay)
    library(parallel)
    library(stringr)
    library(logger)
    library(purrr)
    library(patchwork)
    library(glue)
    library(ggridges)
    library(ggtree)
    #library(extraDistr)
    library(EnsDb.Hsapiens.v86)
})
script_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/code/code_CNV/github/TeaCNV/"
work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
setwd(work_dir)
library(TeaCNV)
source("./Code/Function/funs_evaluation.r")
#work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/code/code_CNV/TeaCNVanalysis/'
#setwd(work_dir)
#load supporing files
refbedfile <- file.path(script_path, "data", "hg38.100Kb.windows.sorted.bed")
bed.ref <- read.table(refbedfile,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

cytoBandFile <- file.path(script_path, "data", "cytoBand_hg38.tsv")
hg38pos <- read.table(cytoBandFile, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
colnames(hg38pos) <- c("chr", "start", "end","arm","band")
hg38pos$arm2 <- substr(hg38pos$arm, 1, 1)
hgpos <- hg38pos %>%
    mutate(arm2 = substr(arm, 1, 1)) %>%
    dplyr::filter(band != "acen")%>%
    dplyr::filter(chr %in% c(paste0("chr",seq(1:22))))%>%
    group_by(chr,arm2) %>%
    mutate(seg_start = min(start), seg_end = max(end))%>%as.data.frame()
hgpos <- unique(hgpos[,c("chr","seg_start","seg_end","arm2")])
sizes_GRCh38 <- file.path(script_path, "data", "sizes.cellranger-GRCh38-1.0.0.txt")
chrom_sizes_hg38 <- read.table(sizes_GRCh38, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
colnames(chrom_sizes_hg38) <- c("chr", "size")

gapFile <- file.path(script_path, "data", "gaps_hg38.rda")
load(gapFile)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
gr_gene <- CollapseToLongestTranscript(ranges = annotations)

# set up paths
outdir = paste0(work_dir,'/EvaluationData')
if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)
data_path = paste0(work_dir, '/AnalysisData/')
sampleID = c("ccRCC1","ccRCC3","ccRCC4")
scCnv_all <- list()

# dna files
dnadir <- data_path
sample_ls <- list.files(dnadir)
sample_ls <- sample_ls[grepl("_SegCNV_res.csv",sample_ls)]


####----------------------####
#(I) Analysis for each sample
####----------------------####
for (i in 1:length(sampleID)){    
    sample = sampleID[i]
    print(paste0("Processing sample: ",sample))
    outdir_sample = file.path(outdir,sample)
    if(!file.exists(outdir_sample)){dir.create(outdir_sample,recursive=T)}
    cnapath = data_path
    ##(1)load dna data
    file_dna <- sample_ls[grepl(sample,sample_ls)]
    dnaFilePath <- file.path(dnadir,file_dna)
    if(sample == "ccRCC1"){
        cnvs_dna <- get_dna_bulk(dnaFilePath,bed.ref,raio.hi=1.1,raio.lo=0.9)
    }else{
       cnvs_dna <- get_dna_bulk(dnaFilePath,bed.ref) 
    }
    
    colnames(cnvs_dna)[2:3] <- c("start","end")
    CNdf <- cnvs_dna
    write.csv(CNdf,paste0(outdir_sample,"/",sample,"_DNA_100kb.csv"))


    #(2)load TeaCNV (scATAC)
    cellinfo <- read.csv(paste0(cnapath, '/cell_info_subClust_',sample,'.csv'),row.names=1)
    normalcell <- rownames(cellinfo)[!cellinfo$group %in% "observed"]
    normalcell <- sapply(strsplit(normalcell,"_"),"[",1)
    # ##load clonal CNAs
    clonalCNApath <- paste0(cnapath,"/final.CNVres_",sample,".rds")
    outres <- readRDS(clonalCNApath)
    clonalSize <- data.frame(table(outres$cellinfo$clone))
    colnames(clonalSize) <- c("clone","count")
    clonalSize$CellProportion <- clonalSize$count/sum(clonalSize$count)

    clonalCN_res <- outres$clonalest
    clonalCNA_bin <- c()
    clonalRatio_bin <- c()
    for(clonei in names(clonalCN_res)){
        seg.dat <-  clonalCN_res[[clonei]]$seg.dat
        bin_dat <- clonalCN_res[[clonei]]$input_BinSegRatio
        bin_dat <- unique(bin_dat[,c("segName","SegMean")])
        seg.dat <- merge(seg.dat,bin_dat,by="segName",all.x=TRUE)
        seg.dat <- seg.dat[,c("chr","start","end","integerCN","SegMean")]
        colnames(seg.dat)[4:5] = c(clonei,paste0(clonei,"_ratio"))
        seg.dat$start <- as.numeric(seg.dat$start)
        seg.dat$end <- as.numeric(seg.dat$end)
        #align TeaCNV to the same bin size as DNA
        TeaCNV <- align_Grange2bin(bed.query=bed.ref,bed.subject=seg.dat)
        clonalCNA <- TeaCNV[,clonei]
        clonalCNA_bin <- cbind(clonalCNA_bin,clonalCNA)
        clonalRatio_bin <- cbind(clonalRatio_bin,TeaCNV[,paste0(clonei,"_ratio")])

    }
    colnames(clonalCNA_bin) <- names(clonalCN_res)
    rownames(clonalCNA_bin) =rownames(clonalRatio_bin) <- paste(bed.ref[,1],bed.ref[,2],bed.ref[,3],sep="_")
    colnames(clonalRatio_bin) <- names(clonalCN_res)
    write.csv(clonalCNA_bin,paste0(cnapath,"/teaCNV_clonalCNA_bin.",sample,".csv"))
    write.csv(clonalRatio_bin,paste0(cnapath,"/teaCNV_clonalRatio_bin.",sample,".csv"))

    # clonalCNA_bin <- read.csv(paste0(cnapath,"/teaCNV_clonalCNA_bin.csv"),row.name=1,check.name=FALSE)
    # clonalRatio_bin <- read.csv(paste0(cnapath,"/teaCNV_clonalRatio_bin.csv"),row.name=1,check.name=FALSE)


    # #save clonal CNA for each bin from TeaCNV
    CNdf <- data.frame(chr=bed.ref[,1],start=bed.ref[,2],end=bed.ref[,3])
    for(cloneNm in names(clonalCN_res)){
        # cloneNm <- "2"
        clonalCNA_bin_mean <- clonalCNA_bin[,cloneNm]
        clonalRatio_bin_mean <- clonalRatio_bin[,cloneNm]

        TeaCNV <- data.frame(chr=bed.ref[,1],start=bed.ref[,2],end=bed.ref[,3],TeaCNV_mean=clonalCNA_bin_mean,TeaCNV_RatioMean=clonalRatio_bin_mean,clone=cloneNm)
        TeaCNV$TeaCNV[TeaCNV$TeaCNV_mean<2]<- 1
        TeaCNV$TeaCNV[TeaCNV$TeaCNV_mean>2 &TeaCNV$TeaCNV_mean<3]<- 3
        TeaCNV$TeaCNV[TeaCNV$TeaCNV_mean>3&!is.na(TeaCNV$TeaCNV_mean)]<- round(TeaCNV$TeaCNV_mean[TeaCNV$TeaCNV_mean>3&!is.na(TeaCNV$TeaCNV_mean)])

        CNdf[,paste0("TeaCNV.",cloneNm)] <- TeaCNV$TeaCNV
        CNdf[,paste0("TeaCNV_mean.",cloneNm)] <-TeaCNV$TeaCNV_mean
        CNdf[,paste0("TeaCNV_RatioMean.",cloneNm)] <-TeaCNV$TeaCNV_RatioMean  
    }
    write.csv(CNdf,paste0(outdir_sample,"/",sample,"_TeaCNV_100kb.csv"))



    ##(3)load epiAneuFinder results
    epiAneu_path <- data_path 
    cna_epiAneu <- read.table(paste0(epiAneu_path, '/epiAneufinder.results_table.',sample,'.tsv'),header=T,row.names=1)
    rownames(cna_epiAneu) <- paste(cna_epiAneu[,1],cna_epiAneu[,2],cna_epiAneu[,3],sep="_")
    cells_epiAneu <- colnames(cna_epiAneu)[4:ncol(cna_epiAneu)]
    cells_epiAneu <- gsub("cell.","",cells_epiAneu)
    cells_epiAneu <- gsub("\\.","-",cells_epiAneu)
    colnames(cna_epiAneu)[4:ncol(cna_epiAneu)] <- cells_epiAneu
    keepIndx <- which(!colnames(cna_epiAneu) %in% normalcell)
    cna_epiAneu <- cna_epiAneu[,c(keepIndx)]

    # generate bulk epiAneuFinder results, whichCNV is determined by the proportion of cells in the gradient
    Psets  <- seq(0.1,1,0.1)
    for(prop in Psets){
        cna_epiAneu_bulk <- get_epiAneufinder_bulk(cna_epiAneu,bed.ref,cnvCellProp.min=prop)
        colnames(cna_epiAneu_bulk)[4]<- "epiAneuRatio"
        write.csv(cna_epiAneu_bulk,paste0(outdir_sample,"/",sample,"_epiAneufinder_100kb_",prop,".csv"))
    }

    #hist(cna_epiAneu_bulk$epiAneuCN,breaks=100,main="epiAneufinder",xlab="bulk epiAneufinder",col="blue",border="black",density=50)

    #(4) load CopyscAT
    copy_path <- data_path
    cnaPath <- paste0(copy_path, '/CopyscAT.cnv_scores.',sample,'.csv')
      
    cna_copyscAT <- get_copyscAT_sc(cnaPath,cytoBandFile)
    cells_copyscAT <- colnames(cna_copyscAT)[4:ncol(cna_copyscAT)]
    cells_copyscAT <- gsub("cell.","",cells_copyscAT)
    cells_copyscAT <- gsub("\\.","-",cells_copyscAT)
    colnames(cna_copyscAT)[4:ncol(cna_copyscAT)] <- cells_copyscAT

    Psets  <- seq(0.1,1,0.1)
    for(prop in Psets){
        cna_copyscAT_bulk_align <- get_copyscAT_bulk(cna_copyscAT,bed.ref,cnvCellProp.min=prop)
        write.csv(cna_copyscAT_bulk_align,paste0(outdir_sample,"/",sample,"_copyscAT_100kb_",prop,".csv"))
    }

    if(sample != "ccRCC1"){
        #(5) load inferCNV
        suppressPackageStartupMessages({
          library(IRanges)
          library(GenomicRanges)
          library(S4Vectors)
          library(EnsDb.Hsapiens.v86)
          library(tidyr)
        })
        obscell <- rownames(cellinfo)[cellinfo$group %in% "observed"]
        inferCNV_path = data_path
        cnvFile1 <- "run.final.infercnv_obj" 
        infercnv_obj = readRDS(paste(inferCNV_path,cnvFile1,sep="/"))
        mtx_cnv <- infercnv_obj@expr.data
        cnaRatio_infercnv_raw <- as.matrix(mtx_cnv[,unlist(infercnv_obj@observation_grouped_cell_indices)])
        cnaRatio_infercnv <- cnaRatio_infercnv_raw[,obscell]
        cnaRatio_infercnv <- as.data.frame(rowMeans(cnaRatio_infercnv,na.rm=TRUE))
        colnames(cnaRatio_infercnv) <- "cnaRatio_infercnv"

        # gene symbol conversion
        genes <- rownames(cnaRatio_infercnv)
        gene_ranges <- gr_gene[gr_gene$gene_name %in% genes,]
        genes_bed <- data.frame(seqnames=seqnames(gene_ranges),start=start(gene_ranges),end=end(gene_ranges),gene_name=gene_ranges$gene_name) %>%
            arrange(as.numeric(seqnames),as.numeric(start)) %>%
            as.data.frame()
        genes_bed$region <- paste(paste0("chr",genes_bed$seqnames),genes_bed$start,genes_bed$end,sep="_")     

        cnaRatio_infercnv <- cnaRatio_infercnv[genes_bed$gene_name,,drop=F]
        cnaRatio_infercnv <- data.frame(chr=paste0("chr",genes_bed$seqnames),start=genes_bed$start,end=genes_bed$end,cnaRatio_infercnv=cnaRatio_infercnv)

        cnaRatio_infercnv_align <- align_Grange2bin(bed.query=bed.ref,bed.subject=cnaRatio_infercnv)
        cnaRatio_infercnv_align <- cnaRatio_infercnv_align[,c("chr","start","end","cnaRatio_infercnv")]
         
        #load inferCNV HMM results
        cnvFile <- "infercnv.20_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt"
        cna_infercnv_raw <- read.table(paste(inferCNV_path,cnvFile,sep="/"))
        colnames(cna_infercnv_raw) <- gsub("\\.","-",colnames(cna_infercnv_raw))
        cna_infercnv <- cna_infercnv_raw[,obscell]
        cna_infercnv <- cna_infercnv[genes_bed$gene_name,,drop=F]
        cna_infercnv_rename <- data.frame(row.names=genes_bed$region,chr=paste0("chr",genes_bed$seqnames),start=genes_bed$start,end=genes_bed$end,cna_infercnv)

        #scCnv_all[[sample]][["infercnv"]] <- cna_infercnv_rename
        Psets  <- seq(0.1,1,0.1)
        for(prop in Psets){
            cna_infercnv_bulk <- get_infercnv_bulk(cna_infercnv_rename,bed.ref,cnvCellProp.min=prop)
            cna_infercnv_bulk$infercnv_ratio <- cnaRatio_infercnv_align$cnaRatio_infercnv
            colnames(cna_infercnv_bulk) <- c("chr","start","end","infercnv_CN","infercnv_ratio")
            rownames(cna_infercnv_bulk) <- paste(cna_infercnv_bulk$chr,cna_infercnv_bulk$start,cna_infercnv_bulk$end,sep="_")    
            write.csv(cna_infercnv_bulk,paste0(outdir_sample,"/",sample,"_infercnv_100kb_",prop,".csv"))
        }

    }

    ####----------------####
    #(6) Combine all results
    ####-----------------####
    dna <- read.csv(paste0(outdir_sample,"/",sample,"_DNA_100kb.csv"),header=TRUE,row.names=1)
    CNdf <- dna
    TeaCNV <- read.csv(paste0(outdir_sample,"/",sample,"_TeaCNV_100kb.csv"),header=TRUE,row.names=1)
    CNdf <- cbind(CNdf,TeaCNV[,4:ncol(TeaCNV)])
    for(prop in Psets){
        epiAneu <- read.csv(paste0(outdir_sample,"/",sample,"_epiAneufinder_100kb_",prop,".csv"),header=TRUE,row.names=1)
        CNdf_p <- CNdf
        CNdf_p <- cbind(CNdf_p,epiAneu[,4:ncol(epiAneu)])

        copyscAT <- read.csv(paste0(outdir_sample,"/",sample,"_copyscAT_100kb_",prop,".csv"),header=TRUE,row.names=1)
        CNdf_p <- cbind(CNdf_p,copyscAT[,4:ncol(copyscAT)])
       
       if(sample != "ccRCC1"){
            infercnv <- read.csv(paste0(outdir_sample,"/",sample,"_infercnv_100kb_",prop,".csv"),header=TRUE,row.names=1)
            CNdf_p <- cbind(CNdf_p,infercnv[,4:ncol(infercnv)])
            row_na_count <- apply(CNdf_p[,c("dna_CNA","TeaCNV.1","epiAneuCN","copyscAT.CN","infercnv_CN")], 1, function(x) sum(is.na(x)))
            CNdf_filt <- CNdf_p[row_na_count < 5, ] #remove rows with more than 5 missing values

        }else{
            row_na_count <- apply(CNdf_p[,c("dna_CNA","TeaCNV.1","epiAneuCN","copyscAT.CN")], 1, function(x) sum(is.na(x)))
            CNdf_filt <- CNdf_p[row_na_count < 4, ] #remove rows with more than 4 missing values

        }    
        write.csv(CNdf_filt, file = paste0(outdir_sample,"/",sample,"_CNdf_clean_clone_",prop,".csv"),row.names=FALSE)
   
    }

    rm(CNdf,dna)
}



### combine cnv-ratio for each method
###Fig.3 histogram
#select the DNA CNV regions
suppressPackageStartupMessages({
    library(ggpubr)
    library(ggpmisc)
    library(ggplot2)
    library(colorspace)
})

methods_set=c('TeaCNV', 'epiAneufinder', 'copyscAT','inferCNV')
cols_fill <- c("#652884","#CC5B45","#6A8EC9","#F5A216","#458A74")
names(cols_fill) <- c("WGS",methods_set)


prop=0.1
for(sam in sampleID){
    print(sam)
    CNdf_all_sam <- read.csv(paste0(outdir,"/",sam,"/",sam,"_CNdf_clean_clone_",prop,".csv"))
    write.csv(CNdf_all_sam,file = paste0(outdir,"/",sam,"_CNdf_clean_clone.csv"),quote = F,row.names = F)

   # CNdf_all_sam <- read.csv(paste0(outdir,"/",sam,"_CNdf_clean_clone.csv"))
    #convert the CNAs dataframe to numeric values
    CNdf_all_sam$epiAneuCN <- CNdf_all_sam$epiAneuCN+1 ##convert epiAneu score to CN leve

    CNdf_all_sam$dna_CNA[CNdf_all_sam$dna_CNA=="amp" & CNdf_all_sam$dna_ratio >1.5] <- 4
    CNdf_all_sam[CNdf_all_sam=='amp'] <- 3
    CNdf_all_sam[CNdf_all_sam=='neu'] <- 2
    CNdf_all_sam[CNdf_all_sam=='del'] <- 1
    CNdf_all_sam[is.na(CNdf_all_sam)] <- 0
    CNdf_all_sam$WGS <- CNdf_all_sam$dna_CNA
    if(sam == "ccRCC3"){
        base_ratio <- mean(CNdf_all_sam$dna_ratio[CNdf_all_sam$dna_CNA==2])
        cn3_ratio <- mean(CNdf_all_sam$dna_ratio[CNdf_all_sam$dna_CNA==3])
        delt <- cn3_ratio-base_ratio
        CNdf_all_sam$WGS[CNdf_all_sam$WGS==1 &CNdf_all_sam$dna_ratio>0.75 ] <- round((CNdf_all_sam$dna_ratio[CNdf_all_sam$WGS==1 &CNdf_all_sam$dna_ratio>0.75 ]-base_ratio)/delt+2,1)
    }


    CNdf_all_sam$chr <- gsub("chr", "", CNdf_all_sam$chr)  
      CNdf_all_sam <- CNdf_all_sam %>%
        dplyr::filter(chr %in% c(1:22)) %>%
        dplyr::filter(as.numeric(as.character(dna_CNA)) > 0) %>%
        as.data.frame()
    CNdf_all_sam <- droplevels(CNdf_all_sam)

    if(sam != "ccRCC1"){
        CNdf_all_sam <- CNdf_all_sam %>%
         dplyr::filter( as.numeric(as.character(infercnv_ratio))>0) %>%
            as.data.frame()
    }   
    CNdf_all_sam <- droplevels(CNdf_all_sam)
    CNdf_all_sam$dna_CNA <- as.numeric(as.character(CNdf_all_sam$dna_CNA))
    CNdf_all_sam$WGS <- as.numeric(as.character(CNdf_all_sam$WGS))

    ###(1) Plot log density histogram
    title_y <- "density (log)"
    ht_dna= hist_custom(CNdf_all_sam,x="WGS",log_y=TRUE,cols_fill=cols_fill[1],
    color_hist=TRUE,color_col="dna_CNA",
    title_x="WGS",title_y=title_y,binwidth=0.1,xlim=c(0.5,4.5))


    df_teaCNV_long <- CNdf_all_sam %>%
        pivot_longer(cols = starts_with("TeaCNV_mean"),  
               names_to = "clone",     
               values_to = "TeaCNV_mean") %>% as.data.frame()
    ht1 <- df_teaCNV_long[as.numeric(as.character(df_teaCNV_long[,"TeaCNV_mean"]))>0,,drop=FALSE] %>%
           hist_custom(x="TeaCNV_mean",log_y=TRUE,scale_factor=1.25,cols_fill=cols_fill[["TeaCNV"]],
           color_hist=TRUE,color_col="dna_CNA",
           title_x="TeaCNV",title_y=title_y,binwidth=0.1,xlim=c(0.5,4.5)) 

    pldata=  CNdf_all_sam[as.numeric(as.character(CNdf_all_sam[,"epiAneuRatio"]))>0,,drop=FALSE]
    colnm="epiAneuRatio"
    CNest <- pldata %>%
        group_by(epiAneuCN) %>%
        summarise(epiAneu_CNratio=mean(epiAneuRatio))%>%
        arrange(epiAneu_CNratio)
    delt <- CNest$epiAneu_CNratio[2]-CNest$epiAneu_CNratio[1]
    lim_x_max <- 1+2.5*delt;lim_x_min <- 1-1.5*delt   
    binWidth <- (lim_x_max - lim_x_min)*0.1/4
    # # lim_x_max <- round(max(as.numeric(as.character(pldata[,colnm])), na.rm = TRUE),2)
    # # lim_x_min <- min(as.numeric(as.character(pldata[,colnm])))
    # # lim_x_max <- 1.4;lim_x_min <- 0.75
    ht2 <- pldata %>%
        hist_custom(x="epiAneuRatio",log_y=TRUE,cols_fill=cols_fill[["epiAneufinder"]],
        color_hist=TRUE,color_col="dna_CNA",
        title_x="epiAneufinder",title_y=title_y,binwidth=binWidth,xlim=c(lim_x_min,lim_x_max),adjust = 3)

    colnm="copyscAT_ratio"
    pldata=  CNdf_all_sam[as.numeric(as.character(CNdf_all_sam[,colnm]))>0,,drop=FALSE]
    ht3 <- pldata %>%
        hist_custom(x="copyscAT_ratio",log_y=TRUE,cols_fill=cols_fill[["copyscAT"]],
        color_hist=TRUE,color_col="dna_CNA",
        binwidth = 0.1,xlim=c(0.5,4.5),adjust = 4,
        title_x="copyscAT",title_y=title_y)


    if(sam != "ccRCC1"){
        colnm="infercnv_ratio"
        pldata=  CNdf_all_sam[as.numeric(as.character(CNdf_all_sam[,colnm]))>0,,drop=FALSE]
        CNest <- pldata %>%
            group_by(infercnv_CN) %>%
            summarise(infercnv_CNratio=mean(infercnv_ratio))%>%
            arrange(infercnv_CNratio)
        delt <- CNest$infercnv_CNratio[2]-CNest$infercnv_CNratio[1]
        lim_x_max <- 1+2.5*delt;lim_x_min <- 1-1.5*delt   
        binWidth <- (lim_x_max - lim_x_min)*0.1/4
        # lim_x_max <- round(max(as.numeric(as.character(pldata[,colnm])), na.rm = TRUE),2)
        # lim_x_min <- min(as.numeric(as.character(pldata[,colnm])))
        ht4 <- pldata %>%
            hist_custom(x="infercnv_ratio",log_y=TRUE,cols_fill=cols_fill[["inferCNV"]],
            color_hist=TRUE,color_col="dna_CNA",
            binwidth = 0.005,xlim=c(lim_x_min,lim_x_max),adjust = 4,
            title_x="inferCNV",title_y=title_y)

        ht_com <- (ht_dna + ht1 + ht2+ht3+ht4)+plot_layout(ncol = 1, heights = c(1,1,1,1,1))
        height = 7
    }else{
        ht_com <- (ht_dna + ht1 + ht2+ht3)+plot_layout(ncol = 1,heights = c(1, 1,1,1))
        height = 7*4/5
    }
    ggsave(paste0(outdir,'/Hist_',sam,'-log.pdf'),ht_com,width = 3.5,height = height)

}



###Fig. 3 Segmentation plots
setwd(work_dir)
for(sam in sampleID){
    outdir_clt <- paste0(work_dir,'/',sam)

    CNVresFile_sc <- paste0("/data/final.CNVres_",sam,".rds")
    outres <- readRDS(CNVresFile_sc)

    #TeaCNV: clonal-CNV
    ##Fig.3b top panel
      #plot combined segmentation and CNV
    pl <- c() 
    pl_data <- c()
    ymax <- quantile(unlist(lapply(outres$clonalest,function(x){x$seg.dat$integerCN})),0.99)
    color_r <- cols_Palette[1:length(outres$clonalest)]
    names(color_r) <- sort(names(outres$clonalest))
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r
    for(cluster in sort(names(outres$clonalest))){
        clonal_res <- outres$clonalest[[cluster]]$input_BinSegRatio
        clonal_res <- clonal_res[,!grepl("relativeCN|integerCN",colnames(clonal_res))]
        integerCNV <- outres$clonalest[[cluster]]$seg.dat
        clonal_res <- merge(clonal_res,integerCNV[,c("relativeCN", "integerCN","segName")],by="segName",all.x=T)
         pl_scRatio <- seg_plot(clonal_res,name.data="",
                        value.bin="integerCN",
                        value.segment="integerCN",
                        ylab="CNAs",
                        plot_dots = F,
                        ylim=c(0,ymax),
                        add_yline = 2,
                        plot_seg = T,
                        seg_col = left_anno_cols[["clone"]][cluster],
                        plotDir=outdir_clt,outPlot=F,
                        color_hist=T,plot_hist=T,
                        color_hist_gradient=F,
                        plot_colors=left_anno_cols[["clone"]][cluster],
                        device = "pdf",disconnected=F)
       p1 <- pl_scRatio$p1
       pl[[cluster]] <- p1
       pl_data[[cluster]] <- pl_scRatio$p1$data
    }
    p1 <- pl[[1]]
    if(length(pl)>1){
        for(i in 2:length(pl_data)){
            cluster <- names(pl_data)[i]
            a <- pl_data[[cluster]]
            dat2<- data.frame(x=a$segEnd[1:(length(a$segEnd)-1)], y=a$value.segment[1:(length(a$segEnd)-1)],
                              xend=a$segStart[2:length(a$segEnd)], yend=a$value.segment[2:length(a$segEnd)])
            p1 <- p1 +
                geom_segment(
                data = dat2,
                mapping = aes(x=x, y=y,
                                xend=xend, yend=yend), 
                colour=left_anno_cols[["clone"]][cluster],na.rm =T,size = 1,
                inherit.aes = FALSE)
        }
    }
    ggsave(paste0(work_dir,"/segPlot_TeaCNV_",sam,".pdf"),p1,height = 2,width=10)

    ##Supplementary Fig. 6
    TeaCNV::plot_combine_seg(outres$clonalest,ylim=NULL,
        outplot_name=paste0("clonalCNV_final_",sam,".pdf"),show_dots=FALSE,
        outdir=outdir_clt)


    ###epiAneufinder
    cna_epiAneu <- read.table(paste0(data_path,"/epiAneufinder.results_table.",sam,".tsv"))
    rownames(cna_epiAneu)<- paste(cna_epiAneu$seq,cna_epiAneu$strat,cna_epiAneu$end,sep="_")
    cnvres <- cna_epiAneu[,4:ncol(cna_epiAneu)]
    #color setting
    max.value <- ceiling(quantile(cnvres,0.98,na.rm=TRUE))
    custom_colors <- list()
    if(max.value>=4){max.value=4}
    colorss <- colorRampPalette(c("#2A72B2","#C1DDF4", "grey95","#eae2b7","#ffba08","#FF6600","#d00000","#9d0208","#6a040f"))(9)
    colors <-colorRamp2(c(0,0.5,1,1.5,2,2.5,3,3.5,max.value),colorss)
    if(max.value<4 &max.value>1){
      colorss <- colorss[1:(2+(2*max.value-1))]
      colors <-colorRamp2(c(0,0.5,seq(1,max.value,by=0.5)),colorss)
    }else if(max.value==1){
        colorss <- colorss[1:3]
        colors <-colorRamp2(c(0,0.5,1),colorss)
    }

    custom_colors$colors <- colors
    custom_colors$color_bars <- "continuous"
    custom_colors$titles <- "Ratio"
    custom_colors$at_brk <- c(0:max.value)
    custom_colors$label_brk <- c(as.character(c(0:(max.value-1))),paste0(max.value,"+"))

    p_scRatio <- heatmap4peakMt(mat=cnvres,
                            meta_info=NULL,
                            sep_by="_",
                            outdir= outdir,value.type="ratio",
                            legend_titles="gain/loss",
                            clust_rows=T,clustering_method_rows = "ward.D2",
                            show_legend_row = T,
                            fileout_name=paste0("heatmap_epiAneufinder_",clt),
                            width=10,height=5,
                            custom_colors = custom_colors,
                            device="pdf")

    #segment plot
    cna_epiAneu_mean <- rowMeans(as.matrix(cnvres),na.rm=TRUE)
    cna_epiAneu_bulk = data.frame(cna_epiAneu[,1:3],epiAneu=cna_epiAneu_mean)
    colnames(cna_epiAneu_bulk)[1] <- c("seq")
    if (!is.data.table(cna_epiAneu_bulk)) {
      cna_epiAneu_bulk <- as.data.table(cna_epiAneu_bulk)
    }
    cna_epiAneu_bulk_df <- as.data.frame(cna_epiAneu_bulk[,c("epiAneu"),drop=F])
    rownames(cna_epiAneu_bulk_df) <- paste(cna_epiAneu_bulk$seq,cna_epiAneu_bulk$start,cna_epiAneu_bulk$end,sep="_")
    cluster_res <- segment4df(cna_epiAneu_bulk_df,cytoBand =cytoBandFile,doFilt=F,outPlot=F,seg.method="gaussian")
    df.seg <- cluster_res$seg_score_segLevel[[1]]
    df.seg <- df.seg[,c("segName","Num_bins")]
    seg_df <- cluster_res$seg_score_binLevel[[1]]
    seg_df <-dplyr::left_join(seg_df,df.seg,by="segName")
    write.csv(seg_df,file.path(work_dir,paste0("epiAneufinder_BulkSegResult_",sam,".csv")))

    plot_ylim <- round(quantile(seg_df$SegMean,c(0,1)),2)
    plot_ylim <- c(max(0,(0.8*plot_ylim[1])),plot_ylim[2]*1.1)
    p1 <- seg_plot(seg_df,name.data=clt,
                genome="hg38", add_yline=NULL,
                    value.bin="binRatio",
                    value.segment="SegMean",
                color_dot =FALSE,plotDir=outdir,ylab = "Ratio",outPlot=FALSE,ylim = plot_ylim,color_hist_gradient=F)
    p2 <- p1$ggarranged_p
    ggsave(paste0(outdir, "/segPlot_epiAneufinder_",sam,".pdf"),p2, width=9.5, height=2,device = pdf,bg="white")



    
    ###copyscAT
    cnvRes <- read.csv(paste0(data_path,"/CopyscAT.cnv_scores.",sam,".csv"),row.names=1)
    cna_copyscAT <- t(cnvRes)

    chr <- sapply(strsplit(rownames(cna_copyscAT),"p|q"),"[",1)
    arm <- ifelse(grepl("p",rownames(cna_copyscAT)),"p","q")
    start <- as.numeric(seg_pos$start[match(paste0(chr,"_",arm),paste0(seg_pos$chr,"_",seg_pos$arm2))])
    end <- as.numeric(seg_pos$end[match(paste0(chr,"_",arm),paste0(seg_pos$chr,"_",seg_pos$arm2))])

    chr_mum <- gsub("chr","",chr)
    chr_mum <- gsub("X|x","23",chr_mum)
    chr_mum <- gsub("Y|y","24",chr_mum)
    cna_copyscAT <- data.frame(chr_mum,chr,start,end,cna_copyscAT)
    cna_copyscAT <- cna_copyscAT %>%filter(chr_mum %in% c(1:22))%>%as.data.frame()

    cna_copyscAT <- merge(cna_copyscAT,seg_pos[,c("chr","start","end")],all=TRUE)
    cna_copyscAT$chr_mum <- gsub("chr","",cna_copyscAT$chr)
    cna_copyscAT <- cna_copyscAT[order(as.numeric(as.character(cna_copyscAT$chr_mum)),as.numeric(cna_copyscAT$start)),]

    cna_copyscAT <- cna_copyscAT[,!grepl("chr_mum",colnames(cna_copyscAT))]
    rownames(cna_copyscAT) <- paste(cna_copyscAT$chr,cna_copyscAT$start,cna_copyscAT$end,sep="_")

    cna_copyscAT <- cna_copyscAT[,4:ncol(cna_copyscAT)]
    cna_copyscAT <- cna_copyscAT[!grepl("chrX|chrY",rownames(cna_copyscAT)),]
    cna_copyscAT[cna_copyscAT<0] <- NA

    cna_copyscAT0 <- cna_copyscAT
    mt <- cna_copyscAT0

    #color setting
    max.value <- ceiling(quantile(mt,0.98,na.rm=TRUE))
    custom_colors <- list()
    if(max.value>=5){max.value=5}
    colorss <- colorRampPalette(c("#2A72B2","#C1DDF4", "grey95","#eae2b7","#ffba08","#FF6600","#d00000","#9d0208","#6a040f"))(9)
    colors <-colorRamp2(c(1,1.5,2,2.5,3,3.5,4,4.5,max.value),colorss)
    if(max.value<5 &max.value>2){
        colorss <- colorss[1:(2*max.value-1)]
        colors <-colorRamp2(c(seq(1,max.value,by=0.5)),colorss)
    }else if(max.value<=3 &max.value>1){
        x.center <- 2
        quantiles = quantile(mt[mt != x.center], c(0.01,0.25,0.75, 0.99),na.rm =TRUE)
        at_brk <- c(quantiles[1],quantiles[2],x.center,quantiles[3],quantiles[4])
        colorss <- colorRampPalette(colors = c("#2A72B2","#C1DDF4", "grey95","#eae2b7","#ffba08"))(length(at_brk))
        colors <- circlize::colorRamp2(at_brk,colorss)
    }else if(max.value==1){
        colorss <- colorss[1:3]
        colors <-colorRamp2(c(0,0.5,1),colorss)
    }

    custom_colors$colors <- colors
    custom_colors$color_bars <- "continuous"
    custom_colors$titles <- "Ratio"
    custom_colors$at_brk <- c(1:max.value)
    custom_colors$label_brk <- c(as.character(c(1:(max.value-1))),paste0(max.value,"+"))

    p_scRatio <- heatmap4peakMt(mat=cna_copyscAT0,
                               meta_info=NULL,
                               sep_by="_",
                               outdir= outdir_clt,value.type="ratio",
                               legend_titles="cnv scores",
                               clust_rows=T,clustering_method_rows = "ward.D2",
                              show_legend_row = T,
                               fileout_name=paste0("heatmap_copyscAT_",clt),
                               custom_colors = custom_colors,
                               width=10,height=5,device="pdf")

    #segment plot
    cna_mean <- rowMeans(as.matrix(cna_copyscAT),na.rm=TRUE)
    cna_bulk = data.frame(CN=cna_mean)

    cluster_res <- segment4df(cna_bulk,cytoBand =cytoBandFile,doFilt=F,outPlot=F,seg.method="gaussian",rmNA=F)
    df.seg <- cluster_res$seg_score_segLevel[[1]]
    df.seg <- df.seg[,c("segName","Num_bins")]
    seg_df <- cluster_res$seg_score_binLevel[[1]]
    seg_df <-dplyr::left_join(seg_df,df.seg,by="segName")
     write.csv(seg_df,file.path(outdir_clt,paste0("copyscAT_BulkSegResult_",clt,".csv")))

    plot_ylim <- round(quantile(seg_df$SegMean,c(0,1),na.rm=T),2)
    plot_ylim <- c(max(0,(0.8*plot_ylim[1])),plot_ylim[2]*1.1)
    p1 <- seg_plot(seg_df,name.data=clt,
                  genome="hg38", add_yline=NULL,
                      value.bin="binRatio",
                      value.segment="SegMean",
                  color_dot =FALSE,plotDir=outdir_clt,ylab = "cnv score",outPlot=FALSE,ylim = plot_ylim,
                  plot_seg=T,
                  color_hist_gradient=F)
    p2 <- p1$ggarranged_p
    ggsave(paste0(outdir_clt, "/segPlot_copyscAT_",clt,".pdf"),p2, width=9.5, height=2,device = pdf,bg="white")

}


    ###infercnv
cytoBandFile <- file.path(script_path, "data", "cytoBand_hg38.tsv")
gtf_file <- "~/Library/Mobile Documents/com~apple~CloudDocs/reference/genes.gtf.gz"
annotations <- import(gtf_file)
annotations <- annotations[annotations$type=="gene",]
gene_order_file="~/Library/Mobile Documents/com~apple~CloudDocs/reference/hg38_filtered_order_gene.txt"
cnvFile <- "run.final.infercnv_obj" 
cnvdir <- data_path
infercnv_obj = readRDS(paste(cnvdir,cnvFile,sep="/"))
mtx_cnv <- infercnv_obj@expr.data
mt_obs <- as.matrix(mtx_cnv[,unlist(infercnv_obj@observation_grouped_cell_indices)])
color.breaks.file=paste0(cnvdir,"/infercnv.heatmap_thresholds.txt")
bk <- read.table(color.breaks.file, header=F)[,1]
cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")


for(sam in c("ccRCC3","ccRCC4")){
    outdir_clt <- paste0(work_dir,'/',sam)

    cellinfo <- read.csv(paste0(data_path, '/cell_info_subClust_',sam,'.csv'),row.names=1)

    cellMeta <- cellinfo[cellinfo$group %in%"observed",,drop=F]
    cellMeta$subCluster <- factor(cellMeta$subCluster,levels=sort(unique(as.numeric(cellMeta$subCluster))))
    cellMeta$clone_merged <- factor(cellMeta$clone_merged,levels=sort(unique(as.numeric(cellMeta$clone_merged))))
    cellMeta <- cellMeta[order(cellMeta$clone_merged,cellMeta$subCluster),,drop=F]
    colnames(cellMeta) <- gsub("clone_merged","clone",colnames(cellMeta))
    left_anno_cols<- c()
    color_r2 <- cols_Palette[1:length(unique(cellMeta$clone))]
    names(color_r2) <- sort(unique(cellMeta$clone))
    left_anno_cols[['clone']] <- color_r2
    cellMeta2 <- cellMeta[,c("clone"),drop=F]

    p_scRatio <- replotInfercnv(expr.mat=mt_obs[,rownames(cellMeta2)],
                            gene.order.file=gene_order_file,
                             metadata=cellMeta2,
                             plotDir= cnvdir,
                             clust_rows=F,
                            show_legend_row = T,
                             names=paste0("heatmap_infercnv_",sam),
                             left_anno_col=left_anno_cols,
                             width=10,height=5)

    #segment plot
    cna_infercnv <- mt_obs[,rownames(cellMeta2)]
    genes <- rownames(cna_infercnv)
    gene_ranges <- annotations[annotations$gene_name %in% genes,]
    genes_bed <- data.frame(seqnames=seqnames(gene_ranges),start=start(gene_ranges),end=end(gene_ranges),gene_name=gene_ranges$gene_name) %>%
        arrange(as.numeric(seqnames),as.numeric(start)) %>%
        as.data.frame()
    genes_bed$region <- paste(genes_bed$seqnames,genes_bed$start,genes_bed$end,sep="_")     
    cna_infercnv <- cna_infercnv[genes_bed$gene_name,,drop=F]
    cna_infercnv_rename <- data.frame(row.names=genes_bed$region,cna_infercnv)

    cna_mean <- rowMeans(as.matrix(cna_infercnv_rename),na.rm=TRUE)
    cna_bulk = data.frame(CN=cna_mean)

    cluster_res <- segment4df(cna_bulk,cytoBand =cytoBandFile,doFilt=F,outPlot=F,seg.method="gaussian",rmNA=F)
    df.seg <- cluster_res$seg_score_segLevel[[1]]
    df.seg <- df.seg[,c("segName","Num_bins")]
    seg_df <- cluster_res$seg_score_binLevel[[1]]
    seg_df <-dplyr::left_join(seg_df,df.seg,by="segName")
    write.csv(seg_df,file.path(outdir_clt,paste0("infercnv_BulkSegResult_",sam,".csv")))



    plot_ylim <- round(quantile(seg_df$SegMean,c(0,1),na.rm=T),2)
    plot_ylim <- c(max(0,(0.8*plot_ylim[1])),plot_ylim[2]*1.2)
    p1 <- seg_plot(seg_df,name.data=clt,
                  genome="hg38", add_yline=NULL,
                      value.bin="binRatio",
                      value.segment="SegMean",
                  color_dot =FALSE,plotDir=outdir_clt,ylab = "cnv score",outPlot=FALSE,ylim = plot_ylim,
                  plot_seg=F,
                  color_hist_gradient=F)
    p2 <- p1$ggarranged_p
    ggsave(paste0(outdir_clt, "/segPlot_infercnv_",sam,".pdf"),p2, width=9.5, height=2,device = pdf,bg="white")

}

