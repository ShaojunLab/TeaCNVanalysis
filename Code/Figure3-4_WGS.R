#Figure3,4
suppressMessages({
    library(ggplot2)
    library(gridExtra)
    library(ggpubr)
    library(dplyr)
    library(TeaCNV)
})

work_dir <- '/code/'
setwd(work_dir)
outdir <- '/results'
cytoBandFile <- file.path("/data/cytoBand_hg38.tsv")
normal_cnv_file <- file.path("/data/2021415KP_WGS_bincounts.csv")
normal_cnv_df <- read.csv(normal_cnv_file,header=T,stringsAsFactors=F,row.names=1)

samples_tumor <- c(ccRCC4="ccRCC4",ccRCC3="ccRCC3",ccRCC1="ccRCC1")
###Fig.3a,3f,4a: WGS
for(sample in samples_tumor){
    cnv_file <- file.path(work_dir,paste0(sample,"_WGS_bincounts.csv"))
    if(file.exists(cnv_file)){
        cnv_df <- read.csv(cnv_file,header=T,stringsAsFactors=F,row.names=1)
    }
    matCount <- cbind(normal_cnv_df,cnv_df[match(rownames(normal_cnv_df),rownames(cnv_df)),,drop=F])
    colnames(matCount)<- c("normal","tumor")
    #normalization
    matCount_norm <-  t(t(matCount)*mean(colSums(matCount))/colSums(matCount))
    matCount$ratio <-    matCount_norm[,"tumor"]/matCount_norm[,"normal"]
    matCount$logratio <- log2(matCount$ratio)
    #segmenation
    seg_res <- segment_by_chr_location(matCount[,"logratio",drop=F],cytoBand=cytoBandFile,
                                            acen_remove=TRUE,
                                            method="gaussian",
                                            segSize_min=5,
                                            penalty=1)

    matCount$changepoint <- seg_res$changepoint_matrix[,1]

    data_ls <- list(changepoint_matrix=matCount[,"changepoint",drop=F],data_matrix=matCount[,"ratio",drop=F],logratio_matrix=matCount[,"logratio",drop=F])
    segScore <- seg_score(data_ls,method="density")
    df.seg <- segScore$seg_score_segLevel[[1]]
    df.seg <- df.seg[,c("segName","Num_bins")]
    seg_df <- segScore$seg_score_binLevel[[1]]
    seg_df <-dplyr::left_join(seg_df,df.seg,by="segName")
    rownames(seg_df) <- seg_df$binID
    write.csv(seg_df,file.path(outdir,paste0(sample,"_SegResult.csv")))

    seg_df <- seg_df[!seg_df$Chromosome %in% c("chrX","chrY"),]

    parm_adj <- 2
    plot_ylim <- c(0.5,1.5)

    if(sample=="ccRCC1"){parm_adj=4}
    if(sample=="ccRCC4"){parm_adj=5}
    pp <- density(round(seg_df$SegMean,2), adjust = parm_adj)
    ploidy <- unlist(sapply(2:(length(pp$x)-1), function(j,pp){
    if (pp$y[j]>=pp$y[j-1]&pp$y[j]>=pp$y[j+1]){
    return(pp$x[j])
    }
    },pp))
    plot(pp)
    ploidy
    diff(ploidy)
    if(sample=="ccRCC3"){
    ploidy <- ploidy[ploidy>1.1 & ploidy<1.5]
    delt <- 0.28
    }else if(sample=="ccRCC1"){
    delt <- 0.19
    }else if(sample=="ccRCC4"){
    delt <- 0.29
    plot_ylim <- c(0.5,1.6)
    }
    baseRatio <- pp$x[which.max(pp$y)]
    cnRatio <- c(baseRatio-delt,baseRatio+delt)

  
    d_mvg <- 1-baseRatio
    seg_df$SegMean_adj <- seg_df$SegMean+d_mvg
    seg_df$binRatio_adj <- seg_df$binRatio+d_mvg
    cnRatio <- cnRatio+d_mvg

    p1 <- seg_plot(seg_df,name.data=sample,
                genome="hg38", add_yline=NULL,
                    value.bin="binRatio_adj",
                    value.segment="SegMean_adj",
                color_dot =FALSE,plotDir=outdir,ylab = "Ratio",outPlot=FALSE,ylim = plot_ylim)
    p2 <- p1$p1
    right_axis <- scale_y_continuous(
    name = "Ratio",
    sec.axis = sec_axis(~ .,breaks = c(cnRatio[1],1,cnRatio[2]),
                                    labels = c("1", "2", "3"),name = "CNA"))

    p1$p2
    p3 <- p2 + labs(x = "", y = "Ratio")+
    geom_hline(yintercept = cnRatio, linetype = "dashed", color = "grey", size = 0.5)+
    right_axis
    #p3
    ggsave(paste0(outdir, "/DNA_segPlot_",sample,".pdf"),p3, width=9.5, height=2,device = pdf,bg="white")

    p4=ggarrange(p3, p1$p2, align ="h",ncol = 2, nrow = 1,widths = c(15,4),heights=3 )
    ggsave(paste0(outdir, "/DNA_segPlot_",sample,"-hist.pdf"),p4, width=15, height=3,device = pdf,bg="white")

    ##Seglength of CNV 
    raio.hi=1.2;raio.lo=0.8
    if(sample=="ccRCC1"){
      delt <- 0.19
      raio.hi <- cnRatio[2]
      raio.lo<- cnRatio[1]
    }
    seg_dna <-  data.frame(seg_df[,c('Chromosome','Start','End','SegMean',"segName","length_seg","Num_bins")])
    #delt <- 0.28
    seg_dna = seg_dna %>%
                mutate(cnv_state = case_when(
                    SegMean > raio.hi ~ 'amp',
                    SegMean < raio.lo ~ 'del',
                    T ~ 'neu')
                )
    seg_dna$dna_CNA <- NA
    seg_dna$dna_CNA[seg_dna$cnv_state=="amp" & seg_dna$SegMean >1.5] <- 4
    seg_dna$dna_CNA[seg_dna$cnv_state=='amp'] <- 3
    seg_dna$dna_CNA[seg_dna$cnv_state=='neu'] <- 2
    seg_dna$dna_CNA[seg_dna$cnv_state=='del'] <- 1
    seg_dna$dna_CNA[is.na(seg_dna$cnv_state)] <- 0

    seg_dna_filt <- seg_dna %>% 
        dplyr::filter(Num_bins > 1)%>%
        as.data.frame()   
    seg_dna2 <- unique(seg_dna_filt[,c('Chromosome',"dna_CNA","segName",'SegMean')])
    seg_dna2$Start <- as.numeric(sapply(strsplit(seg_dna2$segName, "_"), "[", 2))
    seg_dna2$End <- as.numeric(sapply(strsplit(seg_dna2$segName, "_"), "[", 3))
    seg_dna2 <- unique(seg_dna2[,c('Chromosome','Start','End','dna_CNA','SegMean')])
    colnames(seg_dna2) <- c('Chromosome','Start','End','dna_CNA',"dna_ratio")


    merged_df <- seg_dna2 %>%
    # Create a grouping variable to indicate whether the integerCN of adjacent rows is the same
    dplyr::group_by(Chromosome)%>%
    dplyr::mutate(newSeg = cumsum(dna_CNA != lag(dna_CNA, default = first(dna_CNA)))) %>%
    ungroup()%>%
    dplyr::group_by(Chromosome,newSeg, dna_CNA) %>%
    summarise(
      start = min(Start),    # 
      end = max(End)          # 
    ) %>%
    ungroup()%>%
    dplyr::distinct(Chromosome, start, end, .keep_all = TRUE)   %>%
    mutate(segLen=(end-start)/1e6)%>%
    as.data.frame()
    merged_df <- merged_df %>%
    mutate(chr =as.numeric(gsub("chr", "", merged_df$Chromosome))) %>%
    arrange(chr, start, end)%>%as.data.frame()
    
    write.csv(merged_df,file.path(outdir,paste0(sample,"_SegCNV_res.csv")))

}

