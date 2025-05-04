
#' @param dnaFilePath Path to the CNV analysis results output by cnvEsti.r for WGS data.
#' @param bed.ref A bed file or a data frame of query regions
get_dna_bulk = function(dnaFilePath,bed.ref,length_lim=500000,raio.hi=1.2,raio.lo=0.8){
    seg_dna <- read.csv(dnaFilePath,header=T,row.names=1)
    seg_dna <-  data.frame(seg_dna[,c('Chromosome','Start','End','SegMean',"segName","length_seg")])
    #delt <- 0.28
    seg_dna = seg_dna %>%
                mutate(cnv_state = case_when(
                    SegMean > raio.hi ~ 'amp',
                    SegMean < raio.lo ~ 'del',
                    T ~ 'neu')
                )
    seg_dna_filt <- seg_dna %>% 
        dplyr::filter(length_seg > length_lim)%>%
        as.data.frame()   
    seg_dna2 <- unique(seg_dna_filt[,c('Chromosome',"cnv_state","segName",'SegMean')])
    seg_dna2$Start <- sapply(strsplit(seg_dna2$segName, "_"), "[", 2)
    seg_dna2$End <- sapply(strsplit(seg_dna2$segName, "_"), "[", 3)
    seg_dna2 <- unique(seg_dna2[,c('Chromosome','Start','End','cnv_state','SegMean')])
    colnames(seg_dna2) <- c("chr","strat","end","dna_CNA","dna_ratio")
    CNdf <- align_Grange2bin(bed.query=bed.ref,bed.subject=seg_dna2)
    CNdf <-CNdf[,c("chr","start","end","dna_CNA","dna_ratio")]
    return(CNdf)
}

#' @param cna_epiAneu data.frame of epiAneufinder CNA results. The first three columns are chr, start, end, and the remaining columns are the CNA values: 0(loss), 1(neutral), 2(gain).
#' @return A data frame of the same order as the input data frame, with the corresponding bulk-level epiAneufinder CNA values.
epiAneufinder_AssignBulkCN = function(cna_epiAneu,cytoBandFile,Ratio_cutoff_del=0.95,Ratio_cutoff_amp=1.05){
    suppressMessages({
        library(epiAneufinder)
        library(data.table)
        library(parallel)
    })
    cna_epiAneu_mean <- rowMeans(as.matrix(cna_epiAneu[,4:ncol(cna_epiAneu)]),na.rm=TRUE)
    cna_epiAneu_bulk = data.frame(cna_epiAneu[,1:3],epiAneu=cna_epiAneu_mean)
    colnames(cna_epiAneu_bulk)[1] <- c("seq")
    if (!is.data.table(cna_epiAneu_bulk)) {
        cna_epiAneu_bulk <- as.data.table(cna_epiAneu_bulk)
    }
    cna_epiAneu_bulk_df <- as.data.frame(cna_epiAneu_bulk[,c("epiAneu"),drop=F])
    rownames(cna_epiAneu_bulk_df) <- paste(cna_epiAneu_bulk$seq,cna_epiAneu_bulk$start,cna_epiAneu_bulk$end,sep="_")
    cluster_res <- segment4df(cna_epiAneu_bulk_df,cytoBand =cytoBandFile,doFilt=F,outPlot=F,seg.method="gaussian")
    cluster <- cluster_res$seg_score_binLevel$epiAneu$segID

   ## Calculate mean read counts per cluster
    seq_data=cna_epiAneu_bulk$epiAneu
    counts.normal <- seq_data / mean(seq_data)
    counts.normal[counts.normal< 0] <- 0
    qus_global <- quantile(seq_data, c(0.01, 0.98))
    cnmean <- sapply(split(counts.normal,cluster), function(x) {
        qus <- quantile(x, c(0.1, 0.9))
        y <- x[x >= qus[1] & x <= qus[2] & x >= qus_global[1] & x <= qus_global[2]]
        if(sum(y) == 0 | length(y) == 0)
        y <- x
        mean(y)
    })
    cnmean_significance <- dplyr::between(scale(cnmean), -1, 1)
    cnmean.scaled <- cnmean/mean(cnmean[cnmean_significance])
    CN.states <- cnmean.scaled[as.character(cluster)]
    aepiCN <- ifelse(CN.states > Ratio_cutoff_amp, 2, ifelse(CN.states < Ratio_cutoff_del, 0, 1))
    # aepiCN <- assign_gainloss(seq_data=cna_epiAneu_bulk$epiAneu, cluster=cluster, uq=0.99, lq=0.01) #only work for single-cell
    cna_epiAneu_bulk$epiAneuCN <- aepiCN
    colnames(cna_epiAneu_bulk)[1:3] <- c("chr","start","end")
    return(cna_epiAneu_bulk)
}


#' @param cna_mtx data.frame of epiAneufinder CNA results. The first three columns are chr, start, end, and the remaining columns are the CNA values.
#' #' @return A data frame of the same order as the input data frame, with the corresponding bulk-level epiAneufinder CNA values.
get_epiAneufinder_bulk_v0 <- function(cna_mtx,bed.ref,cytoBandFile,Ratio_cutoff_del=0.95,Ratio_cutoff_amp=1.05){
    cna_epiAneu_bulk <- epiAneufinder_AssignBulkCN(cna_mtx,cytoBandFile,Ratio_cutoff_del,Ratio_cutoff_amp)
    cna_epiAneu_bulk <- as.data.frame(cna_epiAneu_bulk)
    cna_epiAneu_bulk_align <- align_Grange2bin(bed.query=bed.ref,bed.subject=cna_epiAneu_bulk)
    cna_epiAneu_bulk_align <- cna_epiAneu_bulk_align[,c("chr","start","end","epiAneu","epiAneuCN")]
    cna_epiAneu_bulk_align <- cna_epiAneu_bulk_align %>%
        mutate(epiAneuCNstate = case_when(
            epiAneuCN ==2 ~ 'amp',
            epiAneuCN ==0 ~ 'del',
            is.na(epiAneuCN) ~ NA,
            T ~ 'neu'
        ))
    return(cna_epiAneu_bulk_align)
}
#' @param cna_mtx data.frame of epiAneufinder CNA results.
#' The first three columns are chr, start, end, and the remaining columns are the CNA values (0=del, 1=neu, 2=amp).
#' #' @return A data frame of the same order as the input data frame, with the corresponding bulk-level epiAneufinder CNA values.
get_epiAneufinder_bulk <- function(cna_mtx,bed.ref,cnvCellProp.min=0.2){
    cna_MT <- as.matrix(cna_mtx[,4:ncol(cna_mtx)])
    cna_epiAneu_mean <- rowMeans(cna_MT,na.rm=TRUE)

    cnv_N_in_row <- apply(cna_MT, 1, function(row) table(factor(row, levels = 0:2)))
    cnv_N_in_row <- t(cnv_N_in_row)
    cnv_Prop_in_row <- cnv_N_in_row/ncol(cna_MT)
    delamp_Prop_in_row <- cnv_Prop_in_row[,c(1,3)]
    epiAneuCN <- apply(delamp_Prop_in_row,1,function(x,cnvCellProp.min){
        CNVlevel <-  c(0,2)[which.max(x)]
        if(x[which.max(x)] < cnvCellProp.min){CNVlevel <- 1}
        return(CNVlevel)
    },cnvCellProp.min)
    cna_epiAneu_bulk <- data.frame(cna_mtx[,1:3],epiAneu=cna_epiAneu_mean,as.data.frame(epiAneuCN))
    cna_epiAneu_bulk_align <- align_Grange2bin(bed.query=bed.ref,bed.subject=cna_epiAneu_bulk)
    cna_epiAneu_bulk_align <- cna_epiAneu_bulk_align[,c("chr","start","end","epiAneu","epiAneuCN")]
    # cna_epiAneu_bulk_align <- cna_epiAneu_bulk_align %>%
    #     mutate(epiAneuCNstate = case_when(
    #         epiAneuCN ==2 ~ 'amp',
    #         epiAneuCN ==0 ~ 'del',
    #         is.na(epiAneuCN) ~ NA,
    #         T ~ 'neu'
    #     ))
    return(cna_epiAneu_bulk_align)
}






#' @param cnaFilePath Path to the "cnv_scores.csv" file from copyscAT analysis. The values are the copy number scores for each segment (score: ～2=neutral).
#' @param cytoBandFile Path to the cytoBand file for the genome build. The first five columns are c("chr", "start", "end","arm","band")
#' @return A data frame of the same order as the input data frame, with the corresponding bulk-level epiAneufinder CNA values.
get_copyscAT_sc = function(cnaFilePath,cytoBandFile){
    hgpos <- read.table(cytoBandFile, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    colnames(hgpos) <- c("chr", "start", "end","arm","band")
    hgpos$arm2 <- substr(hgpos$arm, 1, 1)
    hgpos <- hgpos %>%
        mutate(arm2 = substr(arm, 1, 1)) %>%
        dplyr::filter(band != "acen")%>%
        group_by(chr,arm2) %>%
        mutate(seg_start = min(start), seg_end = max(end))
    seg_pos <- unique(hgpos[,c("chr","seg_start","seg_end","arm2")])

    cna_copyscAT <- read.csv(cnaFilePath,row.names=1)
    cna_copyscAT <- t(cna_copyscAT)
    
    chr <- sapply(strsplit(rownames(cna_copyscAT),"p|q"),"[",1)
    arm <- ifelse(grepl("p",rownames(cna_copyscAT)),"p","q")
    start <- format(seg_pos$seg_start[match(paste0(chr,"_",arm),paste0(seg_pos$chr,"_",seg_pos$arm2))], scientific = FALSE)
    end <- format(seg_pos$seg_end[match(paste0(chr,"_",arm),paste0(seg_pos$chr,"_",seg_pos$arm2))], scientific = FALSE)

    chr_mum <- gsub("chr","",chr)
    chr_mum <- gsub("X|x","23",chr_mum)
    chr_mum <- gsub("Y|y","24",chr_mum)
    cna_copyscAT[cna_copyscAT<0] <- NA  #remove invalid values
    cna_copyscAT <- data.frame(chr_mum,chr,start,end,cna_copyscAT)
    cna_copyscAT <- cna_copyscAT[order(as.numeric(cna_copyscAT$chr_mum),as.numeric(cna_copyscAT$start)),]
    cna_copyscAT <- cna_copyscAT[,2:ncol(cna_copyscAT)]

    return(cna_copyscAT)

}

#' @param cna_copyscAT data.frame of copyscAT CNA results. The first three columns are chr, start, end, and the remaining columns are the CNA values.
#' @return A data frame of the same order as the input data frame, with the corresponding bulk-level copyscAT CNA values.
get_copyscAT_bulk_v0 = function(cna_copyscAT,bed.ref){
    cna_copyscAT_mean <- rowMeans(as.matrix(cna_copyscAT[,4:ncol(cna_copyscAT)]),na.rm=TRUE)
    cna_copyscAT_bulk <- data.frame(cna_copyscAT[,1:3],copyscAT_ratio=cna_copyscAT_mean)

    cna_copyscAT_bulk = cna_copyscAT_bulk %>% 
            mutate(copyscAT.CN = case_when(
                copyscAT_ratio > 2.5 ~ 'amp',
                copyscAT_ratio < 1.5 ~ 'del',
                T ~ 'neu'
            ))
    cna_copyscAT_bulk_align <- align_Grange2bin(bed.query=bed.ref,bed.subject=cna_copyscAT_bulk)
    cna_copyscAT_bulk_align <- cna_copyscAT_bulk_align[,c("chr","start","end","copyscAT_ratio","copyscAT.CN")]

    return(cna_copyscAT_bulk_align)
}
#' @param cna_copyscAT data.frame of copyscAT CNA results. 
#' The first three columns are chr, start, end, and the remaining columns are the CNA values (score: 2=neutral).
#' @return A data frame of the same order as the input data frame, with the corresponding bulk-level copyscAT CNA values.

get_copyscAT_bulk = function(cna_copyscAT,bed.ref,cnvCellProp.min=0.2){
    cna_MT <- as.matrix(cna_copyscAT[,4:ncol(cna_copyscAT)])
    cna_copyscAT_mean <- rowMeans(cna_MT,na.rm=TRUE)
    cnv_N_in_row <- apply(cna_MT, 1, function(row){
        cn_discrete = ifelse(row> 2.5,3,ifelse(row<1.5,1,2))
        table(factor(cn_discrete, levels = 1:3))
    } )
    cnv_N_in_row <- t(cnv_N_in_row)
    cnv_Prop_in_row <- cnv_N_in_row/ncol(cna_MT)
    delamp_Prop_in_row <- cnv_Prop_in_row[,c(1,3)]
    copyscAT.CN <- apply(delamp_Prop_in_row,1,function(x,cnvCellProp.min){
        CNVlevel <-  c(1,3)[which.max(x)]
        if(x[which.max(x)] < cnvCellProp.min){CNVlevel <- 2}
        return(CNVlevel)
    },cnvCellProp.min)
    cna_copyscAT_bulk <- data.frame(cna_copyscAT[,1:3],copyscAT_ratio=cna_copyscAT_mean,as.data.frame(copyscAT.CN))

    cna_copyscAT_bulk_align <- align_Grange2bin(bed.query=bed.ref,bed.subject=cna_copyscAT_bulk)
    cna_copyscAT_bulk_align <- cna_copyscAT_bulk_align[,c("chr","start","end","copyscAT_ratio","copyscAT.CN")]
    return(cna_copyscAT_bulk_align)
}

#' @param cna_infercnv_mt data.frame of infercnv CNA results. 
#' The first three columns are chr, start, end, and the remaining columns are the CNA values (HMM CNV state: 0.5,neu=1,1.5,2...).
#' @return A data frame of the same order as the input data frame, with the corresponding bulk-level copyscAT CNA values.

get_infercnv_bulk = function(cna_infercnv_mt,bed.ref,cnvCellProp.min=0.2){
    cna_MT <- as.matrix(cna_infercnv_mt[,4:ncol(cna_infercnv_mt)])
    cnv_N_in_row <- apply(cna_MT, 1, function(row){
        cn_absolute = 2*row
        table(factor(cn_absolute, levels = 0:5))
    } )
    cnv_N_in_row <- t(cnv_N_in_row)
    cnv_Prop_in_row <- cnv_N_in_row/ncol(cna_MT)
    delamp_Prop_in_row <- cnv_Prop_in_row[,-3]
    infercnv_CN <- apply(delamp_Prop_in_row,1,function(x,cnvCellProp.min){
        CNVlevel <-  c(0:1,3:5)[which.max(x)]
        if(x[which.max(x)] < cnvCellProp.min){CNVlevel <- 2}
        return(CNVlevel)
    },cnvCellProp.min)
    cna_bulk <- data.frame(cna_infercnv_mt[,1:3],as.data.frame(infercnv_CN))

    cna_bulk_align <- align_Grange2bin(bed.query=bed.ref,bed.subject=cna_bulk)
    cna_bulk_align <- cna_bulk_align[,c("chr","start","end","infercnv_CN")]
    return(cna_bulk_align)
}


