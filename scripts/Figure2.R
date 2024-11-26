#Figure2
suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(ggsci)
  library(Signac)
  library(Seurat)
  library(TeaCNV)
 })

work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
setwd(work_dir)
source('./Code/Function/atac_visualizeFun.R')
#plotdir <- '~/Library/Mobile Documents/com~apple~CloudDocs/code/code_CNV/TeaCNVAnalysis'
plotdir <- './'
# setwd(plotdir)
color_celltype <- c(Endothelial="#CCFFCC",Epithelial="#C5A48A",Immne="#CCCCCC",Stromal="#CC99FF")
cols_Palette = c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780"),

cnv_plots_comb <- function(ID,outres,
	cols_Palette = c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780"),
	plotdir="./"
	){
	#(1)
	TeaCNV::plot_combine_seg(outres$clonalest,ylim=NULL,
	    outplot_name=paste0(ID,".clonalCNV_final.pdf"),show_dots=FALSE,
	    outdir=plotdir)
	#(2)
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
	clone_info <- data.frame(row.names=colnames(CNmt),clone=colnames(CNmt))
	color_r <- cols_Palette[1:length(unique(clone_info$clone))]
	    names(color_r) <- sort(unique(clone_info$clone))
	    left_anno_cols <- list()
	    left_anno_cols[["clone"]] <- color_r
	    height <- ifelse(ncol(CNmt)>2,0.35*(ncol(CNmt))+0.5,ifelse(ncol(CNmt)==1,1.2,1.5))
	    p_cloneCN <- heatmap4peakMt(mat=CNmt,
	                          meta_info=clone_info,
	                          sep_by="-",
	                          outdir= plotdir,
	                          value.type="CNV",
	                          clust_rows=F,
	                          show_legend_row = T,
	                          legend_titles="integer CN",
	                          fileout_name=paste0(ID,".heatmap_cloneCNA"),
	                          col_list=left_anno_cols,
	                          column.title = NULL,
	                          width=10,height=height,device="pdf") 

	#(3)
	plt_mt2 <- outres$cellbinRatio_raw
	plt_mt2 <- smooth_and_denoise(plt_mt2,window=60)
	left_anno_cols <- list()
	cellMeta2 <- outres$cellinfo[,c("clone"),drop=F]
	cellMeta2 <- na.omit(cellMeta2)
	color_r2 <- cols_Palette[1:length(unique(cellMeta2$clone))]
	names(color_r2) <- sort(unique(cellMeta2$clone))
	left_anno_cols[['clone']] <- color_r2
	p_scRatio <- heatmap4peakMt(mat=plt_mt2[,rownames(cellMeta2)],
	                         meta_info=cellMeta2,
	                         sep_by="-",
	                         outdir= plotdir,value.type="ratio",
	                         clust_rows=F,clustering_method_rows = "ward.D2",
	                        show_legend_row = T,
	                         fileout_name=paste0(ID,".heatmap_scRatioRaw"),
	                         col_list=left_anno_cols,
	                         width=10,height=5,device="pdf")
}
smooth_median_by_group <- function(df,k = 10,adjust=1) {
	df %>%
	  group_by(integerCN) %>%
	  mutate(
	    smoothed_binRatio = rollmean(binRatio, k = k, fill = NA, align = "right")  # 计算步长为k的均值
	  ) %>%
	  filter(!is.na(smoothed_binRatio)) %>% 
	  # mutate(
	  #   relativeCN_mean = median(smoothed_binRatio, na.rm = TRUE)  # 计算平滑后的中值
	  # )%>%   
	  mutate(
	     relativeCN_mean = {
	      d <- density(smoothed_binRatio,adjust=adjust)
	      max_density_idx <- which.max(d$y)
	      d$x[max_density_idx]
	    }
	  )%>%
	  as.data.frame()
}
DEseg <- function(outres){
	clonalest <- outres$clonalest
	res <- c()
	for(i in 1:length(clonalest)){
	  df_bin <- clonalest[[i]]$input_BinSegRatio
	  df_bin <- df_bin[,c("binID","Chromosome","Start","End","binRatio","length_bin","SegMean","segName")]
	  df_seg <- clonalest[[i]]$seg.dat
	  df_seg <- df_seg[,c("segName","relativeCN","integerCN")]
	  df_bin <- left_join(df_bin,df_seg,by=c("segName"))
	  df_bin$clone <- names(clonalest)[i]
	  df_bin$wei <-  df_bin$length_bin/sum(df_bin$length_bin)
	  res <- rbind(res,df_bin)
	}
	#select differential CNV region accross clones
	res_cnv <- unique(res[,c("Chromosome","segName","integerCN","clone")])
	res_cnv <- na.omit(res_cnv[res_cnv$integerCN != "2",])
	res_cnv$segstart <- as.numeric(sapply(strsplit(res_cnv$segName,"_|-|:"),"[",2))
	res_cnv$segend <- as.numeric(sapply(strsplit(res_cnv$segName,"_|-|:"),"[",3))
	filtered_res_cnv <- res_cnv %>%
	group_by(Chromosome, integerCN) %>%
	filter(n_distinct(clone) != length(unique(res_cnv$clone))) %>%
	ungroup()%>%
	as.data.frame()
	merged_df <- filtered_res_cnv %>%
	    dplyr::group_by(Chromosome)%>%
	    dplyr::mutate(newSeg = cumsum(integerCN != lag(integerCN, default = first(integerCN)) | segstart >= lag(segend + 2e6, default = first(segstart)))) %>%
	    ungroup()%>%
	    dplyr::group_by(Chromosome,newSeg, integerCN) %>%
	    summarise(
	      start = first(segstart),    
	      end = last(segend)          
	    ) %>%
	    ungroup()%>%
	    as.data.frame()
	merged_df <- merged_df[!duplicated(merged_df[, c("Chromosome", "integerCN","start","end")]), ]
return(list(res=res,merged_df=merged_df))
}
custom_theme <-
    theme_classic()+ 
    theme(plot.background=element_blank(),
          legend.position='right',
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(color="black",size=15),
          axis.text=element_text(color="black",size=12),
          axis.line.y.right = element_blank(),
          axis.text.y.right=element_blank(),
          axis.ticks.y.right=element_blank(),
          axis.ticks=element_line(color="black",linewidth=0.5),
          plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_blank()
        )


##a. ccRCC1
ID = "ccRCC1"
rdsFile <- "./Rdata/Kidney_scATAC_2Sample.rds"
obj <- readRDS(rdsFile)
color_list <- list(sampleName=c("#FFCC00","#009999"),
                   CellType = color_celltype)
group_by <- c("sampleName","CellType") 
pd <- DimPlot.multi(obj,group_by = group_by,color_list=color_list,legend_position = "right")
ggsave(paste0(plotdir,"/DimPlot_",group_by[1],".pdf"),pd,width =10,height = 5)


CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)

##c. ccRCC2
ID = "ccRCC2"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)
###d.
clonalest <- outres$clonalest
deseg <- DEseg(outres)
merged_df <- deseg$merged_df
res <-deseg$res 
color_r <- cols_Palette[1:length(clonalest)]
names(color_r) <- sort(unique(res$clone))
left_anno_cols <- list()
left_anno_cols[["clone"]] <- color_r

for(j in 1:nrow(merged_df)){
  chr_plt <- merged_df$Chromosome[j]
  segstart <- merged_df$start[j]
  segend <- merged_df$end[j]
  seg_plt <- paste(chr_plt,segstart,segend,sep="_")
  seg_plt

	res_sub <- res[res$Chromosome %in% chr_plt & res$Start >= segstart & res$End<=segend,]
	adjust.dat =2
	adjust.plt = 2
	k=10
	result <- smooth_median_by_group(res_sub,k=k, adjust=adjust.dat)

  pltdata <- result
  integerCNVcolumn <- grep("integerCN",colnames(pltdata))
  labs <- unique(pltdata[integerCNVcolumn])
  labs <- labs[!is.na(labs[,1]),]

  ##Custermrized
  ht1= ggplot() +  
    geom_density(data = subset(pltdata), aes(x = smoothed_binRatio, color = clone), 
               lwd = 1.5, adjust =adjust.plt)
	ht2 <- ht1+
	  geom_vline(xintercept = unique(pltdata$relativeCN_mean),col="grey",linetype = "dashed")+
	  scale_x_continuous(
	    limits = c(0,3),
	      sec.axis = sec_axis( trans = ~.,breaks=unique(pltdata$relativeCN_mean)[!is.na(unique(pltdata$relativeCN_mean))], labels =labs, name="")
	    )+
	  scale_color_manual(values = left_anno_cols[["clone"]])+
	  labs(title=seg_plt,x="Ratio(Bin)")+custom_theme

  ggsave(paste0(plotdir, "/",ID,".density_seg_ColnalBinRatio_",seg_plt,".pdf"),ht2, width=6, height=5,device = pdf,bg="white")  
}



##f.ccRCC3
ID = "ccRCC3"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)
###h
clonalest <- outres$clonalest
deseg <- DEseg(outres)
merged_df <- deseg$merged_df
res <-deseg$res 
color_r <- cols_Palette[1:length(clonalest)]
names(color_r) <- sort(unique(res$clone))
left_anno_cols <- list()
left_anno_cols[["clone"]] <- color_r
for(j in 1:nrow(merged_df)){
  chr_plt <- merged_df$Chromosome[j]
  segstart <- merged_df$start[j]
  segend <- merged_df$end[j]
  seg_plt <- paste(chr_plt,segstart,segend,sep="_")
  seg_plt

	res_sub <- res[res$Chromosome %in% chr_plt & res$Start >= segstart & res$End<=segend,]
	adjust.dat =1
  adjust.plt = 1.5
  k=1
  limits = c(0,3)
	if(seg_plt %in% c("chr16_46688301_90082747","chr8_231970_43141576")){limits=c(0,2)}
	if(seg_plt=="chr16_9930_31874205"){adjust.plt = 2}  
	result <- smooth_median_by_group(res_sub,k=k, adjust=adjust.dat)

  pltdata <- result
  integerCNVcolumn <- grep("integerCN",colnames(pltdata))
  labs <- unique(pltdata[integerCNVcolumn])
  labs <- labs[!is.na(labs[,1]),]

  ##Custermrized
  ht1= ggplot() +  
    geom_density(data = subset(pltdata), aes(x = smoothed_binRatio, color = clone), 
               lwd = 1.5, adjust =adjust.plt)


  if(seg_plt=="chr8_47260174_145052768"){
    ht1= ggplot() + 
      geom_density(data = subset(pltdata, clone == "3"), aes(x = smoothed_binRatio, color = clone), 
                   lwd = 1.5, adjust =2.5) +
      geom_density(data = subset(pltdata, clone %in%c("1","2") ), aes(x = smoothed_binRatio, color = clone), 
                   lwd = 1.5, adjust =2) 
  }

 if(seg_plt=="chr8_231970_43141576"){
    ht1= ggplot() + 
      geom_density(data = subset(pltdata, clone == "1"), aes(x = smoothed_binRatio, color = clone), 
                   lwd = 1.5, adjust =2) +
      geom_density(data = subset(pltdata, clone %in%c("2","3") ), aes(x = smoothed_binRatio, color = clone), 
                   lwd = 1.5, adjust =1.3)
  }

	ht2 <- ht1+
	geom_vline(xintercept = unique(pltdata$relativeCN_mean),col="grey",linetype = "dashed")+
	scale_x_continuous(
	  limits = limits,
	    sec.axis = sec_axis( trans = ~.,breaks=unique(pltdata$relativeCN_mean)[!is.na(unique(pltdata$relativeCN_mean))], labels =labs, name="")
	  )+
	scale_color_manual(values = left_anno_cols[["clone"]])+
	labs(title=seg_plt,x="Ratio(Bin)")+custom_theme


	ggsave(paste0(plotdir, "/",ID,".density_seg_ColnalBinRatio_",seg_plt,".pdf"),ht2, width=6, height=5,device = pdf,bg="white")  
}


##g.ccRCC4
ID = "ccRCC4"
CNVresFile_sc <- paste0("./AnalysisData/final.CNVres_",ID,".rds")
outres <- readRDS(CNVresFile_sc)
cnv_plots_comb(ID,outres,plotdir=plotdir)
###h
clonalest <- outres$clonalest
deseg <- DEseg(outres)
merged_df <- deseg$merged_df[1:2,]
res <-deseg$res 
color_r <- cols_Palette[1:length(clonalest)]
names(color_r) <- sort(unique(res$clone))
left_anno_cols <- list()
left_anno_cols[["clone"]] <- color_r
for(j in 1:nrow(merged_df)){
  chr_plt <- merged_df$Chromosome[j]
  segstart <- merged_df$start[j]
  segend <- merged_df$end[j]
  seg_plt <- paste(chr_plt,segstart,segend,sep="_")
  seg_plt

	res_sub <- res[res$Chromosome %in% chr_plt & res$Start >= segstart & res$End<=segend,]
	adjust.dat =1
  adjust.plt = 1.5
  k=1
  limits = c(0,3)
	result <- smooth_median_by_group(res_sub,k=k, adjust=adjust.dat)

  pltdata <- result
  integerCNVcolumn <- grep("integerCN",colnames(pltdata))
  labs <- unique(pltdata[integerCNVcolumn])
  labs <- labs[!is.na(labs[,1]),]

  ##Custermrized
  ht1= ggplot() +  
    geom_density(data = subset(pltdata), aes(x = smoothed_binRatio, color = clone), 
               lwd = 1.5, adjust =adjust.plt)

	ht2 <- ht1+
	geom_vline(xintercept = unique(pltdata$relativeCN_mean),col="grey",linetype = "dashed")+
	scale_x_continuous(
	  limits = limits,
	    sec.axis = sec_axis( trans = ~.,breaks=unique(pltdata$relativeCN_mean)[!is.na(unique(pltdata$relativeCN_mean))], labels =labs, name="")
	  )+
	scale_color_manual(values = left_anno_cols[["clone"]])+
	labs(title=seg_plt,x="Ratio(Bin)")+custom_theme

	ggsave(paste0(plotdir, "/",ID,".density_seg_ColnalBinRatio_",seg_plt,".pdf"),ht2, width=6, height=5,device = pdf,bg="white")  
}



