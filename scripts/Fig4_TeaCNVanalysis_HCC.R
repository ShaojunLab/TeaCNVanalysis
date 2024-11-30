## Infering absolute CNV from scATAC using TeaCNV

suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(ggsci)
  library(LaplacesDemon)
  library(Seurat)
  library(Signac)
  library(futile.logger)
  library(Matrix) 
  library(TeaCNV)
})
script_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/code/code_CNV/github/TeaCNV/"
setwd(script_path)
genome = "hg38"
blacklistFile <- file.path(".", "data", "hg38-blacklist.v2.bed")

work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'
setwd(work_dir)
data_path = paste0(work_dir, '/Rdata/')

outdir = paste0(work_dir,'/TeaCNV_HCC_MultiOmic')
if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)

inputdata <- paste0(data_path,"/HCC_MultiOmic_12sam_clean.rds")
obj <- readRDS(inputdata)

cell_anno_raw<- data.frame(row.names=rownames(obj@meta.data),
                           Cluster=obj$group,
                           sampleID = obj$dataset,
                           tissue = obj$tissue,
                           celltype = obj$celltype_peaks,
                           log10_nFrags_peaks=obj$log10_nFrags_peaks)

###Filtering low Cholangiocytes and Hepatocytes cells
obj_sub <- subset(obj, celltype_peaks %in% c("Cholangiocytes","Hepatocytes"))
filtcell_res <- FiltCell.obj(obj_sub,assay = "peaks",filt_by=c("nCount_peaks","nCount_RNA"),
                             filt_per_group ="sampleID",by_zscore.x = F,log_ = TRUE)
pcom <- cowplot::plot_grid(plotlist=filtcell_res$plots,ncol = 4,align="v")
Np <- length(filtcell_res$plots)
p_width <- ifelse(Np>=4,24,6*Np)
ggsave(paste0(outdir, "/TumorCells_filtered.pdf"),pcom, width=p_width, height=4*ceiling(Np/4),device = pdf,bg="white")

Hepato_cells <- filtcell_res$cell_keep
nonHepato_cells <- rownames(cell_anno_raw)[!cell_anno_raw$celltype %in%c("Cholangiocytes","Hepatocytes")]
cells_select <- c(nonHepato_cells,Hepato_cells)
cell_anno <- cell_anno_raw[cells_select,,drop=F]
write.csv(cells_select,paste0(outdir, "/cells_filtered.csv"),row.names = F)
#cells_select <- read.csv(paste0(outdir, "/cells_filtered.csv"))$x

rawmtx <- obj@assays$peaks@counts
sampleIDs <- unique(cell_anno$sampleID)
sampleIDs <- sampleIDs[!sampleIDs %in% c("JS230705J1","JS230725J2")]  #Tumor sample
table(cell_anno$sampleID,cell_anno$celltype,useNA="ifany")
#EffectCorrectByRef
rawmtx_correct_ls <- EffectCorrectByRef(rawmtx,cell_anno = cell_anno_raw,
                                        sampleID_normal = c("JS230705J1","JS230725J2"),
                                        NormalTypeList=c("Endothelial","Lymphoid","Myeloid","Fibroblast"),
                                        sampleIDs_to_correct = sampleIDs)
rawmtx_c <- rawmtx_correct_ls$matrix
rawmtx_c <- as(rawmtx_c, "sparseMatrix") 


#Estimate CNV
cell_anno <- cell_anno_raw[cells_select,,drop=F]
ref_group_names <- "Adjacent.normal"
binSize=1
SegSize_min= ifelse(binSize<6,round(100/binSize),20)
seu_resolution = 1.5
sampleIDs <- unique(cell_anno$sampleID)

for(i in 1:length(sampleIDs)){
    clt <- as.character(sampleIDs[i])
    print(clt)
    print(paste0("Analysing cluster ",clt))
    outdir_clt <- paste0(outdir,'/',clt,"/filtered_ana_binSize",binSize);ifelse(!dir.exists(file.path(outdir_clt)), dir.create(file.path(outdir_clt)), FALSE)
    
    cell_anno_sub <- cell_anno[cell_anno$Cluster %in% c(as.character(clt),ref_group_names),,drop=F ]
    cell_anno_sub <- droplevels(cell_anno_sub)

    #Filtering tumor cells
    cells_ref <-  rownames(cell_anno_sub)[cell_anno_sub$Cluster%in%ref_group_names]
    cell_anno_tumor <- cell_anno_sub[!cell_anno_sub$Cluster %in%ref_group_names,,drop=F]
    cutoff<- quantile(cell_anno_tumor$log10_nFrags_peaks,0.2)
    cells_tumor_filted <- rownames(cell_anno_tumor)[as.numeric(cell_anno_tumor$log10_nFrags_peaks)>= cutoff]
    cell_anno_tumor_filt <- cell_anno_tumor[rownames(cell_anno_tumor)%in%cells_tumor_filted,,drop=F]
    cell_anno_new_filt <- cell_anno_sub[rownames(cell_anno_sub)%in%c(cells_tumor_filted,cells_ref),,drop=F]
    cell_anno_sub <- cell_anno_new_filt[order(cell_anno_new_filt$Cluster,cell_anno_new_filt$sampleID),,drop=F]
    
    mtx_sub <- rawmtx_c[,rownames(cell_anno_sub)]

	cnv_obj <- CreateEiCNVObject(input = mtx_sub,
	                             annotationFile = cell_anno_sub,
	                             ref_group_names = ref_group_names,
	                             ChrRemove = c('chrX', 'chrY', 'chrM'),
	                             FiltCell = FALSE,
	                             cellproplim = 0.05,
	                             count_lim = 4)
	if(clt %in% c("JS230907J2")){
      SegSize_min = 150
      seu_resolution <- 1
    }

    if(clt %in% c("JS230915J1")){
      SegSize_min = 60
    }
    if(clt %in% c("JS230907J3")){
      SegSize_min = 100
      seu_resolution <- 1
    }
    res=tryCatch(runEiCNVs(
        input_obj = cnv_obj,
        binSize = binSize,
        outdir = outdir_clt,
        seu_resolution = seu_resolution,
        SegLen_min = 2e6,
        SegSize_min = SegSize_min,
        seg_method = "PELT",
        segValue.method="median",
        StopStep = 4,
        delt_lim = 0.3,
      ),error=function(e) NA)

}


###re-analysis
cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
genome = "hg38"
min_cells_in_group = 20
FiltSeg= TRUE;
SegLen_min = 2e6;
seg_method = "PELT";
seg.count.lim = 80
penalty.lim = c(0.5,1.5)
ref_group_names = "Adjacent.normal"
binSize=1

sampleIDs <- list.files(outdir)
sampleIDs
for(i in 1:length(sampleIDs)){
  clt <- as.character(sampleIDs[i])
  print(clt)
  outdir_clt <- paste0(outdir,'/',clt,"/filtered_ana_binSize",binSize);ifelse(!dir.exists(file.path(outdir_clt)), dir.create(file.path(outdir_clt)), FALSE)
  input_obj <- readRDS(paste0(outdir_clt,"/EiCNVs.obj"))
  cell_anno_new <- input_obj@cell_anno.filt

  mtx_bin <- input_obj@data.binCount.norm
  SegSize_min <- input_obj@options$SegSize_min

  outdir_sub <- paste0(outdir_clt,"/subCluster")

  input_obj@options$delt_lim
  input_obj@options$minCN.frac

  delt_lim=0.3
  minCN.frac=0.01
  refine_segCN_by_dist=F
  scFactor <- 1
  seg_method = "PELT"
  if(clt %in% c("JS230725J1")){
    delt_lim=0.2

  }
  if(clt %in% c("JS230907J1")){
    delt_lim=0.28
  }
  if(clt %in% c("JS230915J2")){
    delt_lim=0.4
    refine_segCN_by_dist=TRUE
  }
   if(clt %in% c("JS230907J2")){
    delt_lim=0.4
    SegSize_min=150 
  }
  if(clt %in% c("JS230830J2")){
     delt_lim=0.3
     refine_segCN_by_dist=TRUE
  }
  if(clt %in% c("JS230712J1")){
     delt_lim=0.2
     minCN.frac=0.01
  }
  if(clt %in% c("JS230915J1")){
    SegSize_min=60 
  }

  if(clt %in% c("JS230901J1")){
    minCN.frac=0.06
  }
if(clt %in% c("JS230919J1")){
   refine_segCN_by_dist=TRUE
  }

  if(clt %in% c("JS230907J3")){
   minCN.frac <- 0.03
   refine_segCN_by_dist=TRUE
   scFactor <- 1.5
  }

best_clone <- input_obj@options$bestClone 
best_clone

CNVresFile_subC <- paste0(outdir_sub,"/step03.CNV_res_clonal.rds")
clonal_res <- readRDS(CNVresFile_subC)

clone.ref <- best_clone
CNest.ref <- clonal_res[[clone.ref]]$CNest
CNest.ref <- peakIndex(clonal_res[[clone.ref]]$input_BinSegRatio,CNest.ref,index_col="SegMean")
ploidy.ref <- clonal_res[[clone.ref]]$ploidy
seg_dat_ref <- clonal_res[[clone.ref]]$seg.dat

cells_obs <- rownames(cell_anno_new)[cell_anno_new$group %in% "observed"]
cells_ref <-  rownames(cell_anno_new)[cell_anno_new$group%in%"reference"]
cluster1 <- cell_anno_new[cells_obs,"subCluster"]
cluster2 <- cell_anno_new[cells_obs,"clone_merged"]

CNVresFile_subC_update <- paste0(outdir_sub,"/step03.CNV_res_clonal_update.rds")
clonal_res_updated <- readRDS(CNVresFile_subC_update)
clonal_res <- clonal_res_updated



  outres <- celloutput(cells.obs=cells_obs,
    cluster1,cluster2,clone_res=clonal_res,mtx_bin=mtx_bin,
                       CNest.ref= CNest.ref,
                       min_cells_in_group,
                       penalty.lim,
                       seg_method=seg_method,
                       FiltSeg,
                       SegLen_min,
                       SegSize_min*scFactor, #(JS230907J3=SegSize_min*1.5)
                       seg.count.lim,
                       cytoBandFile,
                       outdir=outdir_clt,
                       minCN.frac=minCN.frac,
                       seg_dat_ref=seg_dat_ref,
                       correct_by_dist=refine_segCN_by_dist,
                       delt_lim = delt_lim,
                       segValue_method="median",
                       p.adj_cutoff = 0.05,
                       ratio_diff_cutoff=NULL
                       )

 cellbinCount <- input_obj@data.binCount[,outres$cellinfo$cellname]
  outres$cellbinCount <- cellbinCount
  outres$cellbinRatio_raw <- input_obj@data.binRatio
  outres$cellinfo$clone <- factor(outres$cellinfo$clone,levels=sort(unique(as.numeric(outres$cellinfo$clone))))
  outres$CNest.ref <- CNest.ref

  sm_win <- ifelse(binSize<6,ceiling(60/binSize),10)
  plt_mt2 <- smooth_and_denoise(input_obj@data.binRatio,window=sm_win)
  outres$cellbinRatio_deNoise <- plt_mt2

  outres$cellinfo[1:2,]

  clonal_res_final <- outres$clonalest
  clonal_res_final <- lapply(clonal_res_final,function(x){
    clonal_res <- x$input_BinSegRatio
    clonal_res <- clonal_res[,!grepl("relativeCN|integerCN",colnames(clonal_res))]
    x$input_BinSegRatio<- clonal_res
    return(x)
  })

  if(refine_segCN_by_dist){ 
    clonal_res_final <- refine_segCN.dist_v2(clonal_res_final,Nadjacent=1)
  }
  

  if(nrow(CNest.ref)>1){
    delt.ref <- (CNest.ref$ratio[2]-CNest.ref$ratio[1])/(CNest.ref$CN[2]-CNest.ref$CN[1])
    ratio_diff_cutoff <- round(0.8*delt.ref,2)
  }else{
    ratio_diff_cutoff <-NULL
  }

  

  plot_combine_seg(clonal_res_final,ylim=NULL,
    outplot_name="clonalCNV_final_refine0.pdf",show_dots=TRUE,
    outdir=outdir_clt)

 clonal_res_final <-  lapply(clonal_res_final,function(x){
      seg_dat <- x$seg.dat
      seg_dat <- seg_dat[!is.na(seg_dat$integerCN),,drop=F]
      CNest <- x$CNest
      baseCN <- 2
      baseCN_frac <- sum(seg_dat$w[seg_dat$integerCN==baseCN])/sum(seg_dat$w)
      aneuploidy_score <- 1- baseCN_frac
      x$score$aneuploidy_score <- aneuploidy_score
      return(x)
    })

  outres$clonalest <- clonal_res_final
  aneuploidy_score <- unlist(lapply(outres$clonalest,function(x){
  x$score$aneuploidy_score
  }))
  clone_anno <- data.frame(clone_merged=names(outres$clonalest),aneuploidy_score=aneuploidy_score)
  cell_anno_new<- cell_anno_new[,!grepl("aneuploidy_score",colnames(cell_anno_new))]
  cell_anno_new <- left_join(cell_anno_new,clone_anno,by=intersect(colnames(cell_anno_new),colnames(clone_anno)))
  rownames(cell_anno_new) <- cell_anno_new$row
  
  cell_anno_new_file <- paste0(outdir_clt,"/subCluster/cell_info_subClust.csv")
  write.csv(cell_anno_new,cell_anno_new_file,row.names = T)

  CNVresFile_sc <- paste0(outdir_clt,"/final.CNVres-1015.rds")
  #outres <- readRDS(CNVresFile_sc)  
  saveRDS(outres,CNVresFile_sc)
  saveRDS(input_obj,paste0(outdir_clt,"/EiCNVs.obj"))  
 
  ###clonal paired-comparison

    cell_anno_cl <- outres$cellinfo[,c("clone"),drop=FALSE]
    cell_anno_cl <- cell_anno_cl[!is.na(cell_anno_cl$clone),,drop=F]
    clones <- as.character(sort(as.numeric(as.character(unique(cell_anno_cl$clone)))))
    Nclone <- length(clones)
    outdir_clt4 <- paste0(outdir_clt,"/segRatio_compare_clone")
    ifelse(!dir.exists(file.path(outdir_clt4)), dir.create(file.path(outdir_clt4)), FALSE)
    mtx_bin <- input_obj@data.binCount.norm
    
    if(Nclone>1){
      seg_ls <- list()
      for(col_ref in clones){
        SegScore_total <- seg4clusters(mtx_bin[,rownames(cell_anno_cl)],cell_anno_cl,
                                       cytoBand =cytoBandFile,
                                       col_ref,
                                       doFilt=T,
                                       seg.penalty = c(0.5,1),
                                       seg.method = seg_method,#"gaussian",
                                       segSize.min = SegSize_min,
                                       outplot_path=outdir_clt4,
                                       outname=paste0("segRatio_cloneX_vsC",col_ref,".pdf"))
        
        seg_ls[[col_ref]] <- SegScore_total$seg_score_binLevel
        
      }
      saveRDS(seg_ls,paste0(outdir_sub,"/clonal_RelativeRatio_segment.rds"))
    }



{

   #plot segment (final)
    plot_combine_seg(outres$clonalest,ylim=NULL,
      outplot_name="clonalCNV_final_refine0.pdf",show_dots=TRUE,
      outdir=outdir_clt)
    plot_combine_seg(outres$clonalest,ylim=NULL,
      outplot_name="clonalCNV_final_refine1.pdf",show_dots=FALSE,
      outdir=outdir_clt)

   p_ls <- lapply(sort(names(outres$clonalest)),function(cluster,outdir){
      clonal_res <- outres$clonalest[[cluster]]$input_BinSegRatio
      integerCNV <- outres$clonalest[[cluster]]$seg.dat
      score <-  round(outres$clonalest[[cluster]]$score$score,2)
      ymax <- max(integerCNV$integerCN)
      plt <- cellLineSeg_v3(binRatio.df=clonal_res,integerCNV.df=integerCNV,outdir=outdir_clt,
                            ylab="CNAs",ylim=c(0,ymax),
                            plot_dots = FALSE,
                            plot_hist = FALSE,
                            color_seg_gradient=T, seg.color = "#0d3b66",
                            color_hist_gradient=T,hist_color= NULL,#"#61210f",
                            fname=paste0("C",cluster,": score=",score), ggarrange_width=c(10,5),
                            height = 2,width=10,outPlot=FALSE,
                            color_limit = c(1,6),
                            value.bin="integerCN",
                            value.segment="integerCN",
                            label_CN=FALSE)
      return(plt$p1+ theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5),'lines')))
    },outdir_clt)
    names(p_ls) <- sort(names(outres$clonalest))
    pcom <- cowplot::plot_grid(plotlist=p_ls,ncol = 1,align="v")
    ggsave(paste0(outdir_clt, "/integerCN-1015.pdf"),pcom, width=10, height=2*length(p_ls),device = pdf,bg="white")

   #segment clonal CNAs (grouped)
    pl_grp <- c()
    pl_hist <- c()
    ymin <- min(unlist(lapply(outres$clonalest,function(x){x$seg.dat$integerCN})))-0.5
    ymax <- max(unlist(lapply(outres$clonalest,function(x){x$seg.dat$integerCN})))+0.5
    colorss <- c("#023e8a","grey50", "#cb793a", "#9a031e","#6a040f","#370617")
    for(cluster_indx in 1:length(names(outres$clonalest))){

      cluster <- sort(names(outres$clonalest))[cluster_indx]
      clonal_res <- outres$clonalest[[cluster]]$input_BinSegRatio
      integerCNV <- outres$clonalest[[cluster]]$seg.dat
      if(cluster_indx>1){
        d <- ifelse(cluster_indx%% 2 == 0,0.1*(cluster_indx/2),-0.1*(cluster_indx%/%2))
      }else{d=0}
      integerCNV$integerCN <- integerCNV$integerCN+d

      plt <- cellLineSeg_v3(binRatio.df=clonal_res,integerCNV.df=integerCNV,outdir=outdir_clt,
                            ylab="CNAs",ylim=c(0,ymax),
                            plot_dots = FALSE,
                            plot_hist = FALSE,
                            color_seg_gradient=T, seg.color = "#0d3b66",
                            color_hist_gradient=T,hist_color= NULL,#"#61210f",
                            fname="", ggarrange_width=c(10,5),
                            height = 2,width=10,outPlot=FALSE,
                            color_limit = c(1,6),
                            value.bin="integerCN",
                            value.segment="integerCN",
                            label_CN=FALSE)
      pl_grp[[cluster]] <- plt$p1
    }
    prg <- pl_grp[[1]]
    if(length(pl_grp)>1){
      for(p_i in 2:length(pl_grp)){
        cluster <- names(pl_grp)[p_i]
        a <- pl_grp[[cluster]]$data
       
       prg <- prg+
          geom_segment(
          data = a, 
          mapping = aes(x=segStrat, y=value.segment, xend=segEnd, yend=value.segment,colour=integerCN), 
          na.rm =T,size = 1.5,
          inherit.aes = FALSE)+
        scale_color_gradientn(name = "",colours = colorss,
                              limits  = c(1,6),
                              oob = scales::squish)
      }

    }
    prg <- prg +
          theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), 'lines'), #c(top, right, bottom, left)
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position = "none",plot.title = element_text(size = 15),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 15))
    ggsave(paste0(outdir_clt,"/integerCN_combined_1.pdf"),prg,width=15, height=3.5,device = pdf,bg="white")
   ######
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
        for(pl_i in 2:length(pl_data)){
            cluster <- names(pl_data)[pl_i]
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
    ggsave(paste0(outdir_clt,"/integerCN_combined.pdf"),p1,height = 2,width=10)

    cellMeta <- cell_anno_new[!cell_anno_new$subCluster %in%"reference",,drop=F]
    cellMeta  <- cellMeta[!cellMeta$clone_merged %in%"removed",,drop=F]
    cellMeta$subCluster <- factor(cellMeta$subCluster,levels=sort(unique(as.numeric(cellMeta$subCluster))))
    cellMeta$clone_merged <- factor(cellMeta$clone_merged,levels=sort(unique(as.numeric(cellMeta$clone_merged))))
    cellMeta <- cellMeta[order(cellMeta$clone_merged,cellMeta$subCluster),,drop=F]
    colnames(cellMeta) <- gsub("clone_merged","clone",colnames(cellMeta))
    cellMeta <- cellMeta %>%
      add_count(clone) %>%
      as.data.frame()
    colnames(cellMeta)[colnames(cellMeta)=="n"] <- "Ncells"
    rownames(cellMeta)<- cellMeta$row
    cellMeta$CellProportion <- cellMeta$Ncells/nrow(cellMeta)
    cloneInfo_table <- unique(cellMeta[,c("clone","aneuploidy_score","Ncells","CellProportion")])
    write.csv(cloneInfo_table,paste0(outdir_clt,"/cloneInfo_table.csv"),row.names=F)


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

    # ##scale clonal size
    Ncol <- round(cloneInfo_table$CellProportion*100)
    new_CNmt <- do.call(cbind, lapply(1:ncol(CNmt), function(i) {
      matrix(rep(CNmt[, i], Ncol[i]), ncol = Ncol[i])
    }))
    rownames(new_CNmt) <- rownames(CNmt)
    colnames(new_CNmt) <- seq(1:ncol(new_CNmt))
    clone_info <- data.frame(row.names=colnames(new_CNmt),clone=rep(colnames(CNmt), Ncol))
    color_r <- cols_Palette[1:length(unique(clone_info$clone))]
    names(color_r) <- sort(unique(clone_info$clone))
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r
    p_cloneCN <- heatmap4peakMt(mat=new_CNmt,
                          meta_info=clone_info,
                          sep_by="-",
                          outdir= outdir_clt,
                          value.type="CNV",
                          clust_rows=F,
                          show_legend_row = T,
                          legend_titles="integer CN",
                          fileout_name=paste0("heatmap_cloneCNA_SizeScale"),
                          col_list=left_anno_cols,
                          column.title = NULL,
                          width=10,height=5,device="pdf") #height=5



    clone_info <- data.frame(row.names=colnames(CNmt),clone=colnames(CNmt))
    # color_r <- suppressWarnings(get_group_color_palette("Dark2")(length(unique(clone_info$clone))))
    color_r <- cols_Palette[1:length(unique(clone_info$clone))]
    names(color_r) <- sort(unique(clone_info$clone))
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r
    height <- ifelse(ncol(CNmt)>2,0.35*(ncol(CNmt))+0.5,ifelse(ncol(CNmt)==1,1.2,1.5))
    p_cloneCN <- heatmap4peakMt(mat=CNmt,
                          meta_info=clone_info,
                          sep_by="-",
                          outdir= outdir_clt,
                          value.type="CNV",
                          clust_rows=F,
                          show_legend_row = T,
                          legend_titles="integer CN",
                          fileout_name=paste0("heatmap_cloneCNA_EqualSize"),
                          col_list=left_anno_cols,
                          column.title = NULL,
                          width=10,height=height,device="pdf") #height=5
}

  
  if(length(outres$clonalest)>1){
  	  infoSeg_df <- find_infoSeg(outres$clonalest,seg_len_min=SegLen_min,doPlot=TRUE,outdir=outdir_clt,
	  p.adj_cutoff = 0.05,
	  ratio_diff_cutoff=ratio_diff_cutoff)
	  write.csv(infoSeg_df$merged_df,paste0(outdir_clt,"/infoSeg_",clt,".csv"))
	  write.csv(infoSeg_df$CN_InfoSegs_mt,paste0(outdir_clt,"/CN_InfoSegs_mt_",clt,".csv"))

  }
 
}







