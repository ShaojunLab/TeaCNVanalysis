##Rare clone data
###Simulate scDNA-seq
library(dplyr)
library(TeaCNV)

script_path="./github/TeaCNV/"
cytoBandFile <- file.path(script_path, "data", "cytoBand_hg38.tsv")

wkdir <- "./TeaCNV"
outdir <- paste0(wkdir,"/simulation/simData/scDNA")
if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)

##Load ground truth
cnvpath <- "./TeaCNV"
cnvresFile <- paste0(cnvpath,"/ccRCC4/final.CNVres.rds")
cnvres <- readRDS(cnvresFile)

seg.dat_tr1 <-  cnvres$clonalest[["1"]]$seg.dat
bin_dat_tr1 <- cnvres$clonalest[["1"]]$input_BinSegRatio
seg.dat_tr1 <- unique(seg.dat_tr1[,c("segName","integerCN")])
bin_dat_tr1 <- left_join(bin_dat_tr1,seg.dat_tr1,by="segName")
bin_dat_tr1 <- bin_dat_tr1 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)
dim(bin_dat_tr1)
seg.dat_tr2 <-  cnvres$clonalest[["2"]]$seg.dat
bin_dat_tr2 <- cnvres$clonalest[["2"]]$input_BinSegRatio
seg.dat_tr2 <- unique(seg.dat_tr2[,c("segName","integerCN")])
bin_dat_tr2 <- left_join(bin_dat_tr2,seg.dat_tr2,by="segName")
bin_dat_tr2 <- bin_dat_tr2 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)
dim(bin_dat_tr2)
seg.dat_tr3 <-  cnvres$clonalest[["3"]]$seg.dat
bin_dat_tr3 <- cnvres$clonalest[["3"]]$input_BinSegRatio
seg.dat_tr3 <- unique(seg.dat_tr3[,c("segName","integerCN")])
bin_dat_tr3 <- left_join(bin_dat_tr3,seg.dat_tr3,by="segName")
bin_dat_tr3 <- bin_dat_tr3 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)
dim(bin_dat_tr3)
bins <- intersect(bin_dat_tr2$binID,bin_dat_tr3$binID)
bins <- intersect(bins,bin_dat_tr1$binID)

##dataset 1
RarePct <- 0.1
NcRatio <- c(0.5,0.4,0.1)


n_bins <- length(bins)
size <- 2000
bin_names <- bins
read_depth <- 100
delt_lim <- 0.4

set.seed(123) 

###Simulate three clones (rare 10%)
##clone 1
n_cells <- size*NcRatio[1]
cell_names <- paste0("C1_", 1:n_cells)
bin_dat_tr1 <- bin_dat_tr1[match(bins,bin_dat_tr1$binID),]
true_cnv <- bin_dat_tr1$integerCN

count_matrix <- matrix(nrow = n_bins, ncol = n_cells)
for (i in 1:n_bins) {
  count_matrix[i, ] <- rpois(n_cells, lambda = read_depth * true_cnv[i])
}

noise <- matrix(rnorm(n_bins * n_cells, mean = 1, sd = 10), nrow = n_bins)
count_matrix <- round(count_matrix*noise)
count_matrix[count_matrix < 0] <- 0  # 计数不能为负数

cnv_mt <-matrix(rep(true_cnv, n_cells), ncol = n_cells)

rownames(count_matrix) =rownames(cnv_mt) <- bin_names
colnames(count_matrix) =colnames(cnv_mt) <- cell_names

##clone 2
n_cells2 <- size*NcRatio[2]
cell_names2 <- paste0("C2_", 1:n_cells2)
bin_dat_tr2 <- bin_dat_tr2[match(bins,bin_dat_tr2$binID),]
true_cnv2 <- bin_dat_tr2$integerCN

count_matrix2 <- matrix(nrow = n_bins, ncol = n_cells2)
for (i in 1:n_bins) {
  count_matrix2[i, ] <- rpois(n_cells2, lambda = read_depth * true_cnv2[i])
}
noise <- matrix(rnorm(n_bins * n_cells2, mean = 1, sd = 10), nrow = n_bins)
count_matrix2 <- round(count_matrix2 * noise)
count_matrix2[count_matrix2 < 0] <- 0  
cnv_mt2 <-matrix(rep(true_cnv2, n_cells2), ncol = n_cells2)

rownames(count_matrix2) = rownames(cnv_mt2) <- bin_names
colnames(count_matrix2) = colnames(cnv_mt2) <- cell_names2

#clone 3
n_cells3 <- size*NcRatio[3]
cell_names3 <- paste0("C3_", 1:n_cells3)
bin_dat_tr3 <- bin_dat_tr3[match(bins,bin_dat_tr3$binID),]
true_cnv3 <- bin_dat_tr3$integerCN

count_matrix3 <- matrix(nrow = n_bins, ncol = n_cells3)
for (i in 1:n_bins) {
  count_matrix3[i, ] <- rpois(n_cells3, lambda = read_depth * true_cnv3[i])
}
noise <- matrix(rnorm(n_bins * n_cells3, mean = 1, sd = 10), nrow = n_bins)
count_matrix3 <- round(count_matrix3 * noise)
count_matrix3[count_matrix3 < 0] <- 0  
cnv_mt3 <-matrix(rep(true_cnv3, n_cells3), ncol = n_cells3)

rownames(count_matrix3) = rownames(cnv_mt3)<- bin_names
colnames(count_matrix3) = colnames(cnv_mt3)<- cell_names3

cnv_MT <- cbind(cnv_mt,cnv_mt2,cnv_mt3)
count_mt <- cbind(count_matrix,count_matrix2,count_matrix3)
write.csv(count_mt, file = paste0("simulated_scDNA_count_matrix_Triclonal_Obs_",RarePct,".csv"), quote = FALSE)

###Simulate normal reference
n_bins <- length(bins)
size <- 2000
bin_names <- bins
read_depth <- 100
n_cells <- size
cell_names <- paste0("Ref_", 1:n_cells)
count_matrix_ref <- matrix(nrow = n_bins, ncol = n_cells)
for (i in 1:n_bins) {
  count_matrix_ref[i, ] <- rpois(n_cells, lambda = read_depth * 2)
}
rownames(count_matrix_ref) <- bin_names
colnames(count_matrix_ref) <- cell_names
write.csv(count_matrix_ref, file = paste0("simulated_scDNA_count_matrix_Triclonal_Ref_",RarePct,".csv"), quote = FALSE)

# ref_mean <- rowMeans(count_matrix_ref)
# ratio <- count_mt/ref_mean

##Heatmap

p_scRatio <- heatmap4peakMt(mat=cnv_MT,
                           meta_info=NULL,
                           sep_by="-",
                           outdir= outdir,value.type="CNV",
                           legend_titles="CNV",
                           clust_rows=F,clustering_method_rows = "ward.D2",
                           show_legend_row = T,
                           fileout_name=paste0("heatmap_scDNA_",RarePct),
                           width=10,height=5,
                           device="pdf")


##TeaCNV for 
counts <- cbind(count_matrix_ref,count_mt)
counts <- na.omit(counts)
cell_anno<- data.frame(row.names=colnames(counts),
                               Cluster=c(rep("normal", ncol(count_matrix_ref)), rep("tumor", ncol(count_mt)))
                               )
counts <- as.matrix(counts)
mtx_bin = t(t(counts)*mean(colSums(counts))/colSums(counts))

segRatio_bulk <- CNratioInfer(input_matrix=mtx_bin,
                                  inputcellano=cell_anno,
                                  cytoBand=cytoBandFile,
                                  outdir=outdir,
                                  group2bulkby="Cluster",
                                  min_cells_in_group = 20,
                                  normallabel="normal",
                                  seg_method="PELT",
                                  penalty = c(0.5,1.5),
                                  FiltSeg= TRUE,
                                  SegLen_min = 2e6,
                                  SegSize_min = 110,
                                  plot_ylim = NULL,
                                  outFigure =FALSE,
                                  outFigureName='',
                                  color_dot=FALSE,
                                  output_ref.value=F,
                                  outRDS=FALSE,
                                  RDSname = paste0("segRatio_bulk.rds"))

CNV_res_bulk <- CNV_esti_ini(outdir=outdir,
                         segScore_file=segRatio_bulk,
                         filt_seg = TRUE,
                         length_seg_cutoff =1e6,
                         genome="hg38",
                         true_loc = TRUE,
                         outputFigure = F,
                         outplot_name=paste0("CNV_bulk_initial"))

bulkCN_res <- ploidyRefine(CNV_res_bulk,delt.lim =delt_lim )
saveRDS(bulkCN_res,paste0(outdir, "/pseudobulkCNVres.rds"))
ratiores <- bulkCN_res$tumor$input_BinSegRatio
integerCNVres <- bulkCN_res$tumor$seg.dat
plt <- cellLineSeg_v3(binRatio.df=ratiores,integerCNV.df=integerCNVres,outdir=outdir_clt,
                      #ylim = ylim,
                      fname=paste0("Pseudobulk(scDNA)"),
                      ggarrange_width=c(10,4),
                      height = 2,width=10,outPlot=FALSE)
ggsave(paste0(outdir, "/PseudobulkDNA_CNVs_",RarePct,".pdf"),plt$ggarranged_p, width=10, height=2,device = pdf,bg="white")



