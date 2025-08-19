library(dplyr)
library(ggplot2)
library(reshape2)
library(mclust)
	library(tidyr)
	library(bedtoolsr)
	 library(TeaCNV)

source("./github/TeaCNV/simulation/funs_evaluation_SimData.R")
script_path <- "./github/TeaCNV"
setwd(script_path)
source("./ana/funs_evaluation.r")
refbedfile <- file.path(script_path, "data", "hg38.100Kb.windows.sorted.bed")
bed.ref <- read.table(refbedfile,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

outdir <- "./TeaCNV/simulation/RareCloneEva";if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)

tumor_size <- c(50,seq(100,1000,by=100),2000,5000,10000)

clonal_group <- c("Biclonal","Triclonal","Tetraclonal","Pentaclonal")

##(1)evaluate Rare clone identification
##(2)Rare-clonal CNV events identification
###Load ground truth
cnvpath <- "./TeaCNV/Results/multiOmic/ccRCC/TeaCNV"
cnvresFile <- paste0(cnvpath,"/ccRCC4/final.CNVres.rds")
cnvres <- readRDS(cnvresFile)

seg.dat_tr1 <-  cnvres$clonalest[["1"]]$seg.dat
bin_dat_tr1 <- cnvres$clonalest[["1"]]$input_BinSegRatio
seg.dat_tr1 <- unique(seg.dat_tr1[,c("segName","integerCN")])
bin_dat_tr1 <- left_join(bin_dat_tr1,seg.dat_tr1,by="segName")
bin_dat_tr1 <- bin_dat_tr1 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)%>%
mutate(clone=1)
dim(bin_dat_tr1)
seg.dat_tr2 <-  cnvres$clonalest[["2"]]$seg.dat
bin_dat_tr2 <- cnvres$clonalest[["2"]]$input_BinSegRatio
seg.dat_tr2 <- unique(seg.dat_tr2[,c("segName","integerCN")])
bin_dat_tr2 <- left_join(bin_dat_tr2,seg.dat_tr2,by="segName")
bin_dat_tr2 <- bin_dat_tr2 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)%>%
mutate(clone=2)
dim(bin_dat_tr2)
seg.dat_tr3 <-  cnvres$clonalest[["3"]]$seg.dat
bin_dat_tr3 <- cnvres$clonalest[["3"]]$input_BinSegRatio
seg.dat_tr3 <- unique(seg.dat_tr3[,c("segName","integerCN")])
bin_dat_tr3 <- left_join(bin_dat_tr3,seg.dat_tr3,by="segName")
bin_dat_tr3 <- bin_dat_tr3 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)%>%
mutate(clone=3)
dim(bin_dat_tr3)
seg.dat_tr4 <-  cnvres$clonalest[["4"]]$seg.dat
bin_dat_tr4 <- cnvres$clonalest[["4"]]$input_BinSegRatio
seg.dat_tr4 <- unique(seg.dat_tr4[,c("segName","integerCN")])
bin_dat_tr4 <- left_join(bin_dat_tr4,seg.dat_tr4,by="segName")
bin_dat_tr4 <- bin_dat_tr4 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)%>%
mutate(clone=4)
dim(bin_dat_tr4)

cnvresFile2 <- paste0(cnvpath,"/ccRCC3/final.CNVres.rds")
cnvres2 <- readRDS(cnvresFile2)
seg.dat_tr5 <-  cnvres2$clonalest[["2"]]$seg.dat
bin_dat_tr5 <- cnvres2$clonalest[["2"]]$input_BinSegRatio
seg.dat_tr5 <- unique(seg.dat_tr5[,c("segName","integerCN")])
bin_dat_tr5 <- left_join(bin_dat_tr5,seg.dat_tr5,by="segName")
bin_dat_tr5 <- bin_dat_tr5 %>% dplyr::select(Chromosome,Start,End,integerCN,segName,binID)%>%
mutate(clone=5)
dim(bin_dat_tr5)

bins <- Reduce(intersect, list(bin_dat_tr1$binID,bin_dat_tr2$binID,bin_dat_tr3$binID,bin_dat_tr4$binID,bin_dat_tr5$binID))
trcnv_File <-paste0(outdir,"/GroundTruth_df_diffSize.csv")
if(file.exists(trcnv_File)){ GT <- read.csv(trcnv_File)}else{
	GT <- rbind(bin_dat_tr1[match(bins,bin_dat_tr1$binID),],
	bin_dat_tr2[match(bins,bin_dat_tr2$binID),],
	bin_dat_tr3[match(bins,bin_dat_tr4$binID),],
	bin_dat_tr4[match(bins,bin_dat_tr4$binID),],
	bin_dat_tr5[match(bins,bin_dat_tr5$binID),])
	write.csv(GT,trcnv_File,row.names=FALSE)
}



trFile <-paste0(outdir,"/list_ground_truth_100kb_diffSize.rds")
gt_ls<- list(bin_dat_tr1,bin_dat_tr2,bin_dat_tr3,bin_dat_tr4,bin_dat_tr5)

GT_100k <- lapply(1:5,function(x,gt_ls){
	data_sub <- gt_ls[[x]] %>%
			tibble::column_to_rownames("binID")  %>% 
			dplyr::select("chr"=Chromosome,"start"=Start,"end"=End, segName ,integerCN)
	align <- align_Grange2bin(bed.query=bed.ref,bed.subject=data_sub)
	return(align)
},gt_ls)
saveRDS(GT_100k,trFile)


cellinfo_all <- c()

res_perform_clonal <- c()
res_perform_rare <- c()

for(i in 1:length(clonal_group)){
	cg <- clonal_group[i]
	cnvDatadir <- paste0("./TeaCNV/simulation/TeaCNV_res/",cg)

	for(j in 1:length(tumor_size)){
		size <- tumor_size[j]
		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",cg,"-Size",size,"-",nth,"..."))
			cat("\n")

			cnvFile <- file.path(cnvDatadir,paste0("size",size),nth,"final.CNVres.rds")
			cnvres <- readRDS(cnvFile)
			cellinfo <- cnvres$cellinfo
			cellinfo$GroundTruth <- sapply(strsplit(cellinfo$cellname,"_"),"[",1)
			cellinfo$GroundTruth <- gsub("C","",cellinfo$GroundTruth)

			cellinfo$nth <- nth
			cellinfo$size <- size
			cellinfo$group <- cg

			cellinfo_all <- rbind(cellinfo_all,cellinfo)

			cellinfo <- droplevels(cellinfo)
			evares_ls <- evaluateRareCloneAccuracy(cellinfo,truth_col="GroundTruth",pred_col="clone",rare_prop=0.12,dominance_threshold = 0.8)
			evares <- data.frame(group=cg,size=size,nth=nth,
				accuracy=as.numeric(evares_ls$accuracy),
				precision=as.numeric(evares_ls$precision),
				recall=as.numeric(evares_ls$recall),
				F1=as.numeric(evares_ls$f1),
				BACC=as.numeric(evares_ls$balanced_accuracy))
			res_perform_rare <- rbind(res_perform_rare,evares)

			###评估亚克隆结构鉴定（BACC）
			cellinfo2 <- cellinfo[!is.na(cellinfo$clone),,drop=FALSE]
			truth <- as.character( cellinfo2[["GroundTruth"]])
      pred <- as.character( cellinfo2[["clone"]])
			res_clone <- evaluate_MultiClassification(truth,pred)
 			evares_clone <- data.frame(group=cg,size=size,nth=nth,
				accuracy=as.numeric(res_clone$accuracy),
				precision=as.numeric(res_clone$precision),
				recall=as.numeric(res_clone$recall),
				F1=as.numeric(res_clone$f1),
				BACC=as.numeric(res_clone$balanced_accuracy)
				)
 			res_perform_clonal <- rbind(res_perform_clonal,evares_clone)

		}
	}
}
write.csv(res_perform_rare,paste0(outdir,"/evaRes_TeaCNV_diffSize_RareCloneIdentify.csv"),row.names=FALSE)
write.csv(res_perform_clonal,paste0(outdir,"/evaRes_TeaCNV_diffSize_CloneStructureIdentify.csv"),row.names=FALSE)


saveRDS(cellinfo_all,paste0(outdir,"/cellinfo_all_RareClone.rds"))





###plotRareCloneEvaluation()
library(ggplot2)
library(reshape2)
library(patchwork)
metric_cols = c("accuracy", "precision", "recall", "F1","ARI")
type_col = "type"
size_col = "size"

eval_df <- read.csv(paste0(outdir,"/evaRes_TeaCNV_diffSize_RareCloneIdentify.csv"))
eval_df[is.na(eval_df)] <- 0
long_df <- eval_df %>%
	dplyr::select(all_of(c(type_col, size_col,"nth", metric_cols))) %>%
	reshape2::melt(id.vars = c(size_col,type_col,"nth"), variable.name = "Metric", value.name = "Value")

summary_df <- long_df %>%
	group_by( .data[[size_col]], Metric) %>%
	summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")

color_box <- "#cccccc"   
color_line <- "#0072B2"
plots <- list()
for (metric_name in metric_cols) {
	plot_df <- long_df %>% dplyr::filter(Metric == metric_name)
	summary_plot_df <- summary_df %>% dplyr::filter(Metric == metric_name)

	y_title <- paste0(metric_name," of rare clone assignment")
	x_title <- "The number of cells in rare clone"
	p <- ggplot(plot_df, aes(x = factor(.data[[size_col]]/10), y = Value)) +
	  geom_boxplot(fill = color_box, alpha = 0.4, outlier.shape = NA, width = 0.6) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[[size_col]]/10), y = mean_value),color = color_line,
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[[size_col]]/10), y = mean_value),color = color_line,
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = "", x = x_title, y = y_title) +
      theme_classic(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )
   plots[[metric_name]] <- p
}

final_plot <- wrap_plots(plots, ncol = 2)
ggsave(paste0(outdir,"/accuracy_precision_recall_f1_ARI.pdf"),final_plot,width=12,height=6)


###Statistics of running time and memory
info_total <- c()
for(i in 1:length(clonal_group)){
	cg <- clonal_group[i]
	cnvDatadir <- paste0("./TeaCNV/simulation/TeaCNV_res/",cg)

	for(j in 1:length(tumor_size)){
		size <- tumor_size[j]
		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",cg,"-size",size,"-",nth,"..."))
			cat("\n")
			txtFile <- file.path(cnvDatadir,paste0("size",size),nth,"runtime_log.txt")
			if(file.exists(txtFile)){
				log_t <- extract_runtime_info(txtFile)
			}
			log_t$nth <- nth
	  	log_t$group <- cg
	  	log_t$size <- size
	  	info_total <- rbind(info_total,log_t)
		}
	}
}
info_total$Memory_MB[info_total$Memory_MB<0] <- NA
write.csv(info_total,paste0(outdir,"/runtime_log_total_diffSize.csv"),row.names=FALSE)
##plot
metric_cols <- c('Memory_MB','Runtime_Minutes')
long_df <- info_total %>%
  dplyr::select(all_of(c("size","group", metric_cols))) %>%
  reshape2::melt(id.vars = c( "size","group"), variable.name = "Metric", value.name = "Value")%>%
  as.data.frame()

summary_df <- long_df %>%
  group_by(size, Metric) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
  sd_value = sd(Value, na.rm = TRUE), .groups = "drop")
write.csv(summary_df,paste0(outdir,"/eva_diffSize_runtime.summary.csv"),row.names=FALSE)
long_df %>%
  group_by(Metric) %>%
  summarise(median_value = median(Value, na.rm = TRUE), .groups = "drop")



color_box <- "#cccccc"   
color_line <- "#0072B2"
plots <- list()
for (metric_name in metric_cols) {
	plot_df <- long_df %>% dplyr::filter(Metric == metric_name)
	summary_plot_df <- summary_df %>% dplyr::filter(Metric == metric_name)
	ylim <- c(0,summary_plot_df$mean_value[summary_plot_df$size==10000]*1.2)
	p <- ggplot(plot_df, aes(x = factor(.data[["size"]]), y = Value)) +
	  geom_boxplot(fill = color_box, alpha = 0.4, outlier.shape = NA, width = 0.6) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[["size"]]), y = mean_value),color = color_line,
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[["size"]]), y = mean_value),color = color_line,
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = "", x = "Sample Size", y = metric_name) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )+
      coord_cartesian(ylim = ylim)
	  plots[[metric_name]] <- p
}

final_plot <- wrap_plots(plots, ncol = 1)
ggsave(paste0(outdir,"/eva_runtime_diffSize.pdf"),final_plot,width=7,height=5)



###-----------------------------------###
###(2)###evaluate Rare clone CNVs identification
###-----------------------------------###
outdir <- "./TeaCNV/simulation/RareCloneEva"

tumor_size <- c(50,seq(100,1000,by=100),2000,5000,10000)
clonal_group <- c("Biclonal","Triclonal","Tetraclonal","Pentaclonal")


res_clonal_file_tea <- file.path(outdir,"res_clonalCNV_freq_TeaCNV_diffSize.csv")
# if (file.exists(res_clonal_file_tea)) file.remove(res_clonal_file_tea)
res_rare_file <- file.path(outdir,"res_rareCNVevents_identify_performance_TeaCNV_diffSize.csv")
# if (file.exists(res_rare_file)) file.remove(res_rare_file)



###Ground truth CNV
trcnv_File <-paste0(outdir,"/GroundTruth_df_diffSize.csv")
trFile <-paste0(outdir,"/list_ground_truth_100kb_diffSize.rds")
if(file.exists(trFile)){ GT_100k <- readRDS(trFile)}else{stop(paste0(trFile," NOT exists!"))}
if(file.exists(trcnv_File)){ GT <- read.csv(trcnv_File)}

truth_wide <- GT %>%
	dplyr::select(binID,integerCN,clone)%>%
	dplyr::distinct()%>%
  pivot_wider(
    names_from = clone,
    values_from = integerCN
  )%>%
	tibble::column_to_rownames(var = "binID")


for(i in 1:length(clonal_group)){
	group <- clonal_group[i]
	cnvDatadir <- paste0("./TeaCNV/simulation/TeaCNV_res/",group)

	truth_wide_j <- truth_wide[,1:(i+1)]
	CNVinRow <- apply(truth_wide_j,1,unique)
	cnvRow<- unlist(lapply(CNVinRow,function(x){length(x[!is.na(x)])>1}))
	truth_wide_clonalCNV <- truth_wide_j[cnvRow,,drop=F]
	CNVs_clonal <- rownames(truth_wide_clonalCNV)

	CNVs_clonal_bed <- strings2bed(CNVs_clonal)
	CNVs_clonal_bed <- cbind(CNVs_clonal_bed,truth_wide_clonalCNV)

	rareCNV_index <- apply(truth_wide_clonalCNV,1,function(x,i){
		(x[i+1]!=2) && (!x[i+1] %in% x[1:i])
	},i)
	if(i==1){
		rareCNV_index <- apply(truth_wide_clonalCNV,1,function(x,i){
			(!x[i+1] %in% x[1:i])
		},i)
	}

	CNVs_rare <- rownames(truth_wide_clonalCNV)[rareCNV_index]
	CNVs_rare_bed <- strings2bed(CNVs_rare)
	CNVs_rare_bed100 <- align_Grange2bin(bed.query=bed.ref,bed.subject=CNVs_rare_bed)
	rareCNVRegion_index <- !is.na(CNVs_rare_bed100$chr.subject)




	for(j in 1:length(tumor_size)){
		size <- tumor_size[j]
		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",group,"-Size",size,"-",nth,"..."))
			cat("\n")

			cnvFile <- file.path(cnvDatadir,paste0("size",size),nth,"final.CNVres.rds")
			outres <- readRDS(cnvFile)
			clonalCN_res <- outres$clonalest
			clonalCN_res[["removed"]]<- NULL

			cellinfo <- outres$cellinfo
			cellinfo$truthClone <- sapply(strsplit(cellinfo$cellname,"_"),"[",1)
			cellinfo$truthClone <- gsub("C","",cellinfo$truthClone)
			N_truth <- table(cellinfo$truthClone)
			prop_truth <- N_truth/nrow(cellinfo)
			N_pred <- table(cellinfo$clone)
			prop_pred <- N_pred/nrow(cellinfo)

			##1.1 Counting the proportion of cells with clonal CNVs
			res_TeaCNV<- do.call(rbind,lapply(1:length(CNVs_clonal),function(x,CNVs_clonal,clonalCN_res,N_pred,prop_pred,N_truth,prop_truth){
				clonalCNVs <- CNVs_clonal[x]
				cat(paste0(x,":",clonalCNVs,"\n"))
				cnv_index <- rownames(CNVs_clonal_bed) %in% clonalCNVs
				clonalCNVs_bed <- CNVs_clonal_bed[cnv_index,,drop=FALSE]

				cnv_value <- clonalCNVs_bed[,4:(4+i)]
				cloneTr <- names(clonalCNVs_bed[,4:(4+i)])[which(clonalCNVs_bed[,4:(4+i)]!=2)]
				prop_truth_cloneTr <- sum(prop_truth[cloneTr])
				clonalCNVs_bed$clone_truth <- paste(cloneTr,collapse=",")
		
				eva <- do.call(rbind,lapply(1:length(clonalCN_res),function(cl,clonalCN_res,clonalCNVs_bed,N_pred,prop_pred,N_truth,prop_truth,cloneTr){
					clonei <-names(clonalCN_res)[cl]
				  seg.dat <-  clonalCN_res[[clonei]]$seg.dat
	        bin_dat <- clonalCN_res[[clonei]]$input_BinSegRatio
	        bin_dat <- unique(bin_dat[,c("segName","SegMean")])
	        seg.dat <- merge(seg.dat,bin_dat,by="segName",all.x=TRUE)
	        seg.dat <- seg.dat[,c("chr","start","end","integerCN","segName")]
	        colnames(seg.dat)[4] = c(clonei)
	        seg.dat$start <- as.numeric(seg.dat$start)
	        seg.dat$end <- as.numeric(seg.dat$end)
					seg.dat$clone_pred <- clonei

	        intersect_res <- bt.intersect(a = clonalCNVs_bed, b = seg.dat, wa = TRUE, wb = TRUE)	
	        if(length(intersect_res)>0){
	        	colnames(intersect_res)<- c(colnames(clonalCNVs_bed),c("chr_pred","start_pred","end_pred","integerCN_pred","seg_pred","clone_pred"))
						# if(any(as.numeric(intersect_res$integerCN_pred)%in% as.numeric(intersect_res$integerCN))){
						if(any(as.numeric(intersect_res$integerCN_pred) !=2 )){	
							cnvFreq_results <- data.frame(
			      	clone_pred = clonei,
			      	Ncell_pred = N_pred[[clonei]],
			      	prop_pred = prop_pred[[clonei]],
			      	clone_true = clonalCNVs_bed[["clone_truth"]],
			      	Ncell_truth = sum(N_truth[cloneTr]),
			      	prop_truth = prop_truth_cloneTr,
			      	binID = clonalCNVs)

						}else{
							cnvFreq_results <- c()
						}
	        }else{
							cnvFreq_results <- c()
						}

		     
		      return(cnvFreq_results)
				},clonalCN_res,clonalCNVs_bed,N_pred,prop_pred,N_truth,prop_truth,cloneTr))

				if(!is.null(eva)){
					eva_res <- as.list(rep(NA, ncol(eva)))
					names(eva_res) <- names(eva)
					for (i in seq_along(eva)) {
					  if (i == 1) {
					    eva_res[[i]] <- paste(eva[[i]], collapse = ",")
					  } else if (i %in% c(2,3)){
					  	eva_res[[i]] <- sum(as.numeric(eva[[i]]), na.rm = TRUE)
					  }else if(i %in% c(4,5,6,7)) {
					    eva_res[[i]] <- unique(eva[[i]])
					  } 
					}
					eva_res <- as.data.frame(eva_res, stringsAsFactors = FALSE)
					eva_res$prop_dist <- abs(as.numeric(eva_res$prop_pred)-as.numeric(eva_res$prop_truth))
					}else{
						eva_res <- c()
					}

				return(eva_res)
			},CNVs_clonal,clonalCN_res,N_pred,prop_pred,N_truth,prop_truth))

			res_TeaCNV$method <- "TeaCNV"
			res_TeaCNV$nth <- nth
	  	res_TeaCNV$group <- group
	  	res_TeaCNV$size <- size

		  # res_clonal_all <- rbind(res_clonal_all,res_TeaCNV)

		  write.table(res_TeaCNV,file = res_clonal_file_tea,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_clonal_file_tea),
		    append = TRUE
		  )

		  ##1.2 Evaluating rare CNVs metrics
		  res_rare <- c()
		  for(clonei in names(clonalCN_res)){
		        seg.dat <-  clonalCN_res[[clonei]]$seg.dat
		        bin_dat <- clonalCN_res[[clonei]]$input_BinSegRatio
		        bin_dat <- unique(bin_dat[,c("segName","SegMean")])
		        seg.dat <- merge(seg.dat,bin_dat,by="segName",all.x=TRUE)
		        seg.dat <- seg.dat[,c("chr","start","end","integerCN","segName")]
		        colnames(seg.dat)[4] = c(clonei)
		        seg.dat$start <- as.numeric(seg.dat$start)
		        seg.dat$end <- as.numeric(seg.dat$end)
		        #align TeaCNV to the same bin size as DNA
		        TeaCNV <- align_Grange2bin(bed.query=bed.ref,bed.subject=seg.dat)
		        clonalCNA <- TeaCNV[,clonei]
		        ###Rare CNV events identification accuracy
		        bin_dat_truth <- GT_100k[[i+1]]
		        if(i==1){bin_dat_truth <- GT_100k[[i]]}
						bin_dat_truth$truth <- 0
						bin_dat_truth$truth[rareCNVRegion_index] <- 1

		        # bin_dat_truth$pred <- clonalCNA
		        bin_dat_truth$pred <- ifelse(clonalCNA!=2 & !is.na(clonalCNA), 1, 0)
		        bin_dat_truth$pred[!rareCNVRegion_index] <- 0

		        bin_dat_truth<- bin_dat_truth[!is.na(bin_dat_truth$integerCN),,drop=F]

		        evares_ls <- evaluate_Binary(bin_dat_truth[["truth"]],bin_dat_truth[["pred"]])

		        evares2 <- data.frame(clone_pred=clonei,clone_truth=j,
               	  	accuracy=as.numeric(evares_ls$accuracy),
	                precision=as.numeric(evares_ls$precision),
	                recall=as.numeric(evares_ls$recall),
	                F1=as.numeric(evares_ls$f1))
		        
		        res_rare <- rbind(res_rare,evares2)

		  }
			if(all(!is.na(res_rare$F1))){
			  res_rare <- res_rare[which.max(res_rare$F1),,drop=F]
			}else{
				res_rare <- res_rare[which.max(res_rare$accuracy),,drop=F]
			} 
	  	res_rare$nth <- nth
	  	res_rare$group <- group
	  	res_rare$size <- size
		  res_rare$method <- "TeaCNV"

		   write.table(res_rare,file = res_rare_file,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_rare_file),
		    append = TRUE
		  )
		}
  }
}





####END





