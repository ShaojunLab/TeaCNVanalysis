###模拟数据评估三种方法
suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(Matrix)
    library(data.table)
    library(magrittr)
    library(Signac)
    library(parallel)
    library(stringr)
    library(logger)
    library(purrr)
    library(patchwork)
    library(glue)
    library(ggridges)
    library(EnsDb.Hsapiens.v86)
    library(TeaCNV)
    library(reshape2)
	library(mclust)
	library(tidyr)
	library(bedtoolsr)
})

options(expressions=10000)
script_path <- "./github/TeaCNV"
setwd(script_path)
source("./ana/funs_evaluation.r")
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

source("./github/TeaCNV/simulation/funs_evaluation_SimData.R")
outdir <- "./TeaCNV/simulation/diffCNVEva";if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
cnvdir <- "./TeaCNV/simulation/TeaCNVres_diffCNVpct"
setwd(outdir)
outdir_data <- paste0(outdir,"/data");if(!file.exists(outdir_data)){dir.create(outdir_data,recursive=T)}
outdir_data_TeaCNV <- paste0(outdir_data,"/TeaCNV");if(!file.exists(outdir_data_TeaCNV)){dir.create(outdir_data_TeaCNV,recursive=T)}


simDatadir <- "./TeaCNV/simulation/simData_diffCNVpct"
cnvProp_max <- c(0.3,0.4,0.6,0.8,0.9,1)
clonal_groups=c("Monoclonal","Biclonal","Triclonal","Tetraclonal")

###------------------###
### Part 1
###Global evaluation: precision, recall, and F1-score of CNVs event identification; and overall accuracy of breakpoint identification;
###------------------###

###(1) TeaCNV results
res_bACC_all <- c()
for(i in 1:length(cnvProp_max)){
	pct <- cnvProp_max[i]

	trFile <-paste0(outdir_data,"/list_ground_truth_100kb_CNVpct",pct,".rds")
	if(file.exists(trFile)){ GT_100k <- readRDS(trFile)}else{
		ground_truth <- read.csv(file.path(simDatadir,paste0("GroundTruth_CNVFrac",pct,".csv")))
		clones <- unique(ground_truth$clone)
		GT_100k <- lapply(1:length(clones),function(x,ground_truth){
			cl <- clones[x]
			data_sub <- ground_truth %>%
			dplyr::filter(clone == cl)%>%as.data.frame() %>% 
			column_to_rownames("binID")  %>% 
			dplyr::select("chr"=Chromosome,"start"=Start,"end"=End, segName ,integerCN)

			align <- align_Grange2bin(bed.query=bed.ref,bed.subject=data_sub)
			return(align)
			},ground_truth)
		names(GT_100k) <- clones
		saveRDS(GT_100k,paste0(outdir_data,"/list_ground_truth_100kb_CNVpct",pct,".rds"))

	}

	for(j in 1:4){
		group <- clonal_groups[j]
		GT_100k_j <- GT_100k[1:j] ##ground truth

		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",group,"-CNVpct",pct,"-",nth,"..."))
			cat("\n")
			outdir_obs2 <- paste0(cnvdir,"/CNVpct",pct,"/",group,"/",nth)


			clonalCNApath <- file.path(outdir_obs2,"final.CNVres.rds")
			if(!file.exists(clonalCNApath)){next}
			# cellinfo <- read.csv(file.path(outdir_obs2,"subCluster/cell_info_subClust.csv"),row.names=1)
			# ##load clonal CNAs
	    	outres <- readRDS(clonalCNApath)

	    	clonalCN_res <- outres$clonalest
		    clonalCNA_bin <- c()
		    res_bACC <- c()
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
		        clonalCNA_bin <- cbind(clonalCNA_bin,clonalCNA)

		         ###CNV profile accuracy
			    eva <- do.call(rbind,lapply(1:length(GT_100k_j),function(x,GT_100k_j,clonalCNA,clonei){
	              bin_dat_truthi <- GT_100k_j[[x]]
	              bin_dat_truthi$pred <- as.integer(clonalCNA)
	               bin_dat_truthi<- bin_dat_truthi[!is.na(bin_dat_truthi$integerCN),,drop=F]

               	  evares_ls <- evaluateCNVPerformance(bin_dat_truthi,truth_col = "integerCN",
                                    pred_col = "pred",
                                  normal_value = 2,
                                  evaluate_CNVregion = TRUE,
                                  evaluate_CNVstates = TRUE)
               	  evares2 <- data.frame(clone_pred=clonei,clone_truth=names(GT_100k_j)[x],
               	  	accuracy=as.numeric(evares_ls$CNV_region_detection$accuracy),
	                precision=as.numeric(evares_ls$CNV_region_detection$precision),
	                recall=as.numeric(evares_ls$CNV_region_detection$recall),
	                F1=as.numeric(evares_ls$CNV_region_detection$f1),
	                accuracy.CNVvalue=as.numeric(evares_ls$CNV_states_detection$accuracy),
	                precision.CNVvalue=as.numeric(evares_ls$CNV_states_detection$precision),
	                recall.CNVvalue=as.numeric(evares_ls$CNV_states_detection$recall),
	                F1.CNVvalue=as.numeric(evares_ls$CNV_states_detection$f1),
	                MAE.CNVvalue=as.numeric(evares_ls$CNV_states_detection$MAE),
	                RMSE.CNVvalue = as.numeric(evares_ls$CNV_states_detection$RMSE))
	              return(evares2)
	            },GT_100k_j,clonalCNA,clonei))
			    RMSE_mean <- median(eva$RMSE.CNVvalue,na.rm=TRUE)
		      if(all(!is.na(eva$F1))){
            eva <- eva[which.max(eva$F1),,drop=F]
        	}else{
        		eva <- eva[which.max(eva$accuracy),,drop=F]
        	}
        	eva$RMSE.CNVvalue <- RMSE_mean
            	

        	##3. break point distance
        	cnv_dist <- c()
        	 for(x in 1:length(GT_100k_j)){
	        	 	bin_dat_truthi <- GT_100k_j[[x]]
						  seg_truth <-  unique(bin_dat_truthi[,c("segName","integerCN")])
						  seg_truth <- seg_truth %>%dplyr::filter(integerCN!=2)
						  gr_truth <- strings2bed(seg_truth$segName)
						  seg_pred <- seg.dat %>%
						  	dplyr::filter(.data[[clonei]]!=2)%>%
						  	select(chr,start,end)
						  common_chr <- union(unique(gr_truth[,1]),unique(seg_pred[,1]))	
						  #fill missing chr
						  chr_miss_tr <- common_chr[!common_chr %in% unique(gr_truth[,1])]
						  if(length(chr_miss_tr)>0){
						  	bed_miss <- hgpos[hgpos$chr %in% chr_miss_tr,1:3]
						  	colnames(bed_miss) <- colnames(gr_truth)
						  	gr_truth <- rbind(gr_truth,bed_miss)
						  }
						   chr_miss_pr <- common_chr[!common_chr %in% unique(seg_pred[,1])]
						  if(length(chr_miss_pr)>0){
						  	bed_miss_pr <- hgpos[hgpos$chr %in% chr_miss_pr,1:3]
						  	colnames(bed_miss_pr) <- colnames(seg_pred)
						  	seg_pred <- rbind(seg_pred,bed_miss_pr)
						  }


						  intersect_res <- bt.intersect(a = gr_truth, b = seg_pred, wa = TRUE, wb = TRUE)	
						  colnames(intersect_res)<- c("chrom", "start", "end", "chrom_pr", "start_pr", "end_pr")
						  intersect_res$segName_pred <- paste(intersect_res[,4],intersect_res[,5],intersect_res[,6],sep="_")
						  intersect_res <- intersect_res %>%
						  	mutate(dist_start=abs(start_pr-start),
						  		dist_end=abs(end_pr-end))%>%as.data.frame()
						  dist_res <- intersect_res %>%
						  	group_by(segName_pred)%>%
						  	summarise(dist_start = min(dist_start, na.rm = TRUE),
						  		dist_end = min(dist_end, na.rm = TRUE))
						  dist_res$dist_mean = rowMeans(dist_res[,c('dist_start','dist_end')])
						  
							meandis  <-median(dist_res$dist_mean,na.rm=TRUE)/1e6 ##MB
							cnv_dist <- c(cnv_dist,meandis)
            }

          eva$bp_dist <- min(cnv_dist)
          res_bACC <- rbind(res_bACC,eva)

		    }
		    colnames(clonalCNA_bin) <- names(clonalCN_res)
		    rownames(clonalCNA_bin) <- paste(bed.ref[,1],bed.ref[,2],bed.ref[,3],sep="_")
		   
		    write.csv(clonalCNA_bin,paste0(outdir_data_TeaCNV,"/teaCNV_clonalCNA_bin_CNVpct",pct,"_",group,"_",nth,".csv"))
	    	
		  	res_bACC$nth <- nth
		  	res_bACC$group <- group
		  	res_bACC$CNVFrac <- pct
		  
		  	res_bACC_all <- rbind(res_bACC_all,res_bACC)
		}
	}
}

write.csv(res_bACC_all,paste0(outdir_data,"/evaRes_TeaCNV.csv"),row.names=FALSE)



### ##(2) epiAneuFinder 
epiAneu_path <- "./TeaCNV/simulation/epiAneufinder_diffCNVpct/"
outdir_data_epiAneu <- paste0(outdir_data,"/epiAneu");if(!file.exists(outdir_data_epiAneu)){dir.create(outdir_data_epiAneu,recursive=T)}

res_bACC_all_epiAneu <- c()

for(i in 1:length(cnvProp_max)){
	pct <- cnvProp_max[i]

	trFile <-paste0(outdir_data,"/list_ground_truth_100kb_CNVpct",pct,".rds")
	if(file.exists(trFile)){ GT_100k <- readRDS(trFile)}else{stop(paste0(trFile," NOT exists!"))}

	for(j in 1:4){
		group <- clonal_groups[j]
		GT_100k_j <- GT_100k[1:j] 

		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",group,"-CNVpct",pct,"-",nth,"..."))
			cat("\n")
			epiAneu_path2 <- paste0(epiAneu_path,"/CNVpct",pct,"/",group,"/",nth,"/epiAneufinder_results")

			cna_epiAneu <- read.table(paste0(epiAneu_path2, '/results_table.tsv'),header=T,row.names=1)
		    rownames(cna_epiAneu) <- paste(cna_epiAneu[,1],cna_epiAneu[,2],cna_epiAneu[,3],sep="_")
		    cells_epiAneu <- colnames(cna_epiAneu)[4:ncol(cna_epiAneu)]
		    cells_epiAneu <- gsub("cell.","",cells_epiAneu)
		    colnames(cna_epiAneu)[4:ncol(cna_epiAneu)] <- cells_epiAneu
		    normalcell <- cells_epiAneu[grepl("sample1_",cells_epiAneu)]
		    keepIndx <- which(!colnames(cna_epiAneu) %in% normalcell)
		    cna_epiAneu <- cna_epiAneu[,c(keepIndx)]

		    res_bACC_epi <- c()
		    Psets  <- seq(0.1,1,0.1)
		    for(prop in Psets){
	        cna_epiAneu_bulk <- get_epiAneufinder_bulk(cna_epiAneu,bed.ref,cnvCellProp.min=prop)
	        colnames(cna_epiAneu_bulk)[4]<- "epiAneuRatio"
	        # write.csv(cna_epiAneu_bulk,paste0(outdir_data_epiAneu,"/epiAneufinder_100kb_CNVpct",pct,"_",group,"_n",nth,"_",prop,".csv"))
	   			eva_epi <- do.call(rbind,lapply(1:length(GT_100k_j),function(x,GT_100k_j,cna_epiAneu_bulk){
	            bin_dat_truthi <- GT_100k_j[[x]]
	            ##epiAneuCN==1 : "neu"
	            bin_dat_truthi$pred <- as.integer(cna_epiAneu_bulk$epiAneuCN+1) ##make 2 is neu_value
							# bin_dat_truthi$epiAneuRatio <- as.integer(cna_epiAneu_bulk$epiAneuRatio) ##make 2 is neu_value
	            bin_dat_truthi<- bin_dat_truthi[!is.na(bin_dat_truthi$integerCN),,drop=F]
	         
	            truth <- bin_dat_truthi[["integerCN"]]
 							pred  <- bin_dat_truthi[["pred"]]
 							bin_truth <- ifelse(truth != 2, 1, 0)
      				bin_pred  <- ifelse(pred != 2 & !is.na(pred), 1, 0)

 							TP <- sum(bin_truth == 1 & bin_pred == 1)
					    FP <- sum(bin_truth != 1 & bin_pred == 1)
					    FN <- sum(bin_truth == 1 & bin_pred != 1)
					    TN <- sum(bin_truth != 1 & bin_pred != 1)
	
	       		  evares_ls <- evaluateCNVPerformance(bin_dat_truthi,truth_col = "integerCN",
                                    pred_col = "pred",
                                  normal_value = 2,
                                  evaluate_CNVregion = TRUE,
                                  evaluate_CNVstates = TRUE)
               	  evares2 <- data.frame(method_prop=prop,clone_truth=names(GT_100k_j)[x],
               	  	TP=TP,FP=FP,FN=FN,TN=TN,
               	  	accuracy=as.numeric(evares_ls$CNV_region_detection$accuracy),
	                precision=as.numeric(evares_ls$CNV_region_detection$precision),
	                recall=as.numeric(evares_ls$CNV_region_detection$recall),
	                F1=as.numeric(evares_ls$CNV_region_detection$f1),
	                RMSE.CNVvalue=as.numeric(evares_ls$CNV_states_detection$RMSE)  
	                )
	            return(evares2)
	        },GT_100k_j,cna_epiAneu_bulk))

	   			RMSE_mean <- median(eva_epi$RMSE.CNVvalue,na.rm=TRUE)
          if(all(!is.na(eva_epi$F1))){
          	eva_epi <- eva_epi[which.max(eva_epi$F1),,drop=F]
        	}else{
        		eva_epi <- eva_epi[which.max(eva_epi$accuracy),,drop=F]
        	}
 					eva_epi$RMSE.CNVvalue <- RMSE_mean
	     	
		   		##3. break point distance
	      	cnv_dist <- c()
	      	for(x in 1:length(GT_100k_j)){
	          bin_dat_truthi <- GT_100k_j[[x]]
						seg_truth <-  unique(bin_dat_truthi[,c("segName","integerCN")])
						seg_truth <- seg_truth %>%dplyr::filter(integerCN!=2)
						gr_truth <- strings2bed(seg_truth$segName)

						seg_pred <- cna_epiAneu_bulk %>%
							dplyr::filter(epiAneuCN!=1,!is.na(epiAneuCN))%>%
						  arrange(chr, start) %>%
						  group_by(chr) %>%
						  mutate(
						    group_id = cumsum(
						      # 如果epiAneuCN变化 或 start不是上一行end，则开启新组
						      (epiAneuCN != dplyr::lag(epiAneuCN, default = epiAneuCN[1])) |
						      (start != dplyr::lag(end, default = start[1]))
						    )
						  ) %>%
						  group_by(chr, group_id) %>%
						  summarise(
						    start = dplyr::first(start),
						    end = dplyr::last(end),
						    .groups = "drop"
						  ) %>%
						  select(chr, start, end)%>%as.data.frame()

						
					  common_chr <- union(unique(gr_truth[,1]),unique(seg_pred[,1]))	
					  #fill missing chr
					  chr_miss_tr <- common_chr[!common_chr %in% unique(gr_truth[,1])]
					  if(length(chr_miss_tr)>0){
					  	bed_miss <- hgpos[hgpos$chr %in% chr_miss_tr,1:3]
					  	colnames(bed_miss) <- colnames(gr_truth)
					  	gr_truth <- rbind(gr_truth,bed_miss)
					  }
					   chr_miss_pr <- common_chr[!common_chr %in% unique(seg_pred[,1])]
					  if(length(chr_miss_pr)>0){
					  	bed_miss_pr <- hgpos[hgpos$chr %in% chr_miss_pr,1:3]
					  	colnames(bed_miss_pr) <- colnames(seg_pred)
					  	seg_pred <- rbind(seg_pred,bed_miss_pr)
					  }
					  intersect_res <- bt.intersect(a = gr_truth, b = seg_pred, wa = TRUE, wb = TRUE)	
					  colnames(intersect_res)<- c("chrom", "start", "end", "chrom_pr", "start_pr", "end_pr")
					  intersect_res$segName_pred <- paste(intersect_res[,4],intersect_res[,5],intersect_res[,6],sep="_")
					  intersect_res <- intersect_res %>%
					  	mutate(dist_start=abs(start_pr-start),
					  		dist_end=abs(end_pr-end))%>%as.data.frame()
					  dist_res <- intersect_res %>%
					  	group_by(segName_pred)%>%
					  	summarise(dist_start = min(dist_start, na.rm = TRUE),
					  		dist_end = min(dist_end, na.rm = TRUE))
					  dist_res$dist_mean = rowMeans(dist_res[,c('dist_start','dist_end')])
					  
						meandis  <-median(dist_res$dist_mean,na.rm=TRUE)/1e6 ##MB
						cnv_dist <- c(cnv_dist,meandis)

          }

	        eva_epi$bp_dist <- min(cnv_dist)
	        res_bACC_epi <- rbind(res_bACC_epi,eva_epi)
		    }
		    
		    res_bACC_epi$nth <- nth
		  	res_bACC_epi$group <- group
		  	res_bACC_epi$CNVFrac <- pct
		  
		  	res_bACC_all_epiAneu <- rbind(res_bACC_all_epiAneu,res_bACC_epi)

		}
	}
}
write.csv(res_bACC_all_epiAneu,paste0(outdir_data,"/evaRes_epiAneu.csv"),row.names=FALSE)



###(3)CopyscAT 
CopyscAT_path <- "./TeaCNV/simulation/CopyscAT_diffCNVpct/"
outdir_data_CopyscAT <- paste0(outdir_data,"/CopyscAT");if(!file.exists(outdir_data_CopyscAT)){dir.create(outdir_data_CopyscAT,recursive=T)}
res_bACC_all_CopyscAT <- c()
for(i in 1:length(cnvProp_max)){
	pct <- cnvProp_max[i]

	trFile <-paste0(outdir_data,"/list_ground_truth_100kb_CNVpct",pct,".rds")
	if(file.exists(trFile)){ GT_100k <- readRDS(trFile)}else{stop(paste0(trFile," NOT exists!"))}

	for(j in 1:4){
		group <- clonal_groups[j]
		GT_100k_j <- GT_100k[1:j] 

		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",group,"-CNVpct",pct,"-",nth,"..."))
			cat("\n")

			#(3) load CopyscAT
			copy_path <- paste0(CopyscAT_path,"/CNVpct",pct,"/",group,"/",nth)
			cnaPath <- file.path(copy_path,paste(group,pct,nth,"cnv","scores.csv",sep="_"))
			cna_copyscAT <- get_copyscAT_sc(cnaPath,cytoBandFile)
			cells_copyscAT <- colnames(cna_copyscAT)[4:ncol(cna_copyscAT)]
    		cells_copyscAT <- gsub("cell.","",cells_copyscAT)

    		res_bACC_copyscAT <- c()
    		Psets  <- seq(0.1,1,0.1)
		    for(prop in Psets){
		        cna_copyscAT_bulk_align <- get_copyscAT_bulk(cna_copyscAT,bed.ref,cnvCellProp.min=prop)
		        # write.csv(cna_copyscAT_bulk_align,paste0(outdir_data_CopyscAT,"/copyscAT_100kb_CNVpct",pct,"_",group,"_n",nth,"_",prop,".csv"))
		    
		        eva_copy <- do.call(rbind,lapply(1:length(GT_100k_j),function(x,GT_100k_j,cna_copyscAT_bulk_align){
		            bin_dat_truthi <- GT_100k_j[[x]]
		            bin_dat_truthi$pred <- as.integer(cna_copyscAT_bulk_align$copyscAT.CN)
		            bin_dat_truthi$copyscAT_ratio <- cna_copyscAT_bulk_align$copyscAT_ratio
		            bin_dat_truthi<- bin_dat_truthi[!is.na(bin_dat_truthi$integerCN),,drop=F]

	   	            truth <- bin_dat_truthi[["integerCN"]]
					pred  <- bin_dat_truthi[["pred"]]
					bin_truth <- ifelse(truth != 2, 1, 0)
	  				bin_pred  <- ifelse(pred != 2 & !is.na(pred), 1, 0)

					TP <- sum(bin_truth == 1 & bin_pred == 1)
				    FP <- sum(bin_truth != 1 & bin_pred == 1)
				    FN <- sum(bin_truth == 1 & bin_pred != 1)
				    TN <- sum(bin_truth != 1 & bin_pred != 1)
						   
	       		  	evares_ls <- evaluateCNVPerformance(bin_dat_truthi,truth_col = "integerCN",
	                                pred_col = "pred",
	                              normal_value = 2,
	                              evaluate_CNVregion = TRUE,
	                              evaluate_CNVstates = TRUE)
		           	evares2 <- data.frame(method_prop=prop,clone_truth=names(GT_100k_j)[x],
		           	  	TP=TP,FP=FP,FN=FN,TN=TN,
		           	  	accuracy=as.numeric(evares_ls$CNV_region_detection$accuracy),
		                precision=as.numeric(evares_ls$CNV_region_detection$precision),
		                recall=as.numeric(evares_ls$CNV_region_detection$recall),
		                F1=as.numeric(evares_ls$CNV_region_detection$f1),
		                RMSE.CNVvalue = as.numeric(evares_ls$CNV_states_detection$RMSE)
		              )
		            return(evares2)
		        },GT_100k_j,cna_copyscAT_bulk_align))

		        RMSE_mean <- median(eva_copy$RMSE.CNVvalue,na.rm=TRUE)

		        if(all(!is.na(eva_copy$F1))){
	            	eva_copy <- eva_copy[which.max(eva_copy$F1),,drop=F]
	          	}else{
	          		eva_copy <- eva_copy[which.max(eva_copy$accuracy),,drop=F]
	          	}
		   		eva_copy$RMSE.CNVvalue <- RMSE_mean
	     	 ##3. break point distance
        	 cnv_dist <- c()
        	 for(x in 1:length(GT_100k_j)){
        	 	  bin_dat_truthi <- GT_100k_j[[x]]
							seg_truth <-  unique(bin_dat_truthi[,c("segName","integerCN")])
							seg_truth <- seg_truth %>%dplyr::filter(integerCN!=2)
							gr_truth <- strings2bed(seg_truth$segName)

							seg_pred <- cna_copyscAT_bulk_align %>%
								dplyr::filter(copyscAT.CN!=2,!is.na(copyscAT.CN))%>%
							  arrange(chr, start) %>%
							  group_by(chr) %>%
							  mutate(
							    group_id = cumsum(
							      (copyscAT.CN != dplyr::lag(copyscAT.CN, default = copyscAT.CN[1])) |
							      (start != dplyr::lag(end, default = start[1]))
							    )
							  ) %>%
							  group_by(chr, group_id) %>%
							  summarise(
							    start = dplyr::first(start),
							    end = dplyr::last(end),
							    .groups = "drop"
							  ) %>%
							  select(chr, start, end)%>%as.data.frame()
							common_chr <- union(unique(gr_truth[,1]),unique(seg_pred[,1]))	
						  #fill missing chr
						  chr_miss_tr <- common_chr[!common_chr %in% unique(gr_truth[,1])]
						  if(length(chr_miss_tr)>0){
						  	bed_miss <- hgpos[hgpos$chr %in% chr_miss_tr,1:3]
						  	colnames(bed_miss) <- colnames(gr_truth)
						  	gr_truth <- rbind(gr_truth,bed_miss)
						  }
						   chr_miss_pr <- common_chr[!common_chr %in% unique(seg_pred[,1])]
						  if(length(chr_miss_pr)>0){
						  	bed_miss_pr <- hgpos[hgpos$chr %in% chr_miss_pr,1:3]
						  	colnames(bed_miss_pr) <- colnames(seg_pred)
						  	seg_pred <- rbind(seg_pred,bed_miss_pr)
						  }
						  intersect_res <- bt.intersect(a = gr_truth, b = seg_pred, wa = TRUE, wb = TRUE)	
						  colnames(intersect_res)<- c("chrom", "start", "end", "chrom_pr", "start_pr", "end_pr")
						  intersect_res$segName_pred <- paste(intersect_res[,4],intersect_res[,5],intersect_res[,6],sep="_")
						  intersect_res <- intersect_res %>%
						  	mutate(dist_start=abs(start_pr-start),
						  		dist_end=abs(end_pr-end))%>%as.data.frame()
						  dist_res <- intersect_res %>%
						  	group_by(segName_pred)%>%
						  	summarise(dist_start = min(dist_start, na.rm = TRUE),
						  		dist_end = min(dist_end, na.rm = TRUE))
						  dist_res$dist_mean = rowMeans(dist_res[,c('dist_start','dist_end')])
						  
							meandis  <-median(dist_res$dist_mean,na.rm=TRUE)/1e6 ##MB
						cnv_dist <- c(cnv_dist,meandis)
	      	}

	        eva_copy$bp_dist <- min(cnv_dist)
	        res_bACC_copyscAT <- rbind(res_bACC_copyscAT,eva_copy)
		    }
		    
		    res_bACC_copyscAT$nth <- nth
		  	res_bACC_copyscAT$group <- group
		  	res_bACC_copyscAT$CNVFrac <- pct
		  
		  	res_bACC_all_CopyscAT <- rbind(res_bACC_all_CopyscAT,res_bACC_copyscAT)


		}
	}
}
write.csv(res_bACC_all_CopyscAT,paste0(outdir_data,"/evaRes_CopyscAT.csv"),row.names=FALSE)






###------------------###
### Part 2 ##Statistical clonal composition（BACC）
###Accuracy of subclonal event identification: Accuracy of CNVs precision, recall, and F1-score for the specified subclonal CNVs
###Comparison of the three methods by calculating the absolute difference between the inferred cell proportion and the true proportion
###------------------###
outdir_data <- paste0(outdir,"/data")
cnvdir_tea <- "./TeaCNV/simulation/TeaCNVres_diffCNVpct"
CopyscAT_path <- "./TeaCNV/simulation/CopyscAT_diffCNVpct/"
epiAneu_path <- "./TeaCNV/simulation/epiAneufinder_diffCNVpct/"
cnvProp_max <- c(0.3,0.4,0.6)
clonal_groups=c("Monoclonal","Biclonal","Triclonal","Tetraclonal")
res_perform_clonal <- c()
res_perform_rare <- c()

for(i in 1:length(cnvProp_max)){
	pct <- cnvProp_max[i]

	for(j in 1:4){
			group <- clonal_groups[j]
		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",group,"-CNVpct",pct,"-",nth,"..."))
			cat("\n")

			#load TeaCNV
			cnvdir_tea2 <- paste0(cnvdir_tea,"/CNVpct",pct,"/",group,"/",nth)
			clonalCNApath <- file.path(cnvdir_tea2,"final.CNVres.rds")
			if(!file.exists(clonalCNApath)){print("File dose not exist!");next}
			outres <- readRDS(clonalCNApath)
			clonalCN_res <- outres$clonalest

			cellinfo <- outres$cellinfo
			cellinfo$truthClone <- sapply(strsplit(cellinfo$cellname,"_"),"[",2)
			cellinfo$truthClone <- gsub("C","",cellinfo$truthClone)
			N_truth <- table(cellinfo$truthClone)
			prop_truth <- N_truth/nrow(cellinfo)
			prop_pred <- table(cellinfo$clone)/nrow(cellinfo)
			truth <- as.character( cellinfo[["truthClone"]])
      pred <- as.character( cellinfo[["clone"]])

      res_clone <- evaluate_MultiClassification(truth,pred)
 			evares_clone <- data.frame(group=group,size=2000,CNVFrac=pct,nth=nth,
				accuracy=as.numeric(res_clone$accuracy),
				precision=as.numeric(res_clone$precision),
				recall=as.numeric(res_clone$recall),
				F1=as.numeric(res_clone$f1),
				BACC=as.numeric(res_clone$balanced_accuracy)
				)
 			res_perform_clonal <- rbind(res_perform_clonal,evares_clone)
 			if(j>1){
	      res_rareclone <- evaluateRareCloneAccuracy(cellinfo,truth_col="truthClone",pred_col="clone",rare_prop=0.12)
	      evares_rare <- data.frame(group=group,size=2000,CNVFrac=pct,nth=nth,
					accuracy=as.numeric(res_rareclone$accuracy),
					precision=as.numeric(res_rareclone$precision),
					recall=as.numeric(res_rareclone$recall),
					F1=as.numeric(res_rareclone$f1),
					BACC=as.numeric(res_rareclone$balanced_accuracy)
					)
				res_perform_rare <- rbind(res_perform_rare,evares_rare)
			}
		}
	}
}

write.csv(res_perform_clonal,paste0(outdir_data,"/ "),row.names=FALSE)
write.csv(res_perform_rare,paste0(outdir_data,"/evaRes_TeaCNV_RareClone_assignment.csv"),row.names=FALSE)


###------------------###
### Part 3: Accuracy of rare clonal CNV event identification: Precision, recall, and F1-score for the identified rare clonal CNV events.
### Comparison of the three methods by calculating the absolute difference between the inferred cell proportion and the true proportion for each clonal CNV event.
###------------------###
outdir_data <- paste0(outdir,"/data")
cnvdir_tea <- "./TeaCNV/simulation/TeaCNVres_diffCNVpct"
CopyscAT_path <- "./TeaCNV/simulation/CopyscAT_diffCNVpct/"
epiAneu_path <- "./TeaCNV/simulation/epiAneufinder_diffCNVpct/"
cnvProp_max <- c(0.3,0.4,0.6)
clonal_groups=c("Monoclonal","Biclonal","Triclonal","Tetraclonal")

res_clonal_all_file <- file.path(outdir_data,"res_clonalCNV_freq_TeaCNV.csv")
if (file.exists(res_clonal_all_file)) file.remove(res_clonal_all_file)
res_rare_all_file <- file.path(outdir_data,"res_rareCNVevents_identify_performance_TeaCNV.csv")
if (file.exists(res_rare_all_file)) file.remove(res_rare_all_file)

res_rare_all_file_epi <- file.path(outdir_data,"res_rareCNVevents_identify_performance_epiAneuFinder.csv")
if (file.exists(res_rare_all_file_epi)) file.remove(res_rare_all_file_epi)
res_clonalCNV_freq_file_epi <- file.path(outdir_data,"res_clonalCNV_freq_epiAneuFinder.csv")
if (file.exists(res_clonalCNV_freq_file_epi)) file.remove(res_clonalCNV_freq_file_epi)


res_rare_all_file_copy <- file.path(outdir_data,"res_rareCNVevents_identify_performance_copyscAT.csv")
if (file.exists(res_rare_all_file_copy)) file.remove(res_rare_all_file_copy)
res_clonalCNV_freq_file_copy <- file.path(outdir_data,"res_clonalCNV_freq_copyscAT.csv")
if (file.exists(res_clonalCNV_freq_file_copy)) file.remove(res_clonalCNV_freq_file_copy)


for(i in 1:length(cnvProp_max)){
	pct <- cnvProp_max[i]
	###Ground truth CNV
	trFile <-paste0(outdir_data,"/list_ground_truth_100kb_CNVpct",pct,".rds")
	if(file.exists(trFile)){ GT_100k <- readRDS(trFile)}else{stop(paste0(trFile," NOT exists!"))}
	trcnv_File <-paste0("./TeaCNV/simulation/simData_diffCNVpct/GroundTruth_CNVFrac",pct,".csv")
	if(file.exists(trcnv_File)){ GT <- read.csv(trcnv_File)}

	truth_wide <- GT %>%
		dplyr::select(binID,integerCN,clone)%>%
		dplyr::distinct()%>%
	  pivot_wider(
	    names_from = clone,
	    values_from = integerCN
	  )%>%
 	tibble::column_to_rownames(var = "binID")

	for(j in 2:4){
		group <- clonal_groups[j]

		truth_wide_j <- truth_wide[,1:j]

    
		
	 	CNVinRow <- apply(truth_wide_j,1,unique)
	 	cnvRow<- unlist(lapply(CNVinRow,function(x){length(x[!is.na(x)])>1}))
	 	truth_wide_clonalCNV <- truth_wide_j[cnvRow,,drop=F]
	 	CNVs_clonal <- rownames(truth_wide_clonalCNV)

	 	CNVs_clonal_bed <- strings2bed(CNVs_clonal)
		CNVs_clonal_bed <- cbind(CNVs_clonal_bed,truth_wide_clonalCNV)
		

	 	rareCNV_index <- apply(truth_wide_clonalCNV,1,function(x,j){
	 		(x[j]!=2) && (!x[j] %in% x[1:(j-1)])
	 		},j)
	 	CNVs_rare <- rownames(truth_wide_clonalCNV)[rareCNV_index]
	 	CNVs_rare_bed <- strings2bed(CNVs_rare)
	 	CNVs_rare_bed100 <- align_Grange2bin(bed.query=bed.ref,bed.subject=CNVs_rare_bed)
	 	rareCNVRegion_index <- !is.na(CNVs_rare_bed100$chr.subject)

		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",group,"-CNVpct",pct,"-",nth,"..."))
			cat("\n")


			##(1)TeaCNV
			load TeaCNV
			cnvdir_tea2 <- paste0(cnvdir_tea,"/CNVpct",pct,"/",group,"/",nth)
			clonalCNApath <- file.path(cnvdir_tea2,"final.CNVres.rds")
			if(!file.exists(clonalCNApath)){print("File dose not exist!");next}
			outres <- readRDS(clonalCNApath)
			clonalCN_res <- outres$clonalest

			cellinfo <- outres$cellinfo
			cellinfo$truthClone <- sapply(strsplit(cellinfo$cellname,"_"),"[",2)
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

				cnv_value <- clonalCNVs_bed[,4:(3+j)]
				cloneTr <- names(clonalCNVs_bed[,4:(3+j)])[which(clonalCNVs_bed[,4:(3+j)]!=2)]
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
	  	res_TeaCNV$CNVFrac <- pct

		  # res_clonal_all <- rbind(res_clonal_all,res_TeaCNV)

		  write.table(res_TeaCNV,file = res_clonal_all_file,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_clonal_all_file),
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
		        bin_dat_truth <- GT_100k[[j]]
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
	  	res_rare$CNVFrac <- pct
		  res_rare$method <- "TeaCNV"

		   write.table(res_rare,file = res_rare_all_file,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_rare_all_file),
		    append = TRUE
		  )


		  ###
		  ###(2)	epiAneufinder
		  cna_epiAneu_bulk_align <- cna_epiAneu_bulk_align %>%
        mutate(epiAneuCNstate = case_when(
            epiAneuCN ==2 ~ 'amp',
            epiAneuCN ==0 ~ 'del',
            is.na(epiAneuCN) ~ NA,
            T ~ 'neu'
        ))
       
			epiAneu_path2 <- paste0(epiAneu_path,"/CNVpct",pct,"/",group,"/",nth,"/epiAneufinder_results")
			cna_epiAneu <- read.table(paste0(epiAneu_path2, '/results_table.tsv'),header=T,row.names=1)
	    rownames(cna_epiAneu) <- paste(cna_epiAneu[,1],cna_epiAneu[,2],cna_epiAneu[,3],sep="_")
	    cells_epiAneu <- colnames(cna_epiAneu)[4:ncol(cna_epiAneu)]
	    cells_epiAneu <- gsub("cell.","",cells_epiAneu)
	    colnames(cna_epiAneu)[4:ncol(cna_epiAneu)] <- cells_epiAneu
	    normalcell <- cells_epiAneu[grepl("sample1_",cells_epiAneu)]
	    keepIndx <- which(!colnames(cna_epiAneu) %in% normalcell)
	    cna_epiAneu <- cna_epiAneu[,c(keepIndx)]

	    epiAneu_bed <- data.frame(chr=cna_epiAneu[,1],start=as.numeric(cna_epiAneu[,2]),end=as.numeric(cna_epiAneu[,3]),binID=rownames(cna_epiAneu))
	    cna_epiAneu_mt <- cna_epiAneu[,4:ncol(cna_epiAneu)]




	    # ###2.1 Counting the proportion of cells with clonal CNVs
	    obs_cells <- colnames(cna_epiAneu_mt)
	  	truthClone <- sapply(strsplit(obs_cells,"_"),"[",2)
	  	truthClone <- gsub("C","",truthClone)
	  	cellinfo<- data.frame(cellname=obs_cells,truthClone=truthClone)
	  	N_truth <- table(truthClone)
	  	prop_truth <- N_truth/sum(N_truth)

	    res<- do.call(rbind,lapply(1:length(CNVs_clonal),function(x,CNVs_clonal,N_truth,prop_truth,epiAneu_bed){
				clonalCNVs <- CNVs_clonal[x]
				cat(paste0(x,":",clonalCNVs,"\n"))
				cnv_index <- rownames(CNVs_clonal_bed) %in% clonalCNVs
				clonalCNVs_bed <- CNVs_clonal_bed[cnv_index,,drop=FALSE]

				cnv_value <- clonalCNVs_bed[,4:(3+j)]
				cloneTr <- names(clonalCNVs_bed[,4:(3+j)])[which(clonalCNVs_bed[,4:(3+j)]!=2)]
				prop_truth_cloneTr <- sum(prop_truth[cloneTr])
				clonalCNVs_bed$clone_truth <- paste(cloneTr,collapse=",")

				intersect_res <- bt.intersect(a = clonalCNVs_bed, b = epiAneu_bed, wa = TRUE, wb = TRUE)	
				if(length(intersect_res)>0){
					colnames(intersect_res)<- c(colnames(clonalCNVs_bed),c("chr_pred","start_pred","end_pred","binID_pred"))
					row_index <- which(rownames(cna_epiAneu_mt)%in%intersect_res$binID_pred)

					Npred <- rowSums(cna_epiAneu_mt!=1,na.rm=TRUE)
					Npred <- mean(Npred[row_index])
					prop_pred <- Npred/ncol(cna_epiAneu_mt)
					cnvFreq_results <- data.frame(
			      	Ncell_pred = Npred,
			      	prop_pred = prop_pred,
			      	clone_true = clonalCNVs_bed[["clone_truth"]],
			      	Ncell_truth = sum(N_truth[cloneTr]),
			      	prop_truth = prop_truth_cloneTr,
			      	binID = clonalCNVs)
				}else{
						cnvFreq_results <- c()
				}
			 return(cnvFreq_results)
			 },CNVs_clonal,N_truth,prop_truth,epiAneu_bed))
			res$prop_dist <- abs(as.numeric(res$prop_pred)-as.numeric(res$prop_truth))
			res$method <- "epiAneuFinder"
			res$nth <- nth
	  	res$group <- group
	  	res$CNVFrac <- pct

	  	 write.table(res,file = res_clonalCNV_freq_file_epi,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_clonalCNV_freq_file_epi),
		    append = TRUE
		  )






	    ##2.2 Evaluating rare CNVs metrics
      res_epi <- c()
	    Psets  <- seq(0.1,1,0.1)
	    for(prop in Psets){
        cna_epiAneu_bulk <- get_epiAneufinder_bulk(cna_epiAneu,bed.ref,cnvCellProp.min=prop)
        colnames(cna_epiAneu_bulk)[4]<- "epiAneuRatio"

        bin_dat_truth <- GT_100k[[j]]
        bin_dat_truth$truth <- 0
				bin_dat_truth$truth[rareCNVRegion_index] <- 1

		    ## epiAneuCN=1 is neu_value
		    bin_dat_truth$pred <- ifelse(cna_epiAneu_bulk$epiAneuCN!=1 & !is.na(cna_epiAneu_bulk$epiAneuCN), 1, 0)
		    bin_dat_truth$pred[!rareCNVRegion_index] <- 0



		    bin_dat_truth<- bin_dat_truth[!is.na(bin_dat_truth$integerCN),,drop=F]
	         
				evares_ls <- evaluate_Binary(bin_dat_truth[["truth"]],bin_dat_truth[["pred"]])

        evares2 <- data.frame(method_prop=prop,
        				clone_truth=j,
           	  	accuracy=as.numeric(evares_ls$accuracy),
              precision=as.numeric(evares_ls$precision),
              recall=as.numeric(evares_ls$recall),
              F1=as.numeric(evares_ls$f1))
        
        res_epi <- rbind(res_epi,evares2)

      }

      res_epi$nth <- nth
	  	res_epi$group <- group
	  	res_epi$CNVFrac <- pct
	  	res_epi$method <- "epiAneuFinder"

	  	 write.table(res_epi,file = res_rare_all_file_epi,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_rare_all_file_epi),
		    append = TRUE
		  )





			####
			####(3)copyscAT
	
	  	copy_path <- paste0(CopyscAT_path,"/CNVpct",pct,"/",group,"/",nth)
			cnaPath <- file.path(copy_path,paste(group,pct,nth,"cnv","scores.csv",sep="_"))
			cna_copyscAT <- get_copyscAT_sc(cnaPath,cytoBandFile)
			cells_copyscAT <- colnames(cna_copyscAT)[4:ncol(cna_copyscAT)]
    	cells_copyscAT <- gsub("cell.","",cells_copyscAT)

    	normalcell <- cells_copyscAT[grepl("sample1_",cells_copyscAT)]
	    keepIndx <- which(!colnames(cna_copyscAT) %in% normalcell)
	    cna_copyscAT <- cna_copyscAT[,c(keepIndx)]

	    copyscAT_bed <- data.frame(chr=cna_copyscAT[,1],start=as.numeric(cna_copyscAT[,2]),end=as.numeric(cna_copyscAT[,3]),binID=rownames(cna_copyscAT))
	    cna_copyscAT_mt <- cna_copyscAT[,4:ncol(cna_copyscAT)]

    	###Evaluating rare CNVs metrics
      res_copy <- c()
      Psets  <- seq(0.1,1,0.1)
		  for(prop in Psets){
				cna_copyscAT_bulk_align <- get_copyscAT_bulk(cna_copyscAT,bed.ref,cnvCellProp.min=prop)
		       
				bin_dat_truth <- GT_100k[[j]]
        bin_dat_truth$truth <- 0
				bin_dat_truth$truth[rareCNVRegion_index] <- 1

				bin_dat_truth$pred <- ifelse(cna_copyscAT_bulk_align$copyscAT.CN!=2 & !is.na(cna_copyscAT_bulk_align$copyscAT.CN), 1, 0)
		    bin_dat_truth$pred[!rareCNVRegion_index] <- 0
		    evares_ls <- evaluate_Binary(bin_dat_truth[["truth"]],bin_dat_truth[["pred"]])

        evares3 <- data.frame(method_prop=prop,
        			clone_truth=j,
           	  accuracy=as.numeric(evares_ls$accuracy),
              precision=as.numeric(evares_ls$precision),
              recall=as.numeric(evares_ls$recall),
              F1=as.numeric(evares_ls$f1))
        
        res_copy <- rbind(res_copy,evares3)
      }

      res_copy$nth <- nth
	  	res_copy$group <- group
	  	res_copy$CNVFrac <- pct
	  	res_copy$method <- "copyscAT"

	  	 write.table(res_copy,file = res_rare_all_file_copy,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_rare_all_file_copy),
		    append = TRUE
		  )


	  	##3.1 Counting the proportion of cells with clonal CNVs
	  	# 		# cna_copyscAT_bulk = cna_copyscAT_bulk %>% 
		  #     #       mutate(copyscAT.CN = case_when(
		  #     #           copyscAT_ratio > 2.5 ~ 'amp',
		  #     #           copyscAT_ratio < 1.5 ~ 'del',
		  #     #           T ~ 'neu'
		  #     #       ))
		  #     #
	  	obs_cells <- colnames(cna_copyscAT)[4:ncol(cna_copyscAT)]
	  	truthClone <- sapply(strsplit(obs_cells,"_"),"[",2)
	  	truthClone <- gsub("C","",truthClone)
	  	cellinfo<- data.frame(cellname=obs_cells,truthClone=truthClone)
	  	N_truth <- table(truthClone)
	  	prop_truth <- N_truth/sum(N_truth)


			res<- do.call(rbind,lapply(1:length(CNVs_clonal),function(x,CNVs_clonal,N_truth,prop_truth,copyscAT_bed){
				clonalCNVs <- CNVs_clonal[x]
				cat(paste0(x,":",clonalCNVs,"\n"))
				cnv_index <- rownames(CNVs_clonal_bed) %in% clonalCNVs
				clonalCNVs_bed <- CNVs_clonal_bed[cnv_index,,drop=FALSE]

				cnv_value <- clonalCNVs_bed[,4:(3+j)]
				cloneTr <- names(clonalCNVs_bed[,4:(3+j)])[which(clonalCNVs_bed[,4:(3+j)]!=2)]
				prop_truth_cloneTr <- sum(prop_truth[cloneTr])
				clonalCNVs_bed$clone_truth <- paste(cloneTr,collapse=",")

				intersect_res <- bt.intersect(a = clonalCNVs_bed, b = copyscAT_bed, wa = TRUE, wb = TRUE)	
				if(length(intersect_res)>0){
					colnames(intersect_res)<- c(colnames(clonalCNVs_bed),c("chr_pred","start_pred","end_pred","binID_pred"))
					row_index <- which(rownames(cna_copyscAT_mt)%in%intersect_res$binID_pred)

					Npred <- rowSums(cna_copyscAT_mt>2.5|cna_copyscAT_mt<1.5,na.rm=TRUE)
					Npred <- mean(Npred[row_index])
					prop_pred <- Npred/ncol(cna_copyscAT_mt)
					cnvFreq_results <- data.frame(
			      	Ncell_pred = Npred,
			      	prop_pred = prop_pred,
			      	clone_true = clonalCNVs_bed[["clone_truth"]],
			      	Ncell_truth = sum(N_truth[cloneTr]),
			      	prop_truth = prop_truth_cloneTr,
			      	binID = clonalCNVs)
				}else{
						cnvFreq_results <- c()
				}
			 return(cnvFreq_results)
			 },CNVs_clonal,N_truth,prop_truth,copyscAT_bed))
			res$prop_dist <- abs(as.numeric(res$prop_pred)-as.numeric(res$prop_truth))
			res$method <- "copyscAT"
			res$nth <- nth
	  	res$group <- group
	  	res$CNVFrac <- pct

	  	 write.table(res,file = res_clonalCNV_freq_file_copy,
		    sep = ",",
		    row.names = FALSE,
		    col.names = !file.exists(res_clonalCNV_freq_file_copy),
		    append = TRUE
		  )


}}}











###------------------------###
####Visualization: Part 1
###------------------------###
library(ggthemes)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggpubr)
library(rstatix)
library(dplyr)
outdir <- "./TeaCNV/simulation/diffCNVEva";
setwd(outdir)
outdir_data <- paste0(outdir,"/data");

col_method <- c(TeaCNV='#c85e62',epiAneuFinder='#67a583',copyscAT='#7b95c6')
darken_color <- function(color, factor=0.7){
    col <- col2rgb(color)
    col <- col * factor
    col <- rgb(t(col), maxColorValue=255)
    return(col)
}

res_bACC_all <- read.csv(paste0(outdir_data,"/evaRes_TeaCNV.csv")) 
res_bACC_all <- res_bACC_all %>%
	dplyr::select(-clone_pred,-clone_truth)%>%
  dplyr::group_by(CNVFrac,group,nth)%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE), .names = "{.col}"))%>%
    as.data.frame()


res_bACC_all_CopyscAT <- read.csv(paste0(outdir_data,"/evaRes_CopyscAT.csv"))
res_bACC_all_CopyscAT[is.na(res_bACC_all_CopyscAT)] <- 0
res_bACC_all_CopyscAT <- res_bACC_all_CopyscAT %>%
	dplyr::select(-method_prop,-clone_truth)%>%
  dplyr::group_by(CNVFrac,group,nth)%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE), .names = "{.col}"))%>%
    as.data.frame()


res_bACC_all_epiAneu <- read.csv(paste0(outdir_data,"/evaRes_epiAneu.csv"))
res_bACC_all_epiAneu[is.na(res_bACC_all_epiAneu)] <- 0
res_bACC_all_epiAneu <- res_bACC_all_epiAneu %>%
	dplyr::select(-method_prop,-clone_truth)%>%
  dplyr::group_by(CNVFrac,group,nth)%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE), .names = "{.col}"))%>%
    as.data.frame()



df_tea <- res_bACC_all%>%
	dplyr::select(accuracy,precision,recall,F1,RMSE.CNVvalue,bp_dist, nth,group,CNVFrac)%>%
	mutate(method="TeaCNV")

df_epi <- res_bACC_all_epiAneu%>%
	dplyr::select(accuracy,precision,recall,F1,RMSE.CNVvalue,bp_dist, nth,group,CNVFrac)%>%
	mutate(method="epiAneuFinder")

df_copy <- res_bACC_all_CopyscAT%>%
	dplyr::select(accuracy,precision,recall,F1,RMSE.CNVvalue,bp_dist, nth,group,CNVFrac)%>%
	mutate(method="copyscAT")
df_all <- rbind(df_tea,df_epi,df_copy)
colnames(df_all) <- c('accuracy','precision','recall','F1','RMSE','mean_BP_dist_Mb', 'nth','group','CNVFrac','method')
df_all[['mean_BP_dist_Mb']] <- round(df_all[['mean_BP_dist_Mb']],2)
df_all <- df_all %>%dplyr::filter(CNVFrac<0.8)


color_box <- "#cccccc"   
color_line <- "#0072B2" 
border_colors <- sapply(col_method, darken_color)

metric_cols1 <- c('accuracy','precision','recall','F1')
long_df <- df_all %>%
 #dplyr::filter(CNVFrac<=0.8) %>%
  dplyr::select(all_of(c("method", "CNVFrac","group", metric_cols1))) %>%
  reshape2::melt(id.vars = c("method", "CNVFrac","group"), variable.name = "Metric", value.name = "Value")%>%
  as.data.frame()
summary_df <- long_df %>%
  group_by(method, Metric) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
summary_df
write.csv(summary_df,paste0(outdir,"/eva_diffCNV_3method_F1score.summary.csv"),row.names=FALSE)

summary_2 <- long_df %>% dplyr::filter(method=="TeaCNV" )%>%
  group_by(CNVFrac,Metric) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
write.csv(summary_2,paste0(outdir,"/eva_diffCNV_TeaCNV_F1score_summary.csv"),row.names=FALSE)

###Global qualitative assessment of CNV event identification
my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneuFinder", "TeaCNV"))

plots <- list()
plots2 <- list()
for (metric_name in metric_cols1) {
  plot_df <- long_df %>% dplyr::filter(Metric == metric_name)
  p <- ggplot(plot_df, aes(x = method, y = Value,fill=method)) +
          geom_boxplot(aes(color = method),alpha = 0.7, outlier.shape = NA, width = 0.5) +
           # geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.3), 
           #      show.legend = T,size = 1)+
          stat_compare_means(comparisons = my_comparisons,method = "t.test",method.args = list(alternative = "less"), aes(label=..p.adj..)) +
          labs(title = "", x = "", y = metric_name) +
            theme_classic(base_size = 13) +
            theme(
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black",angle = 45, hjust = 1,vjust=1),
              axis.title = element_text(face = "bold"),
              axis.line = element_line(color = "black")
            )+
    		scale_fill_manual(values =  col_method)+
    		scale_color_manual(values=border_colors)+
      coord_cartesian(ylim = c(0,1), clip = "off")+
     theme(plot.margin = margin(10, 5, 2, 5))

         plots[[metric_name]] <- p


  plot_df2 <- long_df %>% dplyr::filter(Metric == metric_name &method=="TeaCNV" )
  summary_plot_df <- plot_df2 %>%
	  group_by(CNVFrac) %>%
	  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")

	p2 <- ggplot(plot_df2, aes(x = factor(.data[["CNVFrac"]]), y = Value,fill=factor(CNVFrac))) +
	  geom_boxplot( alpha = 0.6, outlier.shape = NA, width = 0.5) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = "", x = "CNV fraction", y = metric_name) +
      theme_classic(base_size = 13) +
      theme(
      	 panel.grid.major = element_line(color = "grey80", size = 0.5), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )+
      scale_fill_manual(values =  c('#f0f0f0', '#cccccc', '#a8a8a8'))+
      coord_cartesian(ylim = c(0,1), clip = "off")+
     theme(plot.margin = margin(5, 5, 5, 5))
  plots2[[metric_name]] <- p2

}

final_plot <- wrap_plots(plots, ncol = 4)
ggsave(paste0(outdir,"/eva_diffCNV_3method_F1score.pdf"),final_plot,width=15,height=4.5)
final_plot2 <- wrap_plots(plots2, ncol = 2)
ggsave(paste0(outdir,"/eva_diffCNV_TeaCNV_F1score.pdf"),final_plot2,width=10,height=6)






###Group comparison RMSE;  Break point distance (Mb)

df_all <- df_all%>%
	group_by(method) %>%
	mutate(BP_dist_variance = var(mean_BP_dist_Mb))

metric_cols2 <- c('RMSE','mean_BP_dist_Mb','BP_dist_variance') 

long_df2 <- df_all %>%
  dplyr::select(all_of(c("method", "CNVFrac","group", metric_cols2))) %>%
  reshape2::melt(id.vars = c("method", "CNVFrac","group"), variable.name = "Metric", value.name = "Value")%>%
  as.data.frame()

summary_df <- long_df2 %>%
  group_by(method, Metric) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
  	max_value = max(Value,na.rm=TRUE), .groups = "drop")
summary_df
write.csv(summary_df,paste0(outdir,"/eva_diffCNV_3method_RMSE_BPdistance.summary.csv"),row.names=FALSE)


my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneuFinder", "TeaCNV"))

plots1 <- list()
plots2 <- list()
for (metric_name in metric_cols2) {
	plot_df <- long_df2 %>% dplyr::filter(Metric == metric_name)

	stat_test <- plot_df %>%
		  # group_by(CNVFrac) %>%
		  rstatix::t_test(Value ~ method,comparisons = my_comparisons) %>%
		    add_significance("p.adj") %>%
  			add_xy_position(x = "method", group = "method",fun="max", dodge = 0.8, step.increase = 0.1)
	stat_test$p_label <- paste0("p=",formatC(stat_test$p.adj, format = "e", digits = 1))

	ylim1 <- c(0,max(stat_test$y.position)*1.3)

	 p <- ggplot(plot_df, aes(x = factor(method), y = Value)) +
          geom_boxplot(aes(fill = method),alpha = 0.7, outlier.shape = NA, width = 0.5,position = position_dodge(0.8)) +
				  stat_pvalue_manual(
				    stat_test, 
				    label = "p_label",  
				    tip.length = 0.01,
				    bracket.size = 0.6,
				    step.increase = 0.1,    
				    hide.ns = TRUE         
				  ) +
          labs(title = "", x = "", y = metric_name) +
            theme_minimal(base_size = 13) +
            theme(
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black",angle = 45, hjust = 1,vjust=1),
              axis.title = element_text(face = "bold"),
              axis.line = element_line(color = "black")
            )+
    		scale_fill_manual(values =  col_method)+
      coord_cartesian(ylim = ylim1)

  # if(metric_name=="RMSE"){
  # 	median_df <- plot_df %>%
	# 	  group_by(CNVFrac) %>%
	# 	  summarise(median_value = median(Value), .groups = "drop")
	#   cor_df <- median_df %>%
	# 	  summarise(
	# 	    cor = cor(CNVFrac, median_value),
	# 	    r_squared = cor^2,
	# 	    label = sprintf("r=%.2f, R²=%.2f", cor, r_squared)
	# 	  )
	#   annotation_text <- paste(cor_df$label, collapse = "\n")
	#   x_levels <- levels(factor(median_df$CNVFrac))
	# 	x_positions <- seq_along(x_levels)
	# 	x_mapping <- data.frame(
	# 	  CNVFrac = median_df$CNVFrac,
	# 	  x_position = x_positions
	# 	)
	# 	median_df <- median_df %>% left_join(x_mapping, by = "CNVFrac")
	# 	regression_lines <- median_df %>%
	# 	  do({
	# 	    model <- lm(median_value ~ x_position, data = .)
	# 	    pred <- predict(model, interval = "confidence")
	# 	    data.frame(., pred)
	# 	  })

	#   p <- p+
	#     geom_line(data = regression_lines, 
  #           aes(x=x_position,y = fit),size = 0.5,color="#0072B2") +
	#      annotate("text", x = 0.5, y = Inf, 
  #          label = annotation_text, 
  #          hjust = 0, vjust = 1, size = 4.5)+
	#      scale_x_discrete(labels = x_levels)

  # }
   plots1[[metric_name]] <- p
  

   plot_df2 <- long_df2 %>% dplyr::filter(Metric == metric_name &method=="TeaCNV" )
  summary_plot_df <- plot_df2 %>%
	  group_by(CNVFrac) %>%
	  summarise(mean_value = mean(Value, na.rm = TRUE),
	  	max_value = max(Value,na.rm=TRUE), .groups = "drop")
	 # summary_plot_df 

	ylim <- c(0,max(summary_plot_df$mean_value)*2)

	p2 <- ggplot(plot_df2, aes(x = factor(.data[["CNVFrac"]]), y = Value,fill=factor(CNVFrac))) +
	  geom_boxplot( alpha = 0.6, outlier.shape = NA, width = 0.5) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = "", x = "CNV fraction", y = metric_name) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )+
      coord_cartesian(ylim = ylim)+
      scale_fill_manual(values =  c('#f0f0f0', '#cccccc', '#a8a8a8'))
  plots2[[metric_name]] <- p2

   
}

final_plot1 <- wrap_plots(plots1, ncol = 2)
ggsave(paste0(outdir,"/eva_diffCNV_3method_RMSE.pdf"),final_plot1,width=8,height=4)

final_plot2 <- wrap_plots(plots2, ncol = 1)
ggsave(paste0(outdir,"/eva_diffCNV_TeaCNV_RMSE.pdf"),final_plot2,width=4,height=5)



###CNV identify precision&recall 
# pldata2 <-df_all %>%
#     dplyr::group_by(method,CNVFrac,group)%>%
#     summarise(across(everything(), mean, .names = "{.col}"))%>%
#     as.data.frame()
P1 = ggplot(df_all, aes(x=recall, y=precision, color = method,size=CNVFrac)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 0.5) + 
    geom_point()+
    xlim(0,1)+
    ylim(0,1)+
    theme_few()+
    scale_color_manual(values =  col_method)+
    scale_size(range = c(0.8, 3))

ggsave(filename = paste0(outdir,"/precision.recall_CNVidentify.pdf"),P1,width =5,height = 3.5)





###Statistics of running time and memory
cnvdir_tea <- "./TeaCNV/simulation/TeaCNVres_diffCNVpct"
CopyscAT_path <- "./TeaCNV/simulation/CopyscAT_diffCNVpct/"
epiAneu_path <- "./TeaCNV/simulation/epiAneufinder_diffCNVpct/"

cnvProp_max <- c(0.3,0.4,0.6,0.8,0.9,1)
clonal_groups=c("Monoclonal","Biclonal","Triclonal","Tetraclonal")
info_total <- c()
for(i in 1:length(cnvProp_max)){
	pct <- cnvProp_max[i]
	for(j in 1:4){
		group <- clonal_groups[j]
		for(nth in 1:10){
			cat("\n")
			print(paste0("Starting ",group,"-CNVpct",pct,"-",nth,"..."))
			cat("\n")
			indir_tea <- paste0(cnvdir_tea,"/CNVpct",pct,"/",group,"/",nth)
			indir_copy <- paste0(CopyscAT_path,"/CNVpct",pct,"/",group,"/",nth)
			indir_epiAneu <- paste0(epiAneu_path,"/CNVpct",pct,"/",group,"/",nth)

			if(file.exists(file.path(indir_tea,"runtime_log.txt"))){
				log_t <- extract_runtime_info(file.path(indir_tea,"runtime_log.txt"))
				log_t$method <- "TeaCNV"
			}
			if(file.exists(file.path(indir_copy,"runtime_log.txt"))){
				log_c <- extract_runtime_info(file.path(indir_copy,"runtime_log.txt"))
				log_c$method <- "copyscAT"
			}
			if(file.exists(file.path(indir_epiAneu,"runtime_log.txt"))){
				log_e <- extract_runtime_info(file.path(indir_epiAneu,"runtime_log.txt"))
				log_e$method <- "epiAneuFinder"
			}
			log_all <- rbind(log_t,log_c,log_e)
		  log_all$nth <- nth
	  	log_all$group <- group
	  	log_all$CNVFrac <- pct
	  	info_total <- rbind(info_total,log_all)


		}
	}
}
info_total$Memory_MB[info_total$Memory_MB<0] <- NA
write.csv(info_total,paste0(outdir,"/runtime_log_total.csv"),row.names=FALSE)

##plot
info_total <- read.csv(paste0(outdir,"/runtime_log_total.csv"))
metric_cols <- c('Memory_MB','Runtime_Minutes')
long_df <- info_total %>%
	dplyr::filter(CNVFrac<0.8)%>%
  dplyr::select(all_of(c("method", "CNVFrac","group", metric_cols))) %>%
  reshape2::melt(id.vars = c("method", "CNVFrac","group"), variable.name = "Metric", value.name = "Value")%>%
  as.data.frame()

summary_df <- long_df %>%
  group_by(method,Metric) %>%
  summarise(median_value = median(Value, na.rm = TRUE), .groups = "drop")
write.csv(summary_df,paste0(outdir,"/eva_diffCNV_3method_runtime.summary.csv"),row.names=FALSE)
summary_df


my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneuFinder", "TeaCNV"))
col_method <- c(TeaCNV='#c85e62',epiAneuFinder='#67a583',copyscAT='#7b95c6')
darken_color <- function(color, factor=0.7){
    col <- col2rgb(color)
    col <- col * factor
    col <- rgb(t(col), maxColorValue=255)
    return(col)
}
border_colors <- sapply(col_method, darken_color)

plots1 <- list()
plots2 <- list()
for (metric_name in metric_cols) {
  plot_df <- long_df %>% dplyr::filter(Metric == metric_name)
  stat_test <- plot_df %>%
		  # group_by(CNVFrac) %>%
		  rstatix::t_test(Value ~ method,comparisons = my_comparisons) %>%
		    add_significance("p.adj") %>%
  			add_xy_position(x = "method", group = "method", dodge = 0.8, step.increase = 0.05)
  stat_test$p_label <- paste0("p=",formatC(stat_test$p.adj, format = "e", digits = 1))

  p <- ggplot(plot_df, aes(x = factor(method), y = Value)) +
          geom_boxplot(aes(fill = method, color = method),alpha = 0.7, outlier.shape = NA, width = 0.5,position = position_dodge(0.8)) +
				  stat_pvalue_manual(
				    stat_test, 
				    label = "p_label", 
				    tip.length = 0.01,
				    bracket.size = 0.6,
				    step.increase = 0.02,     
				    hide.ns = TRUE        
				  ) +
          labs(title = "", x = "", y = metric_name) +
            theme_minimal(base_size = 13) +
            theme(
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black",angle = 45, hjust = 1,vjust=1),
              axis.title = element_text(face = "bold"),
              axis.line = element_line(color = "black")
            )+
    		scale_fill_manual(values =  col_method)+
    		scale_color_manual(values=border_colors)
  plots1[[metric_name]] <- p

  plot_df2 <- long_df %>% dplyr::filter(Metric == metric_name &method=="TeaCNV" )
  summary_plot_df <- plot_df2 %>%
	  group_by(CNVFrac) %>%
	  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")

	ylim <- c(0,max(summary_plot_df$mean_value)*2)
	
	p2 <- ggplot(plot_df2, aes(x = factor(.data[["CNVFrac"]]), y = Value,fill=factor(CNVFrac))) +
	  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = "", x = "CNV fraction", y = metric_name) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )+
      coord_cartesian(ylim = ylim)+
      scale_fill_manual(values =  c('#f0f0f0', '#cccccc', '#a8a8a8'))
  plots2[[metric_name]] <- p2
}

final_plot <- wrap_plots(plots1, ncol = 2)
ggsave(paste0(outdir,"/eva_diffCNV_3method_runtime.pdf"),final_plot,width=8.5,height=4)

final_plot2 <- wrap_plots(plots2, ncol = 1)
ggsave(paste0(outdir,"/eva_diffCNV_TeaCNV_runtime.pdf"),final_plot2,width=4,height=5)







