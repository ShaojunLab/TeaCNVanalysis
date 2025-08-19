# Infering CNV using TeaCNV for simulated scATAC counts (Xupro & Imac)
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
  library(R.utils)
})

get_mem_usage <- function() {
  if (.Platform$OS.type == "windows") {
    memory.size()
  } else {
    sum(gc()[, 2])  # Linux/macOS: # Total memory used in MB
  }
}

script_path <- "./TeaCNV/"
setwd(script_path)
genome = "hg38"
blacklistFile <- file.path(".", "data", "hg38-blacklist.v2.bed")
wkdir <- "./TeaCNV/simulation/TeaCNV_res";
if(!file.exists(wkdir)){dir.create(wkdir,recursive=T)}
clonal_group <- c("Biclonal","Triclonal","Tetraclonal","Pentaclonal")
tumor_size <- c(50,seq(100,1000,by=100),2000,5000,10000)

for(i in 1:length(clonal_group)){
	# cg="Triclonal"
	# size=10000

	cg <- clonal_group[i]
	simDatadir <- paste0("./TeaCNV/simulation/simData/",cg)

	#length(tumor_size)
	for(j in 1:length(tumor_size)){
		size <- tumor_size[j]
		for(nth in 1:10){
			print(paste0("Starting ",cg,"-size",size,"-",nth,"..."))
			refFile <- paste0("./TeaCNV/simulation/simData/refCells/counts_ref",nth,".rds")

			counts_ref <- readRDS(refFile)
			colnames(counts_ref) <- paste0("ref_",colnames(counts_ref))
		
			outdir <- paste0(wkdir,"/",cg,"/size",size,"/",nth);if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
			setwd(outdir)

			obsFile <- paste0(simDatadir,"/size",size,"/",cg,"_counts_size",size,"_",nth,".rds")
			counts_obs <- readRDS(obsFile)
			counts <- cbind(counts_ref,counts_obs)
			cell_anno<- data.frame(row.names=colnames(counts),
		                           Cluster=c(rep("normal", 2000), rep("tumor", size))
		                           )
			counts <- as.matrix(counts)
			rm(counts_ref,counts_obs)
			invisible(gc())

			SegSize_min= 100
			delt_lim = 0.3
			scFactor = 1.2
			seu_resolution = 1
			if(size<=100){seu_resolution = 1.1}
			min_cells_in_group=4
			
			###TeaCNV starting...
			# ==== Script Start: Record start time and memory usage ====
			start_time <- Sys.time()
			mem_before <- get_mem_usage()

			cnv_obj_path <- file.path(outdir, "TeaCNV.obj")
			  cnv_obj <- tryCatch({
			    readRDS(cnv_obj_path)
			  }, error = function(e) {
			    CreateTeaCNVObject(input = counts,
			                       annotationFile = cell_anno,
			                       ref_group_names = "normal",
			                       ChrRemove = c('chrX', 'chrY', 'chrM'),
			                       FiltCell = FALSE,
			                       cellproplim = 0.05,
			                       count_lim = 4)
			  })


			tryCatch(runTeaCNV(
		      input_obj = cnv_obj,
		      binSize = 1,
		      outdir = outdir,
		      SegLen_min = 2e6,
		      SegSize_min = SegSize_min,
		      seg_method = "PELT",
		      seu_resolution=seu_resolution,
		      delt_lim = delt_lim,
		      min_cells_in_group=min_cells_in_group,
		      scFactor = scFactor,
		    ),error=function(e){message("TeaCNV failed: ", e$message)})

			# ==== Script End: Record end time and memory usage ====
			end_time <- Sys.time()
			mem_after <- get_mem_usage()

			run_time_min <- as.numeric(difftime(end_time, start_time, units = "mins"))
			mem_used <- mem_after - mem_before

			log_text <- paste0(
			  "=== R Script Runtime Log ===\n",
			  "Start time: ", start_time, "\n",
			  "End time: ", end_time, "\n",
			  "Total runtime (minutes): ", round(run_time_min, 3), "\n",
			  "Memory used (MB): ", round(mem_used, 3), "\n"
			)

			log_file <- paste0(outdir, "/runtime_log.txt")
			writeLines(log_text, con = log_file)
			cat(log_text)

		}



	}


}




####Bulk analysis
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
  library(circlize)
})

script_path="./github/TeaCNV"
cytoBandFile <- file.path(script_path, "data", "cytoBand_hg38.tsv")
wkdir <- "./TeaCNV/simulation/TeaCNV_res";
clonal_group <- c("Biclonal","Triclonal","Tetraclonal","Pentaclonal")
tumor_size <- c(50,seq(100,1000,by=100),2000,5000,10000)
delt_lim = 0.3

for(i in 1:length(clonal_group)){
	cg <- clonal_group[i]
	for(j in 1:length(tumor_size)){
		size <- tumor_size[j]
		for(nth in 1:10){
			print(paste0("Starting ",cg,"-size",size,"-",nth,"..."))
			outdir_clt <- paste0(wkdir,"/",cg,"/size",size,"/",nth)
			filepath = paste0(outdir_clt,"/TeaCNV.obj")
			if(file.exists(filepath)){
							input_obj <- readRDS(filepath)
				cell_anno_new <- input_obj@cell_anno.filt
				mtx_bin <- input_obj@data.binCount.norm

				segRatio_bulk <- CNratioInfer(input_matrix=mtx_bin,
	                                inputcellano=cell_anno_new,
	                                cytoBand=cytoBandFile,
	                                outdir=outdir_clt,
	                                group2bulkby="group",
	                                min_cells_in_group = 20,
	                                normallabel="reference",
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

	 			CNV_res_bulk <- CNV_esti_ini(outdir=outdir_clt,
	                               segScore_file=segRatio_bulk,
	                               filt_seg = TRUE,
	                               length_seg_cutoff =1e6,
	                               genome="hg38",
	                               true_loc = TRUE,
	                               outputFigure = F,
	                               outplot_name=paste0("CNV_bulk_initial"))

			  bulkCN_res <- ploidyRefine(CNV_res_bulk,delt.lim =delt_lim )
			  saveRDS(bulkCN_res,paste0(outdir_clt, "/bulkCNVres.rds"))
			  ratiores <- bulkCN_res$observed$input_BinSegRatio
			  integerCNVres <- bulkCN_res$observed$seg.dat
			  plt <- cellLineSeg_v3(binRatio.df=ratiores,integerCNV.df=integerCNVres,outdir=outdir_clt,
			                        #ylim = ylim,
			                        fname=paste0("Pseudobulk"),
			                        ggarrange_width=c(10,4),
			                        height = 2,width=10,outPlot=FALSE)
			 ggsave(paste0(outdir_clt, "/bulkCNVs.pdf"),plt$ggarranged_p, width=10, height=2,device = pdf,bg="white")
			}

		}
	}
}



























