#!/usr/bin/env Rscript


library(optparse)
args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("--pct"), type = "double", help = "Downsampling proportion (e.g., 0.3, 0.4,0.6)"),
  make_option(c("--ClonalType"), type = "character", help = "type of clonal structure （Biclonal,Triclonal,Tetraclonal,Pentaclonal)"),
  make_option(c("--nth"), type = "integer", help = "Iteration number (e.g., 1 to 10)")
 )

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


pct <- opt$pct
group <- opt$ClonalType
nth <- opt$nth



print("Running script for epiAneufinder")
get_mem_usage <- function() {
  if (.Platform$OS.type == "windows") {
    memory.size()
  } else {
    sum(gc()[, 2])  # Linux/macOS: # Total memory used in MB
  }
}

library(epiAneufinder)

#BSgenome to use for the analysis. The genome should be already installed in R. In this example we use the UCSC hg38 genome from Bioconductor
genome <- "BSgenome.Hsapiens.UCSC.hg38"
#Default is NULL
exclude <- c('chrX','chrY','chrM')
#Bed file with the blacklisted regions of the genome. This file is genome-version specific and it should be downloaded by the user
blacklist <- "./epiAnuefinder/hg38-blacklist.v2.bed" #Path and file name of the blacklisted regions in bed format. If you use hg38 genome then the blacklisted regions can be found in the sample_data folder.
#Window size for partitioning the genome. Smaller window sizes will result in longer running times. Default is 1e5
windowSize <- 1e5
#Parameter to instruct epiAneufinder to resume from a previous run. Can be set to either True or False
#If certain parameters change, for example minsize, resuming may end in error messages. In such a case change the parameter to False 
#Default in False
reuse.existing=TRUE
#Upper quantile thrshold. Default is 0.9
uq=0.9

#Lower quantile threshold. Default is at 0.1
lq=0.1


#Number of cores to use for the analysis. Default is 4
ncores=4

#Minimum number of fragments for a cell to be included in the analysis. This parameter is only for fragnment files. Default is 20000
minFrags = 5000

#Threshold for filtering bins if the ratio of cells with zero reads is higher than the threshold. Setting it to 0 deactivates the filter. Default is 0.85 
threshold_blacklist_bins=0.85

#Parameter on how many breakpoins to use for the CNV calculation. Default is 1, all breakpoints. If higher than one, the algorithm will calculate every n breakpoints
#Setting it to higher than 1 speeds the process with lower resolution as a result
minsize=1

#Number of segments per chromosomes (2^k). Default value is 3
k=4

#Number of consecutive bins to constitute a CNV
minsizeCNV=0


wkdir = './TeaCNV/simulation'
simDatadir <- "./TeaCNV/simulation/simData_diffCNVpct/Fragment"
outdir <- paste0(wkdir,'/epiAneufinder_diffCNVpct');ifelse(!dir.exists(file.path(outdir)), dir.create(file.path(outdir),recursive =TRUE), FALSE)

clonal_groups=c("Monoclonal","Biclonal","Triclonal","Tetraclonal")
cnvProp_max <- c(0.3,0.4,0.6,0.8,0.9,1)






# for(g in 1:length(clonal_groups)){
#   group <- clonal_groups[g]
#   print(group)
  # for(j in 1:5){
  #   pct <- cnvPct[j]
     simDatadir3 <- file.path(simDatadir, paste0("CNVpct", pct), group)
    # for(nth in 1:10){
			fragPath = paste0(simDatadir3,"/",nth,"/merged_fragments.tsv.gz")
			title_karyo=paste(group,pct,nth,sep="_")


			outdir_clt = paste0(outdir,"/CNVpct",pct,"/",group,"/",nth);
			ifelse(!dir.exists(file.path(outdir_clt)), dir.create(file.path(outdir_clt),recursive =TRUE), FALSE)
			setwd(outdir_clt)
			# ==== Script Start: Record start time and memory usage ====
			start_time <- Sys.time()
			mem_before <- get_mem_usage()


      epiAneufinder::epiAneufinder(input=fragPath, 
                             outdir=outdir_clt, 
                             blacklist=blacklist, 
                             windowSize=windowSize, 
                             genome=genome, 
                             exclude=exclude, 
                             reuse.existing=reuse.existing, 
                             uq=uq, lq=lq, 
                             title_karyo=title_karyo, 
                             ncores=ncores,
                             minFrags=minFrags,
                             minsize=minsize,
                             k=k,
                             threshold_blacklist_bins=threshold_blacklist_bins,
                             minsizeCNV=minsizeCNV)
			}

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

      log_file <- paste0(outdir_clt, "/runtime_log.txt")
      writeLines(log_text, con = log_file)
      cat(log_text)
  #   }
  # }







