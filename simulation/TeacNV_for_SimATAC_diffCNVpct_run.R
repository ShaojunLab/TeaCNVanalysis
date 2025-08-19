# Infering CNV using TeaCNV for simulated scATAC counts
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
  library(dplyr)
})



library(optparse)
args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("--pct"), type = "double", help = "Downsampling proportion (e.g., 0.3, 0.4,0.6)"),
  make_option(c("--ClonalType"), type = "character", help = "type of clonal structure （Biclonal,Triclonal,Tetraclonal,Pentaclonal)"),
  make_option(c("--nth"), type = "integer", help = "Iteration number (e.g., 1 to 10)"),
  make_option(c("--delt_lim"), type = "double",default = 0.4, help = "(default 0.4)"),
  make_option(c("--DEsegCoeff"), type = "double",default = 0.7, help = "(default 0.7)"),
  make_option(c("--CorrectByDist"), type = "logical", default = FALSE, help = "correct CNV of segment by dist [default %default]"),
  make_option(c("--resolution"), type = "double", default = 1, help = "Resolution for clustering"),
  make_option(c("--Zscore_cutoff"), type = "double", default = 1.28, help = "Zscore cutoff for clone merge"),
  make_option(c("--p.adj_CloneMerge"), type = "double", default = 0.05, help = "p.adj_CloneMerge for clustering")
 )
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

pct <- opt$pct
group <- opt$ClonalType
nth <- opt$nth

delt_lim <- opt$delt_lim
DEsegCoeff  <- opt$DEsegCoeff
correct_by_dist <- opt$CorrectByDist
seu_resolution <- opt$resolution
Zscore_cutoff <- opt$Zscore_cutoff
p.adj_CloneMerge <- opt$p.adj_CloneMerge

# pct <- 0.6
# group <- "Triclonal"
# nth <- 9
# delt_lim <- 0.25
# DEsegCoeff=0.6
# correct_by_dist=FALSE
# seu_resolution=1
# Zscore_cutoff = 1.28
# p.adj_CloneMerge=0.7

get_mem_usage <- function() {
  if (.Platform$OS.type == "windows") {
    memory.size()
  } else {
    sum(gc()[, 2])  # Linux/macOS: # Total memory used in MB
  }
}




script_path <- "./github/TeaCNV/"
# folder = "PSSD"

setwd(script_path)
genome = "hg38"
blacklistFile <- file.path(".", "data", "hg38-blacklist.v2.bed")

wkdir <- paste0("./TeaCNV/simulation/TeaCNVres_diffCNVpct")
if(!file.exists(wkdir)){dir.create(wkdir,recursive=T)}

simDatadir <- paste0("./TeaCNV/simulation/simData_diffCNVpct")

clonal_groups=c("Monoclonal","Biclonal","Triclonal","Tetraclonal")


simDatadir2 <- paste0(simDatadir,"/CNVpct",pct)


simDatadir3 <- paste0(simDatadir2,"/",group);

cat("\n\n")
print(paste0("Starting ",group,"-CNVpct",pct,"-",nth,"..."))
cat("\n\n")

refFile <- paste0(simDatadir,"/refCells/counts_ref",nth,".rds")
counts_ref <- readRDS(refFile)
colnames(counts_ref) <- paste0("ref_",colnames(counts_ref))

outdir_obs2 <- paste0(wkdir,"/CNVpct",pct,"/",group,"/",nth);if(!file.exists(outdir_obs2)){dir.create(outdir_obs2,recursive=T)}
setwd(outdir_obs2)
obsFile <- paste0(simDatadir3,"/counts_",nth,".rds")
counts_obs <- readRDS(obsFile)
colnames(counts_obs) <- paste0("obs_",colnames(counts_obs))
counts <- cbind(counts_ref,counts_obs)
cell_anno<- data.frame(row.names=colnames(counts),
                         Cluster=c(rep("normal", ncol(counts_ref)), rep("tumor", ncol(counts_obs)))
                         )
counts <- as.matrix(counts)
rownames(counts) <- gsub("_","-",rownames(counts))

SegSize_min= 100
# delt_lim = 0.4
scFactor = 1.2
# seu_resolution = 1
min_cells_in_group=20
diploidy_pct.max=0.95
# Zscore_cutoff = 3;
# p.adj_CloneMerge=0.1
if(group=="Monoclonal"){Zscore_cutoff=10;p.adj_CloneMerge=0}
if(pct==0.3){SegSize_min=60}

###TeaCNV starting...
# ==== Script Start: Record start time and memory usage ====
start_time <- Sys.time()
mem_before <- get_mem_usage()

###TeaCNV normalize count to count per kb 
cnv_obj <- CreateTeaCNVObject(input = counts,
                                   annotationFile = cell_anno,
                                   ref_group_names = "normal",
                                   ChrRemove = c('chrX', 'chrY', 'chrM'),
                                   FiltCell = TRUE,
                                   CellnCount_quant=c(0.05,0.95),
                                   cellproplim = 0.05,
                                   count_lim = 4,
                                   Correct_by_length=FALSE)

tryCatch(runTeaCNV(
    input_obj = cnv_obj,
    binSize = 1,
    outdir = outdir_obs2,
    seu_resolution = seu_resolution,
    SegLen_min = 2e6,
    SegSize_min = SegSize_min,
    seg_method = "PELT",
    StopStep = 4,
    delt_lim = delt_lim,
    scFactor = scFactor,
    Zscore_cutoff = Zscore_cutoff,
    min_cells_in_group=min_cells_in_group,
    diploidy_pct.max=diploidy_pct.max,
    p.adj_CloneMerge=p.adj_CloneMerge,
    DEsegCoeff=DEsegCoeff,
    correct_by_dist=correct_by_dist
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

writeLines(log_text, con = file.path(outdir_obs2, "runtime_log.txt"))
cat(log_text)


