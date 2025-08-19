##Rare clone data
###Simulate scATAC based on real Tumor clone
suppressMessages({
	library(simATAC)
	library(SingleCellExperiment)
	library(GenomicRanges)
	library(BSgenome.Hsapiens.UCSC.hg38)
	library(SummarizedExperiment)
	library(dplyr)
	library(data.table)
	library(plyranges)
	library(Matrix)

})

wkdir <- "./TeaCNV"
outdir <- paste0(wkdir,"/simulation");if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)
outdir_simData <- paste0(outdir,"/simData");if(!file.exists(outdir_simData)){dir.create(outdir_simData,recursive=T)}

#If no real scATAC-seq data is available in the format of bin by cell matrix，we still can generate synthetic samples using the default parameters that are provided in the simATAC package. The default values are the parameters estimated from the GSE99172 data （K562 cells: 人类白血病细胞系).
##load the bin by cell matrix as seed contains cells having similar biological characteristics. The sparse matrix can be created by Snaptools or or a customozied pipleine to generate bin-by-cell matrix from real data.

###------------------------------------###
## Step 1: generate or load bin-cell matrix
###------------------------------------###
sampleID <- "ccRCC4"
datadir <- "./TeaCNV/Rdata"
inputdata <- paste0(datadir,"/Kidney_scMutiOmic_2Sample.rds")
obj <- readRDS(inputdata)
counts = obj[["peaks"]]@counts
ref_group_names <- c("Myeloid","Tcell","Bcell")
cellAnno = obj@meta.data

cnvpath <- "./TeaCNV"
cells_select <- read.csv(paste0(cnvpath, "/cells_filtered_in.csv"))$x
cellAnno <- cellAnno[cells_select,,drop=F]
cellAnno <- cellAnno[cellAnno$sampleID %in%sampleID,,drop=F]
cells_ref <- rownames(cellAnno)[cellAnno$group %in% ref_group_names]
counts_ref <- counts[,cells_ref]

cnvresFile <- paste0(cnvpath,"/",sampleID,"/final.CNVres.rds")
cnvres <- readRDS(cnvresFile)
clone_info <- cnvres$cellinfo
cells_clone1 <- clone_info$cellname[clone_info$clone=="1"]
counts_clone1 <- counts[,cells_clone1]

cells_clone2 <- clone_info$cellname[clone_info$clone=="2"]
counts_clone2 <- counts[,cells_clone2]

cells_clone3 <- clone_info$cellname[clone_info$clone=="3"]
counts_clone3 <- counts[,cells_clone3]


cells_clone4 <- clone_info$cellname[clone_info$clone=="4"]
counts_clone4 <- counts[,cells_clone4]

rm(obj)

cnvresFile2 <- paste0(cnvpath,"/ccRCC3/final.CNVres.rds")
cnvres2 <- readRDS(cnvresFile2)
clone_info2 <- cnvres2$cellinfo
cells_clone5 <- clone_info2$cellname[clone_info2$clone=="2"]
counts_clone5 <- counts[,cells_clone5]




###------------------------------------###
## Step 2: Simulation
###------------------------------------###
seeds <- sample(1:10000, 10)

#2.1 simulate reference data
outdir_ref <- paste0(outdir_simData,"/refCells");if(!file.exists(outdir_ref)){dir.create(outdir_ref,recursive=T)}
object <- simATACEstimate(counts_ref)
# params <- getParameters(object, c("nBins", "species", "nCells"))
# params
###simulate 10 times
for(nth in 1:10){
	object <- setParameters(object, seed = seeds[nth])
	sim_ref <- simATACSimulate(object,default = FALSE,nCells = 2000)
	count_sim <- counts(sim_ref)
	rownames(count_sim) <- rownames(counts_ref)
	# head(rowData(sim_ref))
	# head(colData(sim_ref))
	# simATACCompare(assay(sim_ref), counts_ref, outdir_ref, 'RefCells')
	saveRDS(count_sim,paste0(outdir_ref,"/counts_ref",nth,".rds"))
}
rm(object)



###------------------------------------###
## Step 3: simulate dataset according to different clonal structure
###------------------------------------###
###3.1 Simulate two clones (90% primary + 10% rare), N=50, 100, 200, 300, 400, 500, ... 1000, gradient 100, evaluate rare clone identification, repeat 10 times
outdir3 <- paste0(outdir_simData,"/Biclonal");if(!file.exists(outdir3)){dir.create(outdir3,recursive=T)}
NcRatio <- c(0.9,0.1)

sizeSet <- c(50,seq(100,1000,by=100),2000,5000,10000)
object1 <- simATACEstimate(counts_clone1)
object2 <- simATACEstimate(counts_clone2)
seeds <- sample(1:10000, 10)

for(size in sizeSet){
	outdir_size <- paste0(outdir3,"/size",size);if(!file.exists(outdir_size)){dir.create(outdir_size,recursive=T)}
	for(nth in 1:10){
		object1 <- setParameters(object1, seed = seeds[nth])
		#sample from clone 1 (Dominant)
		sim1 <- simATACSimulate(object1,default = FALSE,nCells = size*NcRatio[1])
		count_sim1 <- counts(sim1)
		rownames(count_sim1) <- rownames(counts_clone1)
		colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))


		#sample from clone 2
		object2 <- setParameters(object2, seed = seeds[nth]+1000)
		sim2 <- simATACSimulate(object2,default = FALSE,nCells = size*NcRatio[2])
		count_sim2 <- counts(sim2)
		rownames(count_sim2) <- rownames(counts_clone2)
		colnames(count_sim2) <- paste0("C2_",colnames(count_sim2))


		counts_total <- cbind(count_sim1,count_sim2)
		saveRDS(counts_total,paste0(outdir_size,"/Biclonal_counts_size",size,"_",nth,".rds"))
	}
}


###3.2 Simulate three clones (rare 10%), N = 50, 100, 200, 300, 400, 500, ... 1000, with a gradient of 100, evaluate rare clone identification, repeat 10 times
outdir4 <- paste0(outdir_simData,"/Triclonal");if(!file.exists(outdir4)){dir.create(outdir4,recursive=T)}
NcRatio <- c(0.5,0.4,0.1)
sizeSet <- c(50,seq(100,1000,by=100),2000,5000,10000)
object1 <- simATACEstimate(counts_clone1)
object2 <- simATACEstimate(counts_clone2)
object3 <- simATACEstimate(counts_clone3)
seeds <- sample(1:10000, 10)

for(size in sizeSet){
	outdir_size <- paste0(outdir4,"/size",size);if(!file.exists(outdir_size)){dir.create(outdir_size,recursive=T)}
	for(nth in 1:10){
		#sample from clone 1 (Dominant)
		object1 <- setParameters(object1, seed = seeds[nth])
		sim1 <- simATACSimulate(object1,default = FALSE,nCells = size*NcRatio[1])
		count_sim1 <- counts(sim1)
		rownames(count_sim1) <- rownames(counts_clone1)
		colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))


		#sample from clone 2
		object2 <- setParameters(object2, seed = seeds[nth]+1000)
		sim2 <- simATACSimulate(object2,default = FALSE,nCells = size*NcRatio[2])
		count_sim2 <- counts(sim2)
		rownames(count_sim2) <- rownames(counts_clone2)
		colnames(count_sim2) <- paste0("C2_",colnames(count_sim2))

		#sample from clone 3
		object3 <- setParameters(object3, seed = seeds[nth]+2000)
		sim3 <- simATACSimulate(object3,default = FALSE,nCells = size*NcRatio[3])
		count_sim3 <- counts(sim3)
		rownames(count_sim3) <- rownames(counts_clone3)
		colnames(count_sim3) <- paste0("C3_",colnames(count_sim3))

		counts_total <- cbind(count_sim1,count_sim2,count_sim3)
		saveRDS(counts_total,paste0(outdir_size,"/Triclonal_counts_size",size,"_",nth,".rds"))
	}
}


###3.3 Simulate four clones (rare 10%), N = 50, 100, 200, 300, 400, 500, ... 1000, with a gradient of 100, evaluate rare clone identification, repeat 10 times
outdir5 <- paste0(outdir_simData,"/Tetraclonal");if(!file.exists(outdir5)){dir.create(outdir5,recursive=T)}
NcRatio <- c(0.4,0.3,0.2,0.1)
sizeSet <- c(50,seq(100,1000,by=100),2000,5000,10000)
object1 <- simATACEstimate(counts_clone1)
object2 <- simATACEstimate(counts_clone2)
object3 <- simATACEstimate(counts_clone3)
object4 <- simATACEstimate(counts_clone4)
seeds <- sample(1:10000, 10)

for(size in sizeSet){
	outdir_size <- paste0(outdir5,"/size",size);if(!file.exists(outdir_size)){dir.create(outdir_size,recursive=T)}
	for(nth in 1:10){
		#sample from clone 1 (Dominant)
		object1 <- setParameters(object1, seed = seeds[nth])
		sim1 <- simATACSimulate(object1,default = FALSE,nCells = size*NcRatio[1])
		count_sim1 <- counts(sim1)
		rownames(count_sim1) <- rownames(counts_clone1)
		colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))


		#sample from clone 2
		object2 <- setParameters(object2, seed = seeds[nth]+1000)
		sim2 <- simATACSimulate(object2,default = FALSE,nCells = size*NcRatio[2])
		count_sim2 <- counts(sim2)
		rownames(count_sim2) <- rownames(counts_clone2)
		colnames(count_sim2) <- paste0("C2_",colnames(count_sim2))

		#sample from clone 3
		object3 <- setParameters(object3, seed = seeds[nth]+2000)
		sim3 <- simATACSimulate(object3,default = FALSE,nCells = size*NcRatio[3])
		count_sim3 <- counts(sim3)
		rownames(count_sim3) <- rownames(counts_clone3)
		colnames(count_sim3) <- paste0("C3_",colnames(count_sim3))

		#sample from clone 4
		object4 <- setParameters(object4, seed = seeds[nth]+3000)
		sim4 <- simATACSimulate(object4,default = FALSE,nCells = size*NcRatio[4])
		count_sim4 <- counts(sim4)
		rownames(count_sim4) <- rownames(counts_clone4)
		colnames(count_sim4) <- paste0("C4_",colnames(count_sim4))

		counts_total <- cbind(count_sim1,count_sim2,count_sim3,count_sim4)
		saveRDS(counts_total,paste0(outdir_size,"/Tetraclonal_counts_size",size,"_",nth,".rds"))
	}
}


###3.4 Simulate five clones (rare 10%), N = 50, 100, 200, 300, 400, 500, ... 1000, with a gradient of 100, evaluate rare clone identification, repeat 10 times
outdir6 <- paste0(outdir_simData,"/Pentaclonal");if(!file.exists(outdir6)){dir.create(outdir6,recursive=T)}
NcRatio <- c(0.3,0.2,0.2,0.2,0.1)
sizeSet <- c(50,seq(100,1000,by=100),2000,5000,10000)
object1 <- simATACEstimate(counts_clone1)
object2 <- simATACEstimate(counts_clone2)
object3 <- simATACEstimate(counts_clone3)
object4 <- simATACEstimate(counts_clone4)
object5 <- simATACEstimate(counts_clone5)
seeds <- sample(1:10000, 10)
for(size in sizeSet){
	outdir_size <- paste0(outdir6,"/size",size);if(!file.exists(outdir_size)){dir.create(outdir_size,recursive=T)}
	for(nth in 1:10){
		#sample from clone 1 (Dominant)
		object1 <- setParameters(object1, seed = seeds[nth])
		sim1 <- simATACSimulate(object1,default = FALSE,nCells = size*NcRatio[1])
		count_sim1 <- counts(sim1)
		rownames(count_sim1) <- rownames(counts_clone1)
		colnames(count_sim1) <- paste0("C1_",colnames(count_sim1))


		#sample from clone 2
		object2 <- setParameters(object2, seed = seeds[nth]+1000)
		sim2 <- simATACSimulate(object2,default = FALSE,nCells = size*NcRatio[2])
		count_sim2 <- counts(sim2)
		rownames(count_sim2) <- rownames(counts_clone2)
		colnames(count_sim2) <- paste0("C2_",colnames(count_sim2))

		#sample from clone 3
		object3 <- setParameters(object3, seed = seeds[nth]+2000)
		sim3 <- simATACSimulate(object3,default = FALSE,nCells = size*NcRatio[3])
		count_sim3 <- counts(sim3)
		rownames(count_sim3) <- rownames(counts_clone3)
		colnames(count_sim3) <- paste0("C3_",colnames(count_sim3))

		#sample from clone 4
		object4 <- setParameters(object4, seed = seeds[nth]+3000)
		sim4 <- simATACSimulate(object4,default = FALSE,nCells = size*NcRatio[4])
		count_sim4 <- counts(sim4)
		rownames(count_sim4) <- rownames(counts_clone4)
		colnames(count_sim4) <- paste0("C4_",colnames(count_sim4))

		#sample from clone 5
		object5 <- setParameters(object5, seed = seeds[nth]+4000)
		sim5 <- simATACSimulate(object5,default = FALSE,nCells = size*NcRatio[5])
		count_sim5 <- counts(sim5)
		rownames(count_sim5) <- rownames(counts_clone5)
		colnames(count_sim5) <- paste0("C5_",colnames(count_sim5))

		counts_total <- cbind(count_sim1,count_sim2,count_sim3,count_sim4,count_sim5)
		saveRDS(counts_total,paste0(outdir_size,"/Pentaclonal_counts_size",size,"_",nth,".rds"))
	}
}

