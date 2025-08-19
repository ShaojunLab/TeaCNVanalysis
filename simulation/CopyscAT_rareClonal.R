#CopyscAT for simulated fragment

suppressMessages({
  library(CopyscAT)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(compiler)
})
script_path <- "./CopyscAT/"
py_script=paste0(script_path,"/process_fragment_file.py")
dir.ref <- paste0(script_path,"/hg38_references")

# ==== Script Start: Record start time and memory usage ====
start_time <- Sys.time()
mem_before <- get_mem_usage()

if(!file.exists(paste0(dir.ref,"/hg38_chrom_sizes.tsv"))){
  generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38",tileWidth = 1e6,outputDir = dir.ref)
}
##### REGULAR WORKFLOW #####
#initialize the environment (note: this needs to be done with every session in which you run Copy-scAT)
initialiseEnvironment(genomeFile=paste0(dir.ref,"/hg38_chrom_sizes.tsv"),
                      cytobandFile=paste0(dir.ref,"/hg38_1e+06_cytoband_densities_granges.tsv"),
                      cpgFile=paste0(dir.ref,"/hg38_1e+06_cpg_densities.tsv"),
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

###----------------###
# START...
###----------------###

#SET OUTPUT DEFAULT DIRECTORY AND NAME
wkdir = './TeaCNV/simulation'
simDatadir <- "./TeaCNV/simulation/simData/Fragment"
outdir <- paste0(wkdir,'/CopyscAT_rareClone');ifelse(!dir.exists(file.path(outdir)), dir.create(file.path(outdir)), FALSE)

clonal_groups=c("Biclonal","Triclonal","Tetraclonal", "Pentaclonal")
tumor_size <- c(50,seq(100,1000,by=100),2000,5000,10000)


for(g in 1:length(clonal_groups)){
  group <- clonal_groups[g]
  print(group)
  simDatadir2 <- paste0(simDatadir,"/",group)
  for(size in tumor_size){
    for(nth in 1:10){
      fragPath = paste0(simDatadir2,"/size",size,"/",nth,"/merged_fragments.tsv.gz")

      outdir_clt = paste0(outdir,"/",group,"/size",size,"/",nth);
      ifelse(!dir.exists(file.path(outdir_clt)), dir.create(file.path(outdir_clt)), FALSE)
      setwd(outdir_clt)

      setOutputFile(outdir_clt,paste(group,size,nth,sep="_"))
      outputfile <- paste0(outdir_clt,"/fragments_bin.tsv")
      dir.ref2 <- dir.ref
      command <- paste("python", py_script,"-i",fragPath,"-o",outputfile,"-b 1000000","-f",5000,"-g",paste0(dir.ref2,"/hg38_chrom_sizes.tsv"))
      system(command, intern = TRUE)

      scData <- readInputTable(outputfile)
      rownames(scData) <- paste0(rownames(scData),"-1")

      #3.
      scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=200,dividingFactor=1,upperFilterQuantile = 0.98)
      summaryFunction<-cutAverage
      scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 2,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 

      scData_collapse<-filterCells(scData_collapse,minimumSegments = 0,minDensity = 0)
      graphCNVDistribution(scData_collapse,outputSuffix = "_violin") 
      median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)
      #PART 2: ASSESSMENT OF CHROMOSOME-LEVEL CNVs 
      #cleanup step
      subsetSize <- ifelse(length(scData_collapse)>=200,200,length(scData_collapse))
      ###
      candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=subsetSize,fakeCellSD = 0.08, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff= 0.2,medianQuantileCutoff = 0.4)
      candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5) #= 1.5)
      saveRDS(candidate_cnvs_clean,paste0(outdir_clt,"/candidate_cnvs_clean.rds"))
      #OPTION 2: automatically identify non-neoplastic cells and use these for control
      nmf_results<-identifyNonNeoplastic(scData_collapse,methodHclust="ward.D",cutHeight = 0.4)
      write.table(x=rownames_to_column(data.frame(nmf_results$cellAssigns),var="Barcode"),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_clusters.csv"),quote=FALSE,row.names = FALSE,sep=",")

      #' @param minAlteredCells filtering parameter - remove all alterations with a smallest group less than X cells (default is 40)默认filterResults=TRUE
      ###
      final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "",sdCNV = 0.6,filterResults=TRUE,filterRange=0.4)
      ###
      ncells_CNV = candidate_cnvs_clean[[1]] %>% dplyr::filter(str_detect(Cells,"X",negate=TRUE)) %>% gather(chrom,cluster,starts_with("chr")) %>% group_by(chrom,cluster) %>% dplyr::filter(cluster!=0) %>% summarise(num=n()) %>% group_by(chrom) %>% summarise(min=min(num)) 
      ###removed chrom
      # %>% dplyr::filter(min<40)


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

      writeLines(log_text, con = "runtime_log.txt")
      cat(log_text)
    }
  }
}


    










 
