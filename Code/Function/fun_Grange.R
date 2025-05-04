# Find peaks near genes (function from Signac)
#
# Find peaks that are within a given distance threshold to each gene
#
# @param peaks A GRanges object containing peak coordinates
# @param genes A GRanges object containing gene coordinates
# @param distance Distance threshold. Peaks within this distance from the gene
# will be recorded.
# @param sep Separator for peak names when creating results matrix
#
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix sparseMatrix
#' @importFrom GenomicRanges resize
#
# @return Returns a sparse matrix
suppressPackageStartupMessages({
  library(IRanges)
  library(GenomicRanges)
  library(S4Vectors)
  #library(Matrix)
  library(tidyr)
})

find_local_peaks <- function(object,features,
                             peak.assay="ATAC",
                             peak_percentile_min = 0.1,
                             distance=5e+05,
                             toTSS=FALSE,
                             annotations=NULL,
                             gene.id = FALSE,
                             CollapseTranscript=TRUE){

  peakVar <- object[[peak.assay]]@meta.features
  peakVar = peakVar[order(peakVar$percentile,decreasing = T),]
  featurePeak = peakVar[peakVar$percentile>peak_percentile_min,] 
  
  #### generate peak Grange
  gr_peak <- StringToGRanges(row.names(featurePeak))
  gr_peak$peak <- row.names(featurePeak)
  
  if (is.null(x = annotations)) {
    annotations <- Annotation(object = object[[peak.assay]])
    #annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    if (is.null(x = annotations)) {
      stop("Gene annotations not found")
    }
  }
  if(CollapseTranscript){
    gr_gene <- CollapseToLongestTranscript(ranges = annotations)
    }else{
      gr_gene <-annotations
    }
  
  if(gene.id){
    gr_gene.use <- gr_gene[gr_gene$gene_id %in% features,] ###Filtering by features
    gr_gene.use$gene_name <- gr_gene.use$gene_id
  }else{
    gr_gene.use <- gr_gene[gr_gene$gene_name %in% features,] ###Filtering by features
  }
  if (length(x = gr_gene.use) == 0) {
    stop("Could not find gene coordinates for requested genes")
  }
  if(!any(grepl("chr",seqlevels(gr_gene.use)))){
    if(any(grepl("chr",seqlevels(gr_peak)))){
      seqlevels(gr_gene.use) <- paste0("chr", seqlevels(gr_gene.use))
    }
  }
  
  
  if(toTSS){ ###TBD
    peak_distance_matrix <- DistanceToTSS(
      peaks = gr_peak,
      genes = gr_gene.use,
      distance = distance
    )
  }else{
    gr_gene_extend <- Extend(x = gr_gene.use,upstream = distance,downstream = distance)
    gr_gene_region <- paste(seqnames(gr_gene_extend),ranges(gr_gene_extend),sep="-")
    gr_gene_extend$gr_gene_region <- gr_gene_region
    ####find peaks in gene region
    ranges <- subsetByOverlaps(gr_peak, gr_gene_extend)
    hits <- findOverlaps(gr_peak, gr_gene_extend)
    idx_peak <- queryHits(hits)
    idx_gene <- subjectHits(hits)
    ranges_ov <- gr_peak[idx_peak,]
    mcols(ranges_ov) <- c(mcols(ranges_ov), mcols(gr_gene_extend[idx_gene,]))
    head(ranges_ov)
    ranges_ov <- as.data.frame(ranges_ov)
    if(gene.id){
      peaks_in <- ranges_ov[ranges_ov$gene_id %in% features,]
    }else{
      peaks_in <- ranges_ov[ranges_ov$gene_name %in% features,]
    }
    cat(paste0("Find ", nrow(unique(peaks_in[,c("peak","gene_name")]))," peak-gene pairs in ",distance,"kb gene region.\n"))
    cat(paste0("Including ", length(unique(peaks_in[,"peak"]))," peaks and ",length(unique(peaks_in[,"gene_name"]))," genes."))
    
  }
  return(peaks_in)
  
}



##' @param peaks peak file or GRanges object
##' @param genes 	A GenomicRanges object.
##' 
DistanceToTSS <- function(
    peaks,
    genes,
    distance = 200000,
    sep = c("-", "-")
) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
CollapseToLongestTranscript <- function(ranges) {
  library(data.table)
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

#tss.positions <- resize(gene.coords, width = 1, fix = 'start')

#' @param peaks string vector of genomic peaks, e.g. c("chrx-xxx-xxx","chrx-xxx-xxx")
#' @return GRanges object with gene ranges and metadata column of peaks("binID")
map_peaks2gene <- function(peaks,
                           genes=NULL,
                           annotations=NULL,
                           genome="hg38",
                           expand_distance = 0
                           ){
  suppressMessages({
    library(plyranges)
    library(ggpubr)
  })
  
  if (is.null(x = annotations)) {
    if(genome=="hg38"){
      require(EnsDb.Hsapiens.v86)
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    }else if (is.null(x = annotations)) {
      stop("Gene annotations not found")
    }
  }
  bins <- data.frame(peak=peaks)
  chromInfo <- separate(bins, peak, into = c("seqnames", "start","end"), sep = "-")
  bins <- data.frame(bins,chromInfo)
  peak_gr <- GRanges(seqnames = bins$seqnames,
                     ranges = IRanges(start = as.numeric(bins$start), end = as.numeric(bins$end)),
                     strand = "*")
  mcols(peak_gr) <- bins$peak
  
  if(is.null(genes)){
    genes <- unique(annotations$gene_name)
  }
  gene_loc <- annotate_gene(genes,genome=genome,annotations=annotations)
  gene_loc_expan <- Extend(x = gene_loc,upstream = expand_distance,downstream = expand_distance)
  gr_gene_region <- paste(seqnames(gene_loc_expan),ranges(gene_loc_expan),sep="-")
  gene_loc_expan$gene_region <- gr_gene_region
  
  hits <- findOverlaps(gene_loc_expan, peak_gr)
  idx_seg <- subjectHits(hits)
  idx_gene <- queryHits(hits)
  ranges_ov <- gene_loc_expan[idx_gene,]
  mcols(ranges_ov) <- c(mcols(ranges_ov), mcols(peak_gr[idx_seg,]))
  colnames(mcols(ranges_ov))[length(mcols(ranges_ov))] <- "binID"
  return(ranges_ov)
}

#' @param bed.query A bed file or a data frame of query regions
#' @param bed.subject A bed file or a data frame of subject regions
#' @return A data frame of the same order as the reference bed file, with the corresponding bin ID for each query region.
align_Grange2bin = function(bed.query,bed.subject,duplicate_by="binID"){
  suppressMessages({
    library(GenomicRanges)
  })
  if (Reduce("|", is(bed.query) == "character")) {
    if(substr(bed.query, nchar(bed.query)-3, nchar(bed.query)) == ".bed") {
      bed.query <- read.table(bed.query,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    }else if (Reduce("|", is(bed.query ) %in% c("data.frame"))) {
      bed.query  <- as.data.frame(bed.query )
    }
  }
  if (Reduce("|", is(bed.subject ) %in% c("data.frame"))) {
    bed.subject  <- as.data.frame(bed.subject )
  }
  ## check bin bed chr column
  if(!grepl("chr",bed.query[1,1])){
    bed.query[,1]=paste0("chr",bed.query[,1])
  }
  if(!grepl("chr",bed.subject[1,1])){
    bed.subject[,1]=paste0("chr",bed.subject[,1])
  }
  ## convert to GRanges
  query=GRanges(bed.query[,1], IRanges(as.numeric(bed.query[,2]),as.numeric(bed.query[,3])))
  subject=GRanges(bed.subject[,1], IRanges(as.numeric(bed.subject[,2]),as.numeric(bed.subject[,3])))
  ## align
  align=findOverlaps(query,subject)
  idx_query<- queryHits(align)
  idx_subject<- subjectHits(align)
  
  colnames(bed.subject)[1:3] <-c("chr.subject","start.subject","end.subject")
  col_subject_filt <- colnames(bed.subject)[!colnames(bed.subject)%in%colnames(bed.query)]

  bed.query2 <- cbind(bed.query[idx_query,],bed.subject[idx_subject,col_subject_filt])
  colnames(bed.query2)[1:3] = colnames(bed.query)[1:3] <-c("chr","start","end")
  bed.query2$binID <- paste(bed.query2$chr,bed.query2$start,bed.query2$end,sep="_")

  bed.query2 <- bed.query2[!duplicated(bed.query2[,duplicate_by]), ]
  
  bed.query$binID <- paste(bed.query[,1],bed.query[,2],bed.query[,3],sep="_")
  bed_final <- left_join(bed.query,bed.query2,by=intersect(colnames(bed.query),colnames(bed.query2)))
  
  return(as.data.frame(bed_final))
}

# bedtools merge
granges_merge <- function(df, 
                          min_gap = 0, 
                          group_cols = "chr", 
                          count_col = "merged_regions",
                          score_columns = NULL,
                          score_fun = max) {
  # 转换输入为GRanges对象
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  
  # 分组处理
  grouping <- as.list(df[group_cols])
  groups <- interaction(grouping)
  grl <- split(gr, groups)
  
  merged_list <- lapply(grl, function(g) {
    # 合并重叠/邻近区域
    reduced_gr <- reduce(g, min.gapwidth = min_gap, with.revmap = TRUE)
    
    # 初始化元数据列
    mcols(reduced_gr)[[count_col]] <- lengths(mcols(reduced_gr)$revmap)
    
    # 处理统计列
    if (!is.null(score_columns)) {
      # 获取原始区间与合并区间的映射关系
      ov <- findOverlaps(g, reduced_gr)
      
      # 遍历每个统计列
      for (col in score_columns) {
        # 提取原始数据
        original_values <- mcols(g)[[col]]
        
        # 按合并区间分组聚合
        agg_values <- tapply(
          original_values[queryHits(ov)], 
          subjectHits(ov),
          FUN = score_fun
        )
        
        # 添加统计结果列
        col_name <- paste(col, deparse(substitute(score_fun)), sep = "_")
        mcols(reduced_gr)[[col_name]] <- as.vector(agg_values)
      }
    }
    
    # 清理临时列
    mcols(reduced_gr)$revmap <- NULL
    return(reduced_gr)
  })
  
  # 合并所有分组结果
  merged_gr <- unlist(GRangesList(merged_list))
  
  # 转换为数据框
  merged_df <- as.data.frame(merged_gr, row.names = NULL) %>%
    dplyr::rename(chr = seqnames, start = start, end = end) %>%
    dplyr::select(-strand, -width)
  
  return(merged_df)
}

liftover_region <- function(region_str) {
  # chain_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  # download.file(chain_url, destfile = "hg19ToHg38.over.chain.gz")
  # R.utils::gunzip("hg19ToHg38.over.chain.gz")  # 需要安装R.utils包
  # chain <- import.chain("hg19ToHg38.over.chain")
  # 解析输入区域
  parts <- strsplit(region_str, "[:-]|[--]|[__]")[[1]]
  chr <- parts[1]
  start <- as.integer(parts[2])
  end <- as.integer(parts[3])
  
  # 创建GRanges对象
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  
  # 执行LiftOver转换
  converted <- liftOver(gr, chain)
  
  # 处理转换结果
  if (length(converted[[1]]) == 0) {
    return(data.frame(
      Original_Region = region_str,
      Merged_Region = NA_character_,
      chr=NA,
      start = NA,
      end=NA,
      Chromosome_Consistency = "Consistent",
      Mapping_Status = "Failed",
      stringsAsFactors = FALSE
    ))
  }
  # 提取转换后的所有区间
  converted_ranges <- converted[[1]]
   # 检查染色体一致性
  unique_chr <- unique(as.character(seqnames(converted_ranges)))
  if (length(unique_chr) > 1) {
    return(data.frame(
      Original_Region = region_str,
      Merged_Region = NA_character_,
      chr=NA,
      start = NA,
      end=NA,
      Chromosome_Consistency = "Inconsistent",
      Mapping_Status = "Failed",
      stringsAsFactors = FALSE
    ))
  }
   # 合并所有区间为最大范围（允许间隙存在）
  merged_start <- min(start(converted_ranges))
  merged_end <- max(end(converted_ranges))
  merged_coord <- paste0(unique_chr, ":", merged_start, "-", merged_end)
  
  data.frame(
    Original_Region = region_str,
    Merged_Region = merged_coord,
    chr=unique_chr,
    start = as.numeric(merged_start),
    end=as.numeric(merged_end),
    Chromosome_Consistency = "Consistent",
    Mapping_Status = "Success",
    stringsAsFactors = FALSE
  )
}


library(GenomicRanges)
library(dplyr)
library(rtracklayer)
#' Annotate genomic intervals with cytogenetic bands
#'
#' @param bed_df Input BED dataframe, must contain chr/start/end columns
#' @param cytoband_ref cytoband reference (CSV format 
#'        with chr/start/end/Cytoband columns)
#' @param genome_version Genome assembly version, supports hg19/hg38 (default)
#' @param keep_existing Whether to preserve existing Cytoband column values 
#'        (TRUE fills NAs only, FALSE overwrites)
#' @return Tibble with added Cytoband annotations
annotate_cytoband <- function(bed_df, 
                         cytoband_ref,
                         genome_version = "hg38",
                         keep_existing = TRUE){
  required_cols <- c("chr", "start", "end")
  if (!all(required_cols %in% colnames(bed_df))) {
    stop("Input data must contain columns: ", paste(required_cols, collapse = ", "))
  }

     cytoband_ref <- cytoband_ref %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  bed_gr <- bed_df %>%
    mutate(
      chr = factor(chr, levels = paste0("chr", c(1:22, "X", "Y"))),
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end"
    )

  overlaps <- findOverlaps(bed_gr, cytoband_ref)
  
  matched_cytobands <- tibble(
    query_idx = queryHits(overlaps),
    cytoband = mcols(cytoband_ref)$Cytoband[subjectHits(overlaps)]
  ) %>%
    group_by(query_idx) %>%
    summarise(
      new_cytoband = paste(unique(cytoband), collapse = ",")
    )

    bed_df <- bed_df %>%
    mutate(query_idx = row_number()) %>%
    left_join(matched_cytobands, by = "query_idx") %>%
    mutate(
      Cytoband = if_else(Cytoband == "NA", NA_character_, Cytoband),
      Cytoband = case_when(
        keep_existing & !is.na(Cytoband) ~ Cytoband,  # Preserve existing values
        nzchar(new_cytoband) ~ new_cytoband,          # Use new annotations
        TRUE ~ NA_character_                          # Maintain NA if no match
      )
    ) %>%
    dplyr::select(-query_idx, -new_cytoband)
  
  return(bed_df)


}




