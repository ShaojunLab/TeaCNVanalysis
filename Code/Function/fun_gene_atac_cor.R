#' @param gene.coords GRanges object containing coordinates of genes in the
#' expression assay. If NULL, extract from gene annotations stored in the assay.
#' @param n_sample Number of peaks to sample at random when computing the null distribution.

#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom BSgenome.Hsapiens.UCSC.hg38
#' @importMethodsFrom Matrix t
# source("fun_Grange.R")
suppressPackageStartupMessages({
	library(future)
	library(future.apply)
	library(pbapply)
	#library(Matrix)
  library(BSgenome.Hsapiens.UCSC.hg38)
  require("BiocParallel")
})

gene_atac_cor <- function(object,peak.assay,gene.mat,peak.mat,
  links=NULL,
	gene.coords = NULL,
	method = "pearson",
	n_sample = 500,
	pvalue_cutoff = 0.05,
	score_cutoff = 0.05,
  core=6,
	verbose=TRUE){
	if (method == "pearson") {
	  cor_method <- qlcMatrix::corSparse
	} else if (method == "spearman") {
	  cor_method <- SparseSpearmanCor
	} else {
	  stop("method can be one of 'pearson' or 'spearman'.")
	}
	features.match <- c("GC.percent", "count", "sequence.length")

	all.peaks <- rownames(x = peak.mat)
	peak.data <- t(peak.mat)
	genes.use <- rownames(gene.mat)


	if (is.null(x = gene.coords)) {
    annot <- Annotation(object = object[[peak.assay]])
    if (is.null(x = annot)) {
    	annot <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    	
      #stop("Gene annotations not found")
    }
    gene.coords <- CollapseToLongestTranscript(
      ranges = annot
    )
    gene.coords.use <- gene.coords[match(genes.use,gene.coords$gene_name),]
	}
	if(is.null(links)){
	  links <- find_local_peaks(object,genes.use,annotations=Annotation(object = object[[peak.assay]]))
	}
	
  meta.features <- GetAssayData(
    object = object, assay = peak.assay, slot = "meta.features"
  )
  if(!all(grepl(paste(features.match,collapse = "|"),colnames(meta.features)))){
    object <- RegionStats(object = object, 
                       genome = BSgenome.Hsapiens.UCSC.hg38,
                       assay=peak.assay)
    meta.features <- GetAssayData(
      object = object, assay = peak.assay, slot = "meta.features"
    )
  }
  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
  # if (nbrOfWorkers() > 1) {
  #   mylapply <- future_lapply
  # } else {
  #   mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  # }
  mylapply <- bplapply
# run in parallel across genes
res <- mylapply(
X = seq_along(along.with = genes.use),
FUN = function(i) {
  gene.expression <- t(x = gene.mat[genes.use[[i]], , drop = FALSE])
  gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
  peak.use <- colnames(peak.data)
  peak.use <- links$peak[links$gene_name==genes.use[[i]]]
  
  if (length(peak.use) < 2) {
    # no peaks close to gene
    return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
  } else {
    peak.access <- peak.data[, peak.use, drop = FALSE]
    coef.result <- cor_method(
      X = peak.access,
      Y = gene.expression
    )
    rownames(x = coef.result) <- colnames(x = peak.access)
    coef.result <- coef.result[!is.nan(coef.result), , drop = FALSE]
    coef.result <- coef.result[abs(x = coef.result) > score_cutoff, , drop = FALSE]

    if (nrow(x = coef.result) == 0) {
      return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
    } else {

      # select peaks at random with matching GC content and accessibility
      # sample from peaks on a different chromosome to the gene
      peaks.test <- rownames(x = coef.result)
      trans.peaks <- all.peaks[
        !grepl(pattern = paste0("^", gene.chrom), x = all.peaks)
      ]
      meta.use <- meta.features[trans.peaks, ]
      pk.use <- meta.features[peaks.test, ]
      bg.peaks <- lapply(
        X = seq_len(length.out = nrow(x = pk.use)),
        FUN = function(x) {
          MatchRegionStats(
            meta.feature = meta.use,
            query.feature = pk.use[x, , drop = FALSE],
            features.match = features.match,
            n = n_sample,
            verbose = FALSE
          )
        }
                )
      # run background correlations
      bg.access <- peak.data[, unlist(x = bg.peaks), drop = FALSE]
      bg.coef <- cor_method(
        X = bg.access,
        Y = gene.expression
      )
      rownames(bg.coef) <- colnames(bg.access)
      zscores <- vector(mode = "numeric", length = length(x = peaks.test))
      for (j in seq_along(along.with = peaks.test)) {
        coef.use <- bg.coef[(((j - 1) * n_sample) + 1):(j * n_sample), ]
        coef.use <- coef.use[!is.nan(coef.use)]
        z <- (coef.result[j] - mean(x = coef.use)) / sd(x = coef.use)
        zscores[[j]] <- z
      }
      names(x = coef.result) <- peaks.test
      names(x = zscores) <- peaks.test
      zscore.vec <- c(zscore.vec, zscores)
      gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
      coef.vec <- c(coef.vec, coef.result)
    }
    gc(verbose = FALSE)
    pval.vec <- pnorm(q = -abs(x = zscore.vec))
    links.keep <- pval.vec < pvalue_cutoff
    if (sum(x = links.keep) == 0) {
      return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
    } else {
      gene.vec <- gene.vec[links.keep]
      coef.vec <- coef.vec[links.keep]
      zscore.vec <- zscore.vec[links.keep]
      return(list("gene" = gene.vec, "coef" = coef.vec, "zscore" = zscore.vec))
    }
  }
}, 
BPPARAM = MulticoreParam(workers = core)
)

  # combine results
  gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 1))
  coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 2))
  zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 3))
  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(object)
  }
  peak.key <- seq_along(
    along.with = unique(x = names(x = coef.vec))
  )
  names(x = peak.key) <- unique(x = names(x = coef.vec))
  coef.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(x = coef.vec)],
    x = coef.vec,
    dims = c(length(x = genes.use), max(peak.key))
  )
  rownames(x = coef.matrix) <- genes.use
  colnames(x = coef.matrix) <- names(x = peak.key)
  links.new <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  # add zscores
  z.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(x = zscore.vec)],
    x = zscore.vec,
    dims = c(length(x = genes.use), max(peak.key))
  )
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use)
  links.new$zscore <- z.lnk$score
  links.new$pvalue <- pnorm(q = -abs(x = links.new$zscore))
  links.new <- links.new[links.new$pvalue < pvalue_cutoff]
  return(links.new)
}


#subcnv <- read.table("pseudoCNV_200cells_links.txt",header=T,sep="\t")
#colnames(subcnv) <- gsub("\\.","-",colnames(subcnv))
#gene.mat<- read.table("pseudoExp_200cells_links.txt",header=T,sep="\t")
#colnames(gene.mat) <- gsub("\\.","-",colnames(gene.mat))
#peak.mat<- read.table("pseudoPeak_200cells_links.txt",header=T,sep="\t")
#colnames(peak.mat) <- gsub("\\.","-",colnames(peak.mat))

gene_atac_cor.my <- function(gene.mat,peak.mat,links.df=NULL,pvalue_cutoff = 0.05,
                             score_cutoff = 0.05){
  colnames(links.df)[grepl("seqnames|chr",colnames(links.df),ignore.case=TRUE)] <- "chr"
  pearsonR=WGCNA::cor(t(log(gene.mat+1,base=2)),t(log(peak.mat+1,base = 2)),use = "pairwise.complete.obs") ###This will take a long time...
  index1=match(links.df$gene_name,row.names(pearsonR))
  index2=match(links.df$peak,colnames(pearsonR))
  RR<-unlist(lapply(1:dim(links.df)[1], function(i,pearsonR,index1,index2){
    return(pearsonR[index1[i],index2[i]])
  },pearsonR,index1,index2))
  links.df$pearsonR <- RR
  ##correlation null distribution for trans regulation
  ##for each gene, randomly select 10000 peaks on other chromosome and estimated the mean and standard
  links.df$chr= do.call(rbind,strsplit(links.df$gr_gene_region,split="-"))[,1]
  chromosome=unique(links.df$chr)
  chromo_null <- lapply(chromosome, function(x,links.df,Exp,Peak){
    uniqueGene=links.df$gene_name[links.df$chr %in% x]
    uniquepeak=links.df$peak[links.df$chr!=x]
    randompeak=sample(uniquepeak,size=1000)
    gene.mat=Exp[uniqueGene,]
    peak.mat=Peak[randompeak,]
    exp_peak.cor <- WGCNA::cor(t(log(gene.mat+1,base=2)),t(log(peak.mat+1,base = 2)),use = "pairwise.complete.obs")
    return(exp_peak.cor)
  },links.df,gene.mat,peak.mat)
  
  ####estimate empirical P value
  pvalue <- lapply(1:length(chromosome),function(i,chromosome,chromo_null,links.df){
    x=links.df$pearsonR[links.df$chr==chromosome[i]]
    pvalue=pnorm(x,mean=mean(chromo_null[[i]][!is.na(chromo_null[[i]])]),sd=sd(chromo_null[[i]][!is.na(chromo_null[[i]])]),lower.tail=FALSE)
    return(pvalue)
  },chromosome,chromo_null,links.df)
  links.df$pvalue=NA
  for (i in 1:length(chromosome)){
    links.df$pvalue[links.df$chr==chromosome[i]]=pvalue[[i]]
  }
  #write.table(links.df,"totalLinksCor_200kb_VarGene.txt",col.names = T,row.names = F,quote=F,sep="\t")
  #links.df <- read.table("totalLinksCor_200kb_VarGene.txt",header=T,sep="\t")
  
  ###(4)  Significant links
  sigLink=links.df[links.df$pvalue< pvalue_cutoff&!is.na(links.df$pvalue)&abs(links.df$pearsonR)>score_cutoff,]
  sigLink=sigLink[order(sigLink$pearsonR,decreasing = T),]
  return(sigLink)
}




