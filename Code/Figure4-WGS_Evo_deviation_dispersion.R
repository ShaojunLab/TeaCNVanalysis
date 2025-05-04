library(ROCR)
library(pROC)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'

datapath = paste0(work_dir,'/EvaluationData')
if(!file.exists(datapath)){dir.create(datapath,recursive=T)}
setwd(datapath)

method <- c("epiAneu","copyscAT","TeaCNV")
sampleID <- c("ccRCC1","ccRCC3","ccRCC4")


###Figure 2g,2j
###deviation and dispersion
library(sfsmisc)
deviation <- c()
dispersion <- c()
ID="ccRCC1"
est <- read.csv(paste0("./",ID,"_CNdf_clean_clone.csv"))
est <- est[!est$chr %in% c("chrX","chrY"),]

est$groundtruth <- est$dna_ratio
est$groundtruth[est$dna_CNA=='amp'] <- 3
est$groundtruth[est$dna_CNA=="amp" & est$dna_ratio >1.5] <- 4
est$groundtruth[est$dna_CNA=='neu'] <- 2
est$groundtruth[est$dna_CNA=='del'] <- 1
est <- est[!is.na(est$groundtruth),]
base1 <- est$epiAneuRatio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$epiAneuRatio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$epiAneuRatio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$epiAneuRatio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
epiAneu <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$epiAneuRatio[!is.na(est$epiAneuRatio)]))
range <- max(est$epiAneuRatio[!is.na(est$epiAneuRatio)])-min(est$epiAneuRatio[!is.na(est$epiAneuRatio)])
epiAneu <- epiAneu*(max(est$groundtruth)-min(est$groundtruth))/range

lower=min(est$epiAneuRatio[!is.na(est$epiAneuRatio)])-1
upper=max(est$epiAneuRatio[!is.na(est$epiAneuRatio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
dc=density(base3,from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y)
dd$w12=pmin(dd$a,dd$b)
dd$w23=pmin(dd$b,dd$c)
dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
intersection=integrate.xy(dd$x,dd$w23)
overlap=overlap+2*intersection/total23
epiAneu.dis <- (2-overlap)/2

base1 <- est$copyscAT_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$copyscAT_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$copyscAT_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$copyscAT_ratio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
copyscAT <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)]))
range <- (max(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])-min(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)]))
copyscAT <- copyscAT*(max(est$groundtruth)-min(est$groundtruth))/range

lower=min(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])-1
upper=max(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
dc=density(base3,from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y)
dd$w12=pmin(dd$a,dd$b)
dd$w23=pmin(dd$b,dd$c)
dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
intersection=integrate.xy(dd$x,dd$w23)
overlap=overlap+2*intersection/total23
copyscAT.dis <- (2-overlap)/2

d <- est$TeaCNV_mean.1-est$groundtruth
d <- d[!is.na(d)]
TeaCNV <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.1[!is.na(est$TeaCNV_mean.1)])-min(est$TeaCNV_mean.1[!is.na(est$TeaCNV_mean.1)]))
lower=min(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)])-1
upper=max(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)])+1
da=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==3],from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y)
dd$w12=pmin(dd$a,dd$b)
dd$w23=pmin(dd$b,dd$c)
dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
intersection=integrate.xy(dd$x,dd$w23)
overlap=overlap+2*intersection/total23
TeaCNV.dis <- (2-overlap)/2
deviation <- rbind(deviation,data.frame(sample=ID,method=method,RMSE=c(epiAneu,copyscAT,TeaCNV)))
dispersion <- rbind(dispersion,data.frame(sample=ID,method=method,Dispersion=c(epiAneu.dis,copyscAT.dis,TeaCNV.dis)))



ID="ccRCC2"
est <- read.csv(paste0("./",ID,"_CNdf_clean_clone.csv"))
est <- est[!est$chr %in% c("chrX","chrY"),]

est$groundtruth <- est$dna_ratio
est$groundtruth[est$dna_CNA=='amp'] <- 3
est$groundtruth[est$dna_CNA=="amp" & est$dna_ratio >1.5] <- 4
est$groundtruth[est$dna_CNA=='neu'] <- 2
est$groundtruth[est$dna_CNA=='del'] <- 1
est <- est[!is.na(est$groundtruth),]
base1 <- est$epiAneuRatio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$epiAneuRatio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$epiAneuRatio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$epiAneuRatio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
epiAneu <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$epiAneuRatio[!is.na(est$epiAneuRatio)]))
range <- max(est$epiAneuRatio[!is.na(est$epiAneuRatio)])-min(est$epiAneuRatio[!is.na(est$epiAneuRatio)])
epiAneu <- epiAneu*(max(est$groundtruth)-min(est$groundtruth))/range
lower=min(est$epiAneuRatio[!is.na(est$epiAneuRatio)])-1
upper=max(est$epiAneuRatio[!is.na(est$epiAneuRatio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
dc=density(base3,from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y)
dd$w12=pmin(dd$a,dd$b)
dd$w23=pmin(dd$b,dd$c)
dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
intersection=integrate.xy(dd$x,dd$w23)
overlap=overlap+2*intersection/total23
epiAneu.dis <- (2-overlap)/2
base1 <- est$copyscAT_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$copyscAT_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$copyscAT_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$copyscAT_ratio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
copyscAT <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)]))
range <- (max(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])-min(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)]))
copyscAT <- copyscAT*(max(est$groundtruth)-min(est$groundtruth))/range
lower=min(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])-1
upper=max(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
#dc=density(base3,from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y)
dd$w12=pmin(dd$a,dd$b)
#dd$w23=pmin(dd$b,dd$c)
#dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
#total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
#intersection=integrate.xy(dd$x,dd$w23)
#overlap=overlap+2*intersection/total23
copyscAT.dis <- (1-overlap)

####inferCNV
base1 <- est$infercnv_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$infercnv_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$infercnv_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$infercnv_ratio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
inferCNV <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$infercnv_ratio[!is.na(est$infercnv_ratio)]))
range <- (max(est$infercnv_ratio[!is.na(est$infercnv_ratio)])-min(est$infercnv_ratio[!is.na(est$infercnv_ratio)]))
inferCNV <- inferCNV*(max(est$groundtruth)-min(est$groundtruth))/range
lower=min(est$infercnv_ratio[!is.na(est$infercnv_ratio)])-1
upper=max(est$infercnv_ratio[!is.na(est$infercnv_ratio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
#dc=density(base3,from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y)
dd$w12=pmin(dd$a,dd$b)
#dd$w23=pmin(dd$b,dd$c)
#dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
#total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
#intersection=integrate.xy(dd$x,dd$w23)
#overlap=overlap+2*intersection/total23
inferCNV.dis <- (1-overlap)


d <- est$TeaCNV_mean.1-est$groundtruth
d <- d[!is.na(d)]
TeaCNV1 <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.1[!is.na(est$TeaCNV_mean.1)])-min(est$TeaCNV_mean.1[!is.na(est$TeaCNV_mean.1)]))
lower=min(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)])-1
upper=max(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)])+1
da=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==3],from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y)
dd$w12=pmin(dd$a,dd$b)
dd$w23=pmin(dd$b,dd$c)
dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
intersection=integrate.xy(dd$x,dd$w23)
overlap=overlap+2*intersection/total23
TeaCNV.dis1 <- (2-overlap)/2

d <- est$TeaCNV_mean.2-est$groundtruth
d <- d[!is.na(d)]
TeaCNV2 <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.2[!is.na(est$TeaCNV_mean.2)])-min(est$TeaCNV_mean.2[!is.na(est$TeaCNV_mean.2)]))
lower=min(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)])-1
upper=max(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)])+1
da=density(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)&est$groundtruth==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)&est$groundtruth==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)&est$groundtruth==3],from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y)
dd$w12=pmin(dd$a,dd$b)
dd$w23=pmin(dd$b,dd$c)
dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
intersection=integrate.xy(dd$x,dd$w23)
overlap=overlap+2*intersection/total23
TeaCNV.dis2 <- (2-overlap)/2
d <- est$TeaCNV_mean.3-est$groundtruth
d <- d[!is.na(d)]
TeaCNV3 <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.3[!is.na(est$TeaCNV_mean.3)])-min(est$TeaCNV_mean.3[!is.na(est$TeaCNV_mean.3)]))
lower=min(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)])-1
upper=max(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)])+1
da=density(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)&est$groundtruth==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)&est$groundtruth==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)&est$groundtruth==3],from=lower,to=upper)
dd=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y)
dd$w12=pmin(dd$a,dd$b)
dd$w23=pmin(dd$b,dd$c)
dd$w13=pmin(dd$a,dd$c)
total12=integrate.xy(dd$x,dd$a)+integrate.xy(dd$x,dd$b)
intersection=integrate.xy(dd$x,dd$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd$x,dd$b)+integrate.xy(dd$x,dd$c)
intersection=integrate.xy(dd$x,dd$w23)
overlap=overlap+2*intersection/total23
TeaCNV.dis3 <- (2-overlap)/2

TeaCNV <- mean(c(TeaCNV1,TeaCNV2,TeaCNV3))
TeaCNV.dis <- mean(c(TeaCNV.dis1,TeaCNV.dis2,TeaCNV.dis3))
deviation <- rbind(deviation,data.frame(sample=ID,method=c(method,"inferCNV"),RMSE=c(epiAneu,copyscAT,TeaCNV,inferCNV)))
dispersion <- rbind(dispersion,data.frame(sample=ID,method=c(method,"inferCNV"),Dispersion=c(epiAneu.dis,copyscAT.dis,TeaCNV.dis,inferCNV.dis)))


ID="ccRCC3"
est <- read.csv(paste0("./",ID,"_CNdf_clean_clone.csv"))
est <- est[!est$chr %in% c("chrX","chrY"),]
est$groundtruth <- est$dna_ratio
est$groundtruth[est$dna_CNA=='amp'] <- 3
est$groundtruth[est$dna_CNA=="amp" & est$dna_ratio >1.5] <- 4
est$groundtruth[est$dna_CNA=='neu'] <- 2
est$groundtruth[est$dna_CNA=='del'] <- 1
est <- est[!is.na(est$groundtruth),]
base1 <- est$epiAneuRatio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$epiAneuRatio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$epiAneuRatio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$epiAneuRatio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
epiAneu <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$epiAneuRatio[!is.na(est$epiAneuRatio)]))
range <- max(est$epiAneuRatio[!is.na(est$epiAneuRatio)])-min(est$epiAneuRatio[!is.na(est$epiAneuRatio)])
epiAneu <- epiAneu*(max(est$groundtruth)-min(est$groundtruth))/range
lower=min(est$epiAneuRatio[!is.na(est$epiAneuRatio)])-1
upper=max(est$epiAneuRatio[!is.na(est$epiAneuRatio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
dc=density(base3,from=lower,to=upper)
dd=density(base4,from=lower,to=upper)
dd1=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y,d=dd$y)
dd1$w12=pmin(dd1$a,dd1$b)
dd1$w23=pmin(dd1$b,dd1$c)
dd1$w34=pmin(dd1$c,dd1$d)
total12=integrate.xy(dd1$x,dd1$a)+integrate.xy(dd1$x,dd1$b)
intersection=integrate.xy(dd1$x,dd1$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd1$x,dd1$b)+integrate.xy(dd1$x,dd1$c)
intersection=integrate.xy(dd1$x,dd1$w23)
overlap=overlap+2*intersection/total23
total34=integrate.xy(dd1$x,dd1$c)+integrate.xy(dd1$x,dd1$d)
intersection=integrate.xy(dd1$x,dd1$w34)
overlap=overlap+2*intersection/total34
epiAneu.dis <- (3-overlap)/3
base1 <- est$copyscAT_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$copyscAT_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$copyscAT_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$copyscAT_ratio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
copyscAT <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)]))
range <- (max(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])-min(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)]))
copyscAT <- copyscAT*(max(est$groundtruth)-min(est$groundtruth))/range
lower=min(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])-1
upper=max(est$copyscAT_ratio[!is.na(est$copyscAT_ratio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
dc=density(base3,from=lower,to=upper)
dd=density(base4,from=lower,to=upper)
dd1=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y,d=dd$y)
dd1$w12=pmin(dd1$a,dd1$b)
dd1$w23=pmin(dd1$b,dd1$c)
dd1$w34=pmin(dd1$c,dd1$d)
total12=integrate.xy(dd1$x,dd1$a)+integrate.xy(dd1$x,dd1$b)
intersection=integrate.xy(dd1$x,dd1$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd1$x,dd1$b)+integrate.xy(dd1$x,dd1$c)
intersection=integrate.xy(dd1$x,dd1$w23)
overlap=overlap+2*intersection/total23
total34=integrate.xy(dd1$x,dd1$c)+integrate.xy(dd1$x,dd1$d)
intersection=integrate.xy(dd1$x,dd1$w34)
overlap=overlap+2*intersection/total34
copyscAT.dis <- (3-overlap)/3
####inferCNV
base1 <- est$infercnv_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
d1 <- base1-mean(base1)
base2 <- est$infercnv_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
d2 <- base2-mean(base2)
base3 <- est$infercnv_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
d3 <- base3-mean(base3)
base4 <- est$infercnv_ratio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
d4 <- base4-mean(base4)
inferCNV <- sqrt((sum(d1^2)+sum(d2^2)+sum(d3^2)+sum(d4^2))/length(est$infercnv_ratio[!is.na(est$infercnv_ratio)]))
range <- (max(est$infercnv_ratio[!is.na(est$infercnv_ratio)])-min(est$infercnv_ratio[!is.na(est$infercnv_ratio)]))
inferCNV <- inferCNV*(max(est$groundtruth)-min(est$groundtruth))/range
lower=min(est$infercnv_ratio[!is.na(est$infercnv_ratio)])-1
upper=max(est$infercnv_ratio[!is.na(est$infercnv_ratio)])+1
da=density(base1,from=lower,to=upper)
db=density(base2,from=lower,to=upper)
dc=density(base3,from=lower,to=upper)
dd=density(base4,from=lower,to=upper)
dd1=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y,d=dd$y)
dd1$w12=pmin(dd1$a,dd1$b)
dd1$w23=pmin(dd1$b,dd1$c)
dd1$w34=pmin(dd1$c,dd1$d)
total12=integrate.xy(dd1$x,dd1$a)+integrate.xy(dd1$x,dd1$b)
intersection=integrate.xy(dd1$x,dd1$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd1$x,dd1$b)+integrate.xy(dd1$x,dd1$c)
intersection=integrate.xy(dd1$x,dd1$w23)
overlap=overlap+2*intersection/total23
inferCNV.dis <- (3-overlap)/3

d <- est$TeaCNV_mean.1-est$groundtruth
d <- d[!is.na(d)]
TeaCNV1 <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.1[!is.na(est$TeaCNV_mean.1)])-min(est$TeaCNV_mean.1[!is.na(est$TeaCNV_mean.1)]))
lower=min(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)])-1
upper=max(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)])+1
da=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==3],from=lower,to=upper)
dd=density(est$TeaCNV_RatioMean.1[!is.na(est$TeaCNV_mean.1)&est$groundtruth==4],from=lower,to=upper)
dd1=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y,d=dd$y)
dd1$w12=pmin(dd1$a,dd1$b)
dd1$w23=pmin(dd1$b,dd1$c)
dd1$w34=pmin(dd1$c,dd1$d)
total12=integrate.xy(dd1$x,dd1$a)+integrate.xy(dd1$x,dd1$b)
intersection=integrate.xy(dd1$x,dd1$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd1$x,dd1$b)+integrate.xy(dd1$x,dd1$c)
intersection=integrate.xy(dd1$x,dd1$w23)
overlap=overlap+2*intersection/total23
total34=integrate.xy(dd1$x,dd1$c)+integrate.xy(dd1$x,dd1$d)
intersection=integrate.xy(dd1$x,dd1$w34)
overlap=overlap+2*intersection/total34
TeaCNV.dis1 <- (3-overlap)/3
d <- est$TeaCNV_mean.2-est$groundtruth
d <- d[!is.na(d)]
TeaCNV2 <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.2[!is.na(est$TeaCNV_mean.2)])-min(est$TeaCNV_mean.2[!is.na(est$TeaCNV_mean.2)]))
lower=min(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)])-1
upper=max(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)])+1
da=density(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)&est$groundtruth==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)&est$groundtruth==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)&est$groundtruth==3],from=lower,to=upper)
dd=density(est$TeaCNV_RatioMean.2[!is.na(est$TeaCNV_mean.2)&est$groundtruth==4],from=lower,to=upper)
dd1=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y,d=dd$y)
dd1$w12=pmin(dd1$a,dd1$b)
dd1$w23=pmin(dd1$b,dd1$c)
dd1$w34=pmin(dd1$c,dd1$d)
total12=integrate.xy(dd1$x,dd1$a)+integrate.xy(dd1$x,dd1$b)
intersection=integrate.xy(dd1$x,dd1$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd1$x,dd1$b)+integrate.xy(dd1$x,dd1$c)
intersection=integrate.xy(dd1$x,dd1$w23)
overlap=overlap+2*intersection/total23
total34=integrate.xy(dd1$x,dd1$c)+integrate.xy(dd1$x,dd1$d)
intersection=integrate.xy(dd1$x,dd1$w34)
overlap=overlap+2*intersection/total34
TeaCNV.dis2 <- (3-overlap)/3
d <- est$TeaCNV_mean.3-est$groundtruth
d <- d[!is.na(d)]
TeaCNV3 <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.3[!is.na(est$TeaCNV_mean.3)])-min(est$TeaCNV_mean.3[!is.na(est$TeaCNV_mean.3)]))
lower=min(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)])-1
upper=max(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)])+1
da=density(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)&est$groundtruth==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)&est$groundtruth==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)&est$groundtruth==3],from=lower,to=upper)
dd=density(est$TeaCNV_RatioMean.3[!is.na(est$TeaCNV_mean.3)&est$groundtruth==4],from=lower,to=upper)
dd1=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y,d=dd$y)
dd1$w12=pmin(dd1$a,dd1$b)
dd1$w23=pmin(dd1$b,dd1$c)
dd1$w34=pmin(dd1$c,dd1$d)
total12=integrate.xy(dd1$x,dd1$a)+integrate.xy(dd1$x,dd1$b)
intersection=integrate.xy(dd1$x,dd1$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd1$x,dd1$b)+integrate.xy(dd1$x,dd1$c)
intersection=integrate.xy(dd1$x,dd1$w23)
overlap=overlap+2*intersection/total23
total34=integrate.xy(dd1$x,dd1$c)+integrate.xy(dd1$x,dd1$d)
intersection=integrate.xy(dd1$x,dd1$w34)
overlap=overlap+2*intersection/total34
TeaCNV.dis3 <- (3-overlap)/3
d <- est$TeaCNV_mean.4-est$groundtruth
d <- d[!is.na(d)]
TeaCNV4 <- sqrt(sum(d^2)/length(d))/(max(est$TeaCNV_mean.4[!is.na(est$TeaCNV_mean.4)])-min(est$TeaCNV_mean.4[!is.na(est$TeaCNV_mean.4)]))
lower=min(est$TeaCNV_RatioMean.4[!is.na(est$TeaCNV_mean.4)])-1
upper=max(est$TeaCNV_RatioMean.4[!is.na(est$TeaCNV_mean.4)])+1
da=density(est$TeaCNV_RatioMean.4[!is.na(est$TeaCNV_mean.4)&est$TeaCNV_mean.4==1],from=lower,to=upper)
db=density(est$TeaCNV_RatioMean.4[!is.na(est$TeaCNV_mean.4)&est$TeaCNV_mean.4==2],from=lower,to=upper)
dc=density(est$TeaCNV_RatioMean.4[!is.na(est$TeaCNV_mean.4)&est$TeaCNV_mean.4==3],from=lower,to=upper)
dd=density(est$TeaCNV_RatioMean.4[!is.na(est$TeaCNV_mean.4)&est$TeaCNV_mean.4==4],from=lower,to=upper)
dd1=data.frame(x=da$x,a=da$y,b=db$y,c=dc$y,d=dd$y)
dd1$w12=pmin(dd1$a,dd1$b)
dd1$w23=pmin(dd1$b,dd1$c)
dd1$w34=pmin(dd1$c,dd1$d)
total12=integrate.xy(dd1$x,dd1$a)+integrate.xy(dd1$x,dd1$b)
intersection=integrate.xy(dd1$x,dd1$w12)
overlap=2*intersection/total12
total23=integrate.xy(dd1$x,dd1$b)+integrate.xy(dd1$x,dd1$c)
intersection=integrate.xy(dd1$x,dd1$w23)
overlap=overlap+2*intersection/total23
total34=integrate.xy(dd1$x,dd1$c)+integrate.xy(dd1$x,dd1$d)
intersection=integrate.xy(dd1$x,dd1$w34)
overlap=overlap+2*intersection/total34
TeaCNV.dis4 <- (3-overlap)/3
TeaCNV <- mean(c(TeaCNV1,TeaCNV2,TeaCNV3,TeaCNV4))
TeaCNV.dis <- mean(c(TeaCNV.dis1,TeaCNV.dis2,TeaCNV.dis3,TeaCNV.dis4))
deviation <- rbind(deviation,data.frame(sample=ID,method=c(method,"inferCNV"),RMSE=c(epiAneu,copyscAT,TeaCNV,inferCNV)))
dispersion <- rbind(dispersion,data.frame(sample=ID,method=c(method,"inferCNV"),Dispersion=c(epiAneu.dis,copyscAT.dis,TeaCNV.dis,inferCNV.dis)))

##Fig 3f: RMSE
my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneu", "TeaCNV"))
deviation$metrics = 'RMSE'
P3 = ggplot(deviation[deviation$method!="inferCNV",], aes(x=method, y=RMSE)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3),size = 2, 
               aes(color=sample), show.legend = T)+
    #facet_wrap(~metrics, scale="free",ncol = 5)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    stat_compare_means(comparisons = my_comparisons,method = "t.test",method.args = list(alternative = "greater"))+
    theme_few()+ scale_color_nejm() + xlab('')

ggsave(plot = P3,filename = paste0('WGS.RMSE.pdf'),path = './',
       width = 4,height = 3)


dispersion$sample <- factor(dispersion$sample,levels = unique(dispersion$sample[order(dispersion$Dispersion,decreasing = T)]))
# p <- ggplot(dispersion, aes(x = sample, y = Dispersion, fill = method)) + 
#     geom_bar(stat="identity",position=position_dodge()) + 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# pdf(file=paste0("WGS.dispersion_",prop,".pdf"),width = 4,height = 4)
# p
# dev.off()

write.csv(deviation,file = paste0("WGS.deviation.csv"),quote = F,row.names = F)
write.csv(dispersion,file = paste0("WGS.dispersion.csv"),quote = F,row.names = F)

deviation %>% group_by(method) %>%summarise(RMSE_MEAN=mean(RMSE))


###consistent
CN.dispersion <- list()

ID="ccRCC3"
est <- read.csv(paste0("./",ID,"_CNdf_clean_clone.csv"))
est <- est[!est$chr %in% c("chrX","chrY"),]


est$groundtruth <- est$dna_ratio
est$groundtruth[est$dna_CNA=='amp'] <- 3
est$groundtruth[est$dna_CNA=="amp" & est$dna_ratio >1.5] <- 4
est$groundtruth[est$dna_CNA=='neu'] <- 2
est$groundtruth[est$dna_CNA=='del'] <- 1
est <- est[!is.na(est$groundtruth),]
base1 <- est$epiAneuRatio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$epiAneuRatio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$epiAneuRatio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])
epiAneuFinder.dispersion <- 1-overlap/length(c(base1,base3))

base1 <- est$copyscAT_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$copyscAT_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$copyscAT_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])
copyscAT.dispersion <- 1-overlap/length(c(base1,base3))

base1 <- est$infercnv_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$infercnv_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$infercnv_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])
inferCNV.dispersion <- 1-overlap/length(c(base1,base3))

base1 <- est$TeaCNV_mean.1[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.1[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.1[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion1 <- 1-overlap/length(c(base1,base3))
base1 <- est$TeaCNV_mean.2[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.2[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.2[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion2 <- 1-overlap/length(c(base1,base3))
base1 <- est$TeaCNV_mean.3[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.3[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.3[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion3 <- 1-overlap/length(c(base1,base3))
TeaCNV.dispersion <- median(c(TeaCNV.dispersion1,TeaCNV.dispersion2,TeaCNV.dispersion3))
dispersion <- c(inferCNV.dispersion,epiAneuFinder.dispersion,copyscAT.dispersion,TeaCNV.dispersion)
names(dispersion) <- c("inferCNV",method)
barplot(dispersion,ylab = "Dispersion",ylim=c(0,1))
CN.dispersion[[1]] <- dispersion
names(CN.dispersion)[1] <- ID


ID="ccRCC4"
est <- read.csv(paste0("./",ID,"_CNdf_clean_clone.csv"))
est <- est[!est$chr %in% c("chrX","chrY"),]

est$groundtruth <- est$dna_ratio
est$groundtruth[est$dna_CNA=='amp'] <- 3
est$groundtruth[est$dna_CNA=="amp" & est$dna_ratio >1.5] <- 4
est$groundtruth[est$dna_CNA=='neu'] <- 2
est$groundtruth[est$dna_CNA=='del'] <- 1
est <- est[!is.na(est$groundtruth),]
base1 <- est$epiAneuRatio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$epiAneuRatio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$epiAneuRatio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
base4 <- est$epiAneuRatio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
minbase4 <- min(base4)
maxbase4 <- max(base4)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])+length(base4[base4<=maxbase3])++length(base3[base3>=minbase4&base3>maxbase2])
epiAneuFinder.dispersion <- 1-overlap/length(c(base1,base3,base4))

base1 <- est$copyscAT_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$copyscAT_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$copyscAT_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
base4 <- est$epiAneuRatio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
minbase4 <- min(base4)
maxbase4 <- max(base4)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])+length(base4[base4<=maxbase3])+length(base3[base3>=minbase4&base3>maxbase2])
copyscAT.dispersion <- 1-overlap/length(c(base1,base3,base4))

base1 <- est$infercnv_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$infercnv_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$infercnv_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
base4 <- est$infercnv_ratio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
minbase4 <- min(base4)
maxbase4 <- max(base4)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])+length(base4[base4<=maxbase3])+length(base3[base3>=minbase4&base3>maxbase2])
inferCNV.dispersion <- 1-overlap/length(c(base1,base3,base4))

base1 <- est$TeaCNV_mean.1[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.1[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.1[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
base4 <- est$TeaCNV_mean.1[est$groundtruth==4]
base4 <- base3[!is.na(base4)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion1 <- 1-overlap/length(c(base1,base3,base4))
base1 <- est$TeaCNV_mean.2[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.2[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.2[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
base4 <- est$TeaCNV_mean.2[est$groundtruth==4]
base4 <- base3[!is.na(base4)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion2 <- 1-overlap/length(c(base1,base3,base4))
base1 <- est$TeaCNV_mean.3[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.3[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.3[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
base4 <- est$TeaCNV_mean.3[est$groundtruth==4]
base4 <- base3[!is.na(base4)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion3 <- 1-overlap/length(c(base1,base3,base4))
base1 <- est$TeaCNV_mean.4[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.4[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.4[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
base4 <- est$TeaCNV_mean.4[est$groundtruth==4]
base4 <- base3[!is.na(base4)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion4 <- 1-overlap/length(c(base1,base3,base4))
TeaCNV.dispersion <- median(c(TeaCNV.dispersion1,TeaCNV.dispersion2,TeaCNV.dispersion3,TeaCNV.dispersion4))
dispersion <- c(inferCNV.dispersion,epiAneuFinder.dispersion,copyscAT.dispersion,TeaCNV.dispersion)
names(dispersion) <- c("inferCNV",method)
barplot(dispersion,ylab = "Dispersion",ylim=c(0,1))
CN.dispersion[[2]] <- dispersion
names(CN.dispersion)[2] <- ID


ID="ccRCC1"
est <- read.csv(paste0("./",ID,"_CNdf_clean_clone.csv"))
est <- est[!est$chr %in% c("chrX","chrY"),]
est$groundtruth <- est$dna_ratio
est$groundtruth[est$dna_CNA=='amp'] <- 3
est$groundtruth[est$dna_CNA=="amp" & est$dna_ratio >1.5] <- 4
est$groundtruth[est$dna_CNA=='neu'] <- 2
est$groundtruth[est$dna_CNA=='del'] <- 1
est <- est[!is.na(est$groundtruth),]
base1 <- est$epiAneuRatio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$epiAneuRatio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$epiAneuRatio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
base4 <- est$epiAneuRatio[est$groundtruth==4]
base4 <- base4[!is.na(base4)]
minbase4 <- min(base4)
maxbase4 <- max(base4)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])
epiAneuFinder.dispersion <- 1-overlap/length(c(base1,base3))

base1 <- est$copyscAT_ratio[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
minbase1 <- min(base1)
maxbase1 <- max(base1)
base2 <- est$copyscAT_ratio[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
minbase2 <- min(base2)
maxbase2 <- max(base2)
base3 <- est$copyscAT_ratio[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
minbase3 <- min(base3)
maxbase3 <- max(base3)
overlap <- length(base1[base1 >= minbase2])+length(base3[base3<=maxbase2])
copyscAT.dispersion <- 1-overlap/length(c(base1,base3))

base1 <- est$TeaCNV_mean.1[est$groundtruth==1]
base1 <- base1[!is.na(base1)]
base2 <- est$TeaCNV_mean.1[est$groundtruth==2]
base2 <- base2[!is.na(base2)]
base3 <- est$TeaCNV_mean.1[est$groundtruth==3]
base3 <- base3[!is.na(base3)]
base4 <- est$TeaCNV_mean.1[est$groundtruth==4]
base4 <- base3[!is.na(base4)]
overlap <- length(base1[base1 >1])+length(base3[base3<3])
TeaCNV.dispersion <- 1-overlap/length(c(base1,base3))
dispersion <- c(epiAneuFinder.dispersion,copyscAT.dispersion,TeaCNV.dispersion)
names(dispersion) <- method
barplot(dispersion,ylab = "Dispersion",ylim=c(0,1))
CN.dispersion[[3]] <- dispersion
names(CN.dispersion)[3] <- ID

saveRDS(CN.dispersion,file = paste0('./Figure2g.CNV_dirpersion.RDS'))

###Fig 2g: dispersion
pdf(file = paste0('./WGS.CNV.dispersion.pdf'),width=6.5,height = 2.5)
par(mfrow=c(1,3))
barplot(CN.dispersion[[3]],ylab = "Dispersion",ylim=c(0,1),main = 'ccRCC1',las=2,
        col = c('#00A087FF','#8491B4FF', "#B09C85FF"))
barplot(CN.dispersion[[2]],ylab = "Dispersion",ylim=c(0,1),main = 'ccRCC3',las=2,
        col = c('#F39B7FFF','#00A087FF','#8491B4FF', "#B09C85FF"))
barplot(CN.dispersion[[1]],ylab = "Dispersion",ylim=c(0,1),main = 'ccRCC4',las=2,
        col = c('#F39B7FFF','#00A087FF','#8491B4FF', "#B09C85FF"))
dev.off()

##Mean
tb <- do.call(rbind,CN.dispersion[1:2])

tb <- rbind(tb,c(NA,CN.dispersion[[3]]))
rownames(tb)[3] <- names(CN.dispersion)[3]
write.csv(tb,paste0('./Figure2j.CNV_dirpersion.table.csv'),row.names=FALSE)
colMeans(tb,na.rm=TRUE)

