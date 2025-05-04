library(ROCR)
library(pROC)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(ggpubr)
library(ggthemes)
library(ggsci)

work_dir <- '~/Library/Mobile Documents/com~apple~CloudDocs/TeaCNVmanuscript/Code_Rdata_upload/'

datapath = paste0(work_dir,'/EvaluationData')
if(!file.exists(datapath)){dir.create(datapath,recursive=T)}
setwd(datapath)

method <- c("epiAneu","copyscAT","TeaCNV")
sampleID <- c("ccRCC1","ccRCC3","ccRCC4")

Psets  <- seq(0.1,1,0.1)

for(j in 1:length(Psets)){
    prop <- Psets[j]
    print(paste0("Process: ", prop))
    ##CNV accuracy
    power.est <- do.call(rbind,lapply(sampleID,function(ID,method,prop){
        est <- read.csv(paste0(datapath,"/",ID,"/",ID,"_CNdf_clean_clone_",prop,".csv"))
        est <- est[!est$chr %in% c("chrX","chrY"),]
        est$epiAneuCN[est$epiAneuCN==1] <- "neu"
        est$copyscAT.CN[est$copyscAT.CN==2] <- "neu"
        if(ID != "ccRCC1"){est$infercnv_CN[est$infercnv_CN==2] <- "neu"}
        
        groundtruth <- est[,1:5] 
        power.est <- do.call(rbind,lapply(method, function(m,est,groundtruth,ID){
            subres <- est[,grep(m,colnames(est))]
            if (m == "TeaCNV"){
                subres <- subres[,grep("ean",colnames(subres))]
                subres <- cbind(groundtruth,subres[,dim(subres)[2]:1])
            }else{
                subres <- cbind(groundtruth,subres)
            }
            subres <- subres[!is.na(subres$dna_CNA),]
            subres$label <- subres$dna_CNA
            subres$label[subres$label!="neu"] <- 1
            subres$label[subres$label=="neu"] <- 0
            if (m == "TeaCNV"){
                kk <- (dim(subres)[2]-6)/2
                power.est <- do.call(rbind,lapply(1:kk, function(i,subres){
                    #subres <- subres[!is.na(subres[,(2*i-1)+5]),]
                    base <- mean(subres[subres$dna_CNA=="neu",(2*i-1)+5][!is.na(subres[subres$dna_CNA=="neu",(2*i-1)+5])])
                    subres$score <- abs(subres[,(2*i-1)+5]-base)
                    subres$predict <- subres[,(2*i)+5]
                    subres$predict[subres[,(2*i)+5]!=2] <- 1
                    subres$predict[subres[,(2*i)+5]==2] <- 0
                    pred <- prediction(subres$score[!is.na(subres$score)], subres$label[!is.na(subres$score)])
                    perfB <- performance(pred,"tpr","fpr")
                    aucB  <- performance(pred, 'auc')
                    aucB  <- unlist(slot(aucB,"y.values"));
                    confusion_matrix <- table(as.character(subres$label),as.character(subres$predict))
                    precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
                    recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
                    f1_score <- 2 * precision * recall/(precision + recall)
                    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
                    power.est <- c(ID,m,aucB,precision,recall,f1_score,accuracy)
                    return(power.est)
                },subres))
            }else{
                #subres <- subres[!is.na(subres[,6]),]
                base <- mean(subres[subres$dna_CNA=="neu",6][!is.na(subres[subres$dna_CNA=="neu",6])])
                subres$score <- abs(subres[,6]-base)
                subres$predict <- subres[,7]
                subres$predict[subres$predict!="neu"] <- 1
                subres$predict[subres$predict=="neu"] <- 0
                pred <- prediction(subres$score[!is.na(subres$score)], subres$label[!is.na(subres$score)])
                perfB <- performance(pred,"tpr","fpr")
                aucB  <- performance(pred, 'auc')
                aucB  <- unlist(slot(aucB,"y.values"));
                confusion_matrix <- table(subres$label,subres$predict)
                if(ncol(confusion_matrix)<2){
                    col_miss <- c("0","1")[!c("0","1") %in% colnames(confusion_matrix)]
                    confusion_matrix <- as.matrix(confusion_matrix)
                    confusion_matrix <- cbind(confusion_matrix,c(0,0))
                    colnames(confusion_matrix)[2] <- col_miss
                }
                precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
                precision[is.na(precision)] <- 0
                recall <- confusion_matrix[2, 2] / sum(subres$label==1)
                f1_score <- 2 * precision * recall/(precision + recall)

                accuracy <- sum(diag(confusion_matrix)) / sum(table(subres$label))               
                power.est <- c(ID,m,aucB,precision,recall,f1_score,accuracy)
            }
            return(power.est)
        },est,groundtruth,ID))
        return(power.est)
    },method,prop))

    power.est <- data.frame(sample = power.est[,1],method = power.est[,2],AUC = as.numeric(power.est[,3]),precision = as.numeric(power.est[,4]),recall = as.numeric(power.est[,5]),F1 = as.numeric(power.est[,6]),ACC = as.numeric(power.est[,7]))
    power.bulk <- power.est[power.est$method!="TeaCNV",]
    for (ID in sampleID){
        subres <- power.est[power.est$method=="TeaCNV"&power.est$sample==ID,]
        subres <- subres[which.max(subres$ACC),]
        power.bulk <- rbind(power.bulk,subres)
    }
    power.est$method <- factor(power.est$method,levels = c("copyscAT","epiAneu","TeaCNV"))




    #####precision and recall
    precision.res <- c()

    ID="ccRCC1"
    est <- read.csv(paste0("./",ID,"/",ID,"_CNdf_clean_clone_",prop,".csv"))
    est <- est[!est$chr %in% c("chrX","chrY"),]
    est$epiAneuCN[est$epiAneuCN==1] <- "neu"
    est$copyscAT.CN[est$copyscAT.CN==2] <- "neu"
    groundtruth <- est[,1:5]
    m="epiAneu"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))/sum(subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))
    recall <- sum(subres$dna_CNA!="neu"&subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="copyscAT"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))/sum(subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))
    recall <- sum(subres$dna_CNA!="neu"&subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="TeaCNV"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$TeaCNV_mean.1!=2&!is.na(subres$TeaCNV_mean.1))/sum(subres$TeaCNV_mean.1!=2&!is.na(subres$TeaCNV_mean.1))
    recall <- sum(subres$dna_CNA!="neu"&subres$TeaCNV_mean.1!=2&!is.na(subres$TeaCNV_mean.1))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    rm(est,ID)

    ID="ccRCC3"
    est <- read.csv(paste0("./",ID,"/",ID,"_CNdf_clean_clone_",prop,".csv"))
    est <- est[!est$chr %in% c("chrX","chrY"),]
    est$epiAneuCN[est$epiAneuCN==1] <- "neu"
    est$copyscAT.CN[est$copyscAT.CN==2] <- "neu"
    est$infercnv_CN[est$infercnv_CN==2] <- "neu"
    groundtruth <- est[,1:5]
    m="epiAneu"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))/sum(subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))
    recall <- sum(subres$dna_CNA!="neu"&subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="copyscAT"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))/sum(subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))
    recall <- sum(subres$dna_CNA!="neu"&subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="infercnv"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$infercnv_CN!="neu"&!is.na(subres$infercnv_CN))/sum(subres$infercnv_CN!="neu"&!is.na(subres$infercnv_CN))
    recall <- sum(subres$dna_CNA!="neu"&subres$infercnv_CN!="neu"&!is.na(subres$infercnv_CN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="TeaCNV"
    subres <- est[,grep(m,colnames(est))]
    predict1 <- apply(subres,1,function(x){
        y <- x[c(2,5,8)]
        y <- y[!is.na(y)]
        if (length(y)!=0){
            if (max(y)==2 & min(y)==2){
                return("neu")
            }else{
                return("cnv")
            }
        }else{
            return(NA)
        }
    })
    predict2 <- apply(subres,1,function(x){
        y <- x[c(2,5,8)]
        y <- y[!is.na(y)]
        if (length(y)!=0){
            if (max(y)!=2 & min(y)!=2){
                return("cnv")
            }else{
                return("neu")
            }
        }else{
            return(NA)
        }
    })
    subres <- cbind(groundtruth,predict1,predict2)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$predict2!="neu"&!is.na(subres$predict2))/sum(subres$predict2!="neu"&!is.na(subres$predict2))
    recall <- sum(subres$dna_CNA!="neu"&subres$predict1!="neu"&!is.na(subres$predict1))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    rm(est,ID)

    ID="ccRCC4"
    est <- read.csv(paste0("./",ID,"/",ID,"_CNdf_clean_clone_",prop,".csv"))
    est <- est[!est$chr %in% c("chrX","chrY"),]
    est$epiAneuCN[est$epiAneuCN==1] <- "neu"
    est$copyscAT.CN[est$copyscAT.CN==2] <- "neu"
    est$infercnv_CN[est$infercnv_CN==2] <- "neu"
    groundtruth <- est[,1:5]
    m="epiAneu"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))/sum(subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))
    recall <- sum(subres$dna_CNA!="neu"&subres$epiAneuCN!="neu"&!is.na(subres$epiAneuCN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="copyscAT"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))/sum(subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))
    recall <- sum(subres$dna_CNA!="neu"&subres$copyscAT.CN!="neu"&!is.na(subres$copyscAT.CN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="infercnv"
    subres <- est[,grep(m,colnames(est))]
    subres <- cbind(groundtruth,subres)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$infercnv_CN!="neu"&!is.na(subres$infercnv_CN))/sum(subres$infercnv_CN!="neu"&!is.na(subres$infercnv_CN))
    recall <- sum(subres$dna_CNA!="neu"&subres$infercnv_CN!="neu"&!is.na(subres$infercnv_CN))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    m="TeaCNV"
    subres <- est[,grep(m,colnames(est))]
    predict1 <- apply(subres,1,function(x){
        y <- x[c(2,5,8)]
        y <- y[!is.na(y)]
        if (length(y)!=0){
            if (max(y)==2 & min(y)==2){
                return("neu")
            }else{
                return("cnv")
            }
        }else{
            return(NA)
        }
    })
    predict2 <- apply(subres,1,function(x){
        y <- x[c(2,5,8)]
        y <- y[!is.na(y)]
        if (length(y)!=0){
            if (max(y)!=2 & min(y)!=2){
                return("cnv")
            }else{
                return("neu")
            }
        }else{
            return(NA)
        }
    })
    subres <- cbind(groundtruth,predict1,predict2)
    subres <- subres[!is.na(subres$dna_CNA),]
    precision <- sum(subres$dna_CNA!="neu"&subres$predict2!="neu"&!is.na(subres$predict2))/sum(subres$predict2!="neu"&!is.na(subres$predict2))
    recall <- sum(subres$dna_CNA!="neu"&subres$predict1!="neu"&!is.na(subres$predict1))/sum(subres$dna_CNA!="neu")
    precision.res <- rbind(precision.res,data.frame(ID = ID,method=m,precision=precision,recall=recall))
    rm(subres,precision,recall)
    rm(est,ID)

    library(ggpubr)
    library(ggthemes)
    library(ggsci)

    write.csv(precision.res,paste0("WGS.precision.recall_",prop,".csv"),row.names = F,quote = F)

    ID1 <- paste0(precision.res$ID,precision.res$method)
    ID2 <- paste0(power.bulk$sample,power.bulk$method)
    index <- match(ID2,ID1)
    power.bulk$precision <- precision.res$precision[index]
    power.bulk$recall <- precision.res$recall[index]
    power.estlong <- melt(setDT(power.bulk), id.vars = c("sample","method"), variable.name = "metrics")
    my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneu", "TeaCNV"))
    write.csv(power.bulk,file = paste0("WGS.power_",prop,".csv"),quote = F,row.names = F)

}

###select prop by 'precision' and 'recall' 
Psets  <- seq(0.1,1,0.1)
power.bulk_all <- c()
for(j in 1:length(Psets)){
    prop <- Psets[j]
    res <- read.csv(paste0("WGS.power_",prop,".csv"))
    res2 <- read.csv(paste0("WGS.precision.recall_",prop,".csv"))
    res2 <- res2[res2$method == "infercnv",]
    res <- full_join(res,res2,by=c('sample'='ID','method','precision','recall'))
    res$prop <- prop
    power.bulk_all <- rbind(power.bulk_all,res)
}
write.csv(power.bulk_all,file = paste0("WGS.power_all.csv"),quote = F,row.names = F)


power.bulk_all <- read.csv(paste0("WGS.power_all.csv"))

power.bulk_all[is.na(power.bulk_all)] <- 0
for(m in c("epiAneu","copyscAT")){
   # m="epiAneu"
    pldata <-power.bulk_all %>%
        dplyr::filter(method==m)%>%
        mutate(prop = factor(prop))%>%
        tidyr::pivot_longer(cols = c("precision","recall"), 
               names_to = "Variable", 
               values_to = "Value")%>%
        as.data.frame()
    pldata <- unique(pldata[,c("sample","prop","Variable","Value")])


    p1 <- ggplot(pldata,aes(x=prop, y=Value,fill=sample))+
        geom_bar(stat = "identity")+
        facet_grid(Variable~sample,switch = "y")+
        theme(
            legend.position="none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        scale_fill_nejm() 
    ggsave(plot = p1,filename = paste0(m,"_ACC_precision_recall_prop.pdf"),path = './',
       width =6,height = 4)

 
}


pldata2 <-power.bulk_all %>%
    dplyr::group_by(sample,method)%>%
    summarise(across(everything(), mean, .names = "{.col}"))%>%
    as.data.frame()
pldata2[pldata2==0] <- NA   
pldata2 <- pldata2[,1:7]
write.csv(pldata2,file = paste0("WGS.power.csv"),quote = F,row.names = F)

#Figure 2h
P1 = ggplot(pldata2, aes(x=recall, y=precision, color = method)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 0.5) + 
    geom_point(size = 2)+
    xlim(0,1)+
    ylim(0,1)+
    theme_few()+
    scale_color_manual(values =  c('#8491B4FF','#00A087FF','#F39B7FFF', "#B09C85FF"))
ggsave(filename = paste0("WGS.precision.recall_propMean.pdf"),P1,
       width =4,height = 2.5)


power.estlong <- melt(setDT(pldata2), id.vars = c("sample","method"), variable.name = "metrics")
power.estlong <- power.estlong[power.estlong$method !="infercnv",,drop=F]
my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneu", "TeaCNV"))

###Fig2i: F1 score
p2 = ggplot(power.estlong, aes(x=method, y=value)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
               aes(color=sample), show.legend = T,size = 2)+
    facet_wrap(~metrics, scale="free",ncol = 5)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    stat_compare_means(comparisons = my_comparisons,method = "t.test",method.args = list(alternative = "less")) +
    theme_test()+
    scale_color_nejm() 
ggsave(plot = p2,filename = paste0("WGS.accuracy.pdf"),path = './',
   width =12,height = 3)


###reProcess data
prop=0.1
for(ID in sampleID){
    est1 <- read.csv(paste0("./",ID,"/",ID,"_CNdf_clean_clone_",prop,".csv"))
    write.csv(est1,file = paste0("./",ID,"_CNdf_clean_clone.csv"),quote = F,row.names = F)
 
}


