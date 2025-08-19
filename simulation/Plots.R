##Visualization of BACC for clonal structure identifization and rare clone identification
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)
library(ggthemes)
library(ggsci) 

outdir <- "./TeaCNV/simulation/RareCloneEva";if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
setwd(outdir)

##(1) rare clone identify
#load data
Plot_for <- "RareClone"
###diffSize
eval_df <- read.csv(paste0(outdir,"/evaRes_TeaCNV_diffSize_RareCloneIdentify.csv"))
###diffCNV
res_perform_rare <- read.csv(paste0("./TeaCNV/simulation/diffCNVEva/data/evaRes_TeaCNV_RareClone_assignment.csv"))
res_perform_rare <- res_perform_rare %>%
 select(-CNVFrac)

df_all <- rbind(eval_df,res_perform_rare)

metric_cols = c("accuracy", "precision", "recall", "F1","BACC")
size_col <- "size"

df_all[is.na(df_all)] <- 0
long_df <- df_all %>%
	dplyr::select(all_of(c("group", size_col,"nth", metric_cols))) %>%
	reshape2::melt(id.vars = c("group", size_col,"nth"), variable.name = "Metric", value.name = "Value")
summary_df <- long_df %>%
	group_by( .data[[size_col]], Metric) %>%
	summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")%>%as.data.frame()
write.csv(summary_df,paste0(outdir,"/summary_RareCloneAssignment.csv"),row.names=FALSE)	

color_box <- "#cccccc"   
color_line <- "#0072B2"
plots <- list()
for (metric_name in metric_cols) {
	plot_df <- long_df %>% dplyr::filter(Metric == metric_name)
	summary_plot_df <- summary_df %>% dplyr::filter(Metric == metric_name)
	y_title <- paste0(metric_name," of rare clone assignment")
	x_title <- "The number of cells in rare clone"
	p <- ggplot(plot_df, aes(x = factor(.data[[size_col]]/10), y = Value)) +
	  geom_boxplot(fill = color_box, alpha = 0.4, outlier.shape = NA, width = 0.6) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[[size_col]]/10), y = mean_value),color = color_line,
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[[size_col]]/10), y = mean_value),color = color_line,
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = "", x = x_title, y = y_title) +
      theme_classic() +
      theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )
   plots[[metric_name]] <- p
}

final_plot <- wrap_plots(plots, ncol = 2)
ggsave(paste0(outdir,"/accuracy_precision_recall_f1_RareClone_Assignment.pdf"),final_plot,width=12,height=9)


P1 = ggplot(df_all, aes(x=recall, y=precision, color = size)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey70", size = 0.5) + 
    geom_point(size = 3, alpha = 0.8)+
    xlim(0,1)+
    ylim(0,1)+
    theme_few()+
    scale_color_viridis_c()+
    labs(title = "Rare clone assignment") +
   theme(
    legend.position = "right",  
    panel.grid.major = element_line(color = "grey95")
  )

ggsave(paste0(outdir,"/precision.recall_RareClone_Assignment.pdf"),P1,width =5,height = 3.5)



###评估rare CNVs事件鉴定指标
library(ggthemes)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggpubr)
library(rstatix)
outdir <- "./TeaCNV/simulation/RareCloneEva"
Plot_for <- "RareCNVidentify"
outdir_data <- "./TeaCNV/simulation/diffCNVEva/data"
res_rare_all_file <- file.path(outdir_data,"res_rareCNVevents_identify_performance_TeaCNV.csv")
res_rare_all_file_epi <- file.path(outdir_data,"res_rareCNVevents_identify_performance_epiAneuFinder.csv")
res_rare_all_file_copy <- file.path(outdir_data,"res_rareCNVevents_identify_performance_copyscAT.csv")

tea <- read.csv(res_rare_all_file)
epi <- read.csv(res_rare_all_file_epi)
copy <- read.csv(res_rare_all_file_copy)

tea <- tea %>% dplyr::select(CNVFrac,group,nth,method,accuracy,precision,recall,F1)

##load diffSize result
res_clonal_file_tea <- "./TeaCNV/simulation/RareCloneEva/res_rareCNVevents_identify_performance_TeaCNV_diffSize.csv"
res_clonal_tea_diffSize <- read.csv(res_clonal_file_tea)
res_clonal_tea_diffSize <- res_clonal_tea_diffSize %>% dplyr::select(CNVFrac=size,group,nth,method,accuracy,precision,recall,F1)


epi[is.na(epi)] <- 0
epi <- epi %>%
	dplyr::select(-method_prop,-clone_truth)%>%
  dplyr::group_by(CNVFrac,group,nth,method)%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE), .names = "{.col}"))%>%
    as.data.frame()

copy[is.na(copy)] <- 0
copy <- copy %>%
	dplyr::select(-method_prop,-clone_truth)%>%
  dplyr::group_by(CNVFrac,group,nth,method)%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE), .names = "{.col}"))%>%
    as.data.frame()
df_all <- rbind(tea,epi,copy,res_clonal_tea_diffSize)

#Plot
col_method <- c(TeaCNV='#c85e62',epiAneuFinder='#67a583',copyscAT='#7b95c6')
darken_color <- function(color, factor=0.7){
    col <- col2rgb(color)
    col <- col * factor
    col <- rgb(t(col), maxColorValue=255)
    return(col)
}
color_box <- "#cccccc"   
color_line <- "#0072B2" 
border_colors <- sapply(col_method, darken_color)

metric_cols1 <- c('accuracy','precision','recall','F1')
long_df <- df_all %>%
  dplyr::select(all_of(c("method", "CNVFrac","group", metric_cols1))) %>%
  reshape2::melt(id.vars = c("method", "CNVFrac","group"), variable.name = "Metric", value.name = "Value")%>%
  as.data.frame()
summary_df <- long_df %>%
  group_by(method, Metric) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
summary_df
write.csv(summary_df,paste0(outdir,"/eva_3method_RareCNVidentify_F1score.summary_all.csv"),row.names=FALSE)

###Global assessment of CNV event identification
my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneuFinder", "TeaCNV"))

plots <- list()
plots2 <- list()
for (metric_name in metric_cols1) {
  plot_df <- long_df %>% dplyr::filter(Metric == metric_name)
  y_title <- paste0(metric_name," of rare CNVs identification")
  p <- ggplot(plot_df, aes(x = method, y = Value,fill=method)) +
          geom_boxplot(aes(color = method),alpha = 0.7, outlier.shape = NA, width = 0.5) +
           # geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.3), 
           #      show.legend = T,size = 1)+
          stat_compare_means(comparisons = my_comparisons,method = "t.test",method.args = list(alternative = "less"), aes(label=..p.adj..)) +
          labs(title = "", x = "", y = y_title) +
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


  plot_df2 <- long_df %>% dplyr::filter(Metric == metric_name,method=="TeaCNV",CNVFrac>1 )
  summary_plot_df <- plot_df2 %>%
	  group_by(CNVFrac) %>%
	  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")

	p2 <- ggplot(plot_df2, aes(x = factor(.data[["CNVFrac"]]), y = Value)) +
	  geom_boxplot(fill= color_box,alpha = 0.6, outlier.shape = NA, width = 0.5) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[["CNVFrac"]]), y = mean_value),color = "#0072B2",
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = " ", x = "The number of cells", y = y_title) +
      theme_classic(base_size = 13) +
      theme(
      	 panel.grid.major = element_line(color = "grey80", size = 0.5), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )+
      # scale_fill_manual(values =  c('#f0f0f0', '#cccccc', '#a8a8a8'))+
      coord_cartesian(ylim = c(0,1), clip = "off")+
     theme(plot.margin = margin(5, 5, 5, 5))
  plots2[[metric_name]] <- p2

}

final_plot <- wrap_plots(plots, ncol = 4)
ggsave(paste0(outdir,"/eva_3method_F1score_",Plot_for,".pdf"),final_plot,width=15,height=4.5)
final_plot2 <- wrap_plots(plots2, ncol = 2)
ggsave(paste0(outdir,"/eva_TeaCNV_F1score_",Plot_for,".pdf"),final_plot2,width=12,height=6)







rm(df_all)

##(2) clonal structure identify
Plot_for <- "ClonalStructure"
###diffSize
df1 <- read.csv(paste0(outdir,"/evaRes_TeaCNV_diffSize_CloneStructureIdentify.csv"))
###diffCNV
res_perform_clonal <- read.csv(paste0("./TeaCNV/simulation/diffCNVEva/data/evaRes_TeaCNV_clonalStructure_assignment.csv"))
res_perform_clonal <- res_perform_clonal %>%
 select(-CNVFrac)
df_all <- rbind(df1,res_perform_clonal)
metric_cols = c("accuracy", "precision", "recall", "F1","BACC")
size_col <- "size"

df_all[is.na(df_all)] <- 0
long_df <- df_all %>%
	dplyr::select(all_of(c("group", size_col,"nth", metric_cols))) %>%
	reshape2::melt(id.vars = c("group", size_col,"nth"), variable.name = "Metric", value.name = "Value")
summary_df <- long_df %>%
	group_by( .data[[size_col]], Metric) %>%
	summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")%>%as.data.frame()
write.csv(summary_df,paste0(outdir,"/summary_",Plot_for,"_Assignment.csv"),row.names=FALSE)	
color_box <- "#cccccc"   
color_line <- "#0072B2"
plots <- list()
for (metric_name in metric_cols) {
	plot_df <- long_df %>% dplyr::filter(Metric == metric_name)
	summary_plot_df <- summary_df %>% dplyr::filter(Metric == metric_name)
	y_title <- paste0(metric_name," of \nclonal structure identification")
	x_title <- "The number of cells"
	p <- ggplot(plot_df, aes(x = factor(.data[[size_col]]), y = Value)) +
	  geom_boxplot(fill = color_box, alpha = 0.4, outlier.shape = NA, width = 0.6) +
	  geom_line(data = summary_plot_df, 
	            aes(x = factor(.data[[size_col]]), y = mean_value),color = color_line,
	            size = 0.7,group = 1, inherit.aes = FALSE) +
	  geom_point(data = summary_plot_df,
	             aes(x = factor(.data[[size_col]]), y = mean_value),color = color_line,
	             size = 1.2, inherit.aes = FALSE) +
	  labs(title = "", x = x_title, y = y_title) +
      theme_classic() +
      theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")
      )
   plots[[metric_name]] <- p
}

final_plot <- wrap_plots(plots, ncol = 2)
ggsave(paste0(outdir,"/accuracy_precision_recall_f1_",Plot_for,"_identification.pdf"),final_plot,width=12,height=9)


P1 = ggplot(df_all, aes(x=recall, y=precision, color = size)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey70", size = 0.5) + 
    geom_point(size = 3, alpha = 0.8)+
    xlim(0,1)+
    ylim(0,1)+
    theme_few()+
    scale_color_viridis_c()+
    labs(title = "Clone Structure Identification") +
   theme(
    legend.position = "right",  
    panel.grid.major = element_line(color = "grey95")
  )

ggsave(paste0(outdir,"/precision.recall_",Plot_for,"_identification.pdf"),P1,width =5,height = 3.5)


###Scatter plot
### Clonal CNVs frequency
library(ggthemes)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggpubr)
library(rstatix)
library(ggExtra) 
library(dplyr)
library(gghalves) 


outdir <- "./TeaCNV/simulation/RareCloneEva"
Plot_for <- "ClonalCNVsfrequency"
outdir_data <- "./TeaCNV/simulation/diffCNVEva/data"
res_tea1_fn <- file.path(outdir_data,"res_clonalCNV_freq_TeaCNV.csv")
res_tea2_fn <- file.path("./TeaCNV/simulation/RareCloneEva/res_clonalCNV_freq_TeaCNV_diffSize.csv")

res_tea1 <- read.csv(res_tea1_fn)
res_tea2 <- read.csv(res_tea2_fn)
colnames(res_tea2)[ncol(res_tea2)] <- "CNVFrac"
res_tea <- rbind(res_tea1,res_tea2)
res_tea <- res_tea[!is.na(res_tea$group),]
res_tea <- droplevels(res_tea)
res_tea <- res_tea%>%
dplyr::filter(!is.na(group))%>%
mutate(size=ifelse(CNVFrac<1,2000,CNVFrac))

plt_dat <- na.omit(res_tea)
write.csv(plt_dat,paste0(outdir_data,"/res_clonalCNV_freq_TeaCNV_total.csv"),row.names=FALSE)




###Comparison of three methods (frequency of clonal CNVs events, difference between prediction and reality)
library(ggpubr)       # stat_compare_means
library(rstatix)      # pairwise comparisons helper

outdir <- "./TeaCNV/simulation/RareCloneEva"
Plot_for <- "ClonalCNVsfrequency_3method"
outdir_data <- "./TeaCNV/simulation/diffCNVEva/data"

res_tea <- read.csv(paste0(outdir_data,"/res_clonalCNV_freq_TeaCNV_total.csv"))
res_copy <- read.csv(paste0(outdir_data,"/res_clonalCNV_freq_copyscAT.csv"))
res_epi <- read.csv(paste0(outdir_data,"/res_clonalCNV_freq_epiAneuFinder.csv"))

res_tea$clone_pred <- NULL
res_copy$size=2000
res_epi$size=2000

res_all <- rbind(res_tea,res_copy,res_epi)
res_all$prop_dist <- as.numeric(res_all$prop_dist)

col_method <- c(TeaCNV='#c85e62',epiAneuFinder='#67a583',copyscAT='#7b95c6')
darken_color <- function(color, factor=0.7){
    col <- col2rgb(color)
    col <- col * factor
    col <- rgb(t(col), maxColorValue=255)
    return(col)
}
border_colors <- sapply(col_method, darken_color)

my_comparisons <- list(c("copyscAT", "TeaCNV"), c("epiAneuFinder", "TeaCNV"))
# Kruskal-Wallis test
kruskal_res <- kruskal.test(prop_dist ~ method, data = res_all)

p <- ggplot(res_all, aes(x = method, y = prop_dist,fill=method)) +
      geom_boxplot(aes(color = method),alpha = 0.7, outlier.shape = NA, width = 0.5) +
	  # geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.3), 
	  #       show.legend = F,size = 1)+
      stat_compare_means(comparisons = my_comparisons,
      	 method = "wilcox.test",        
	    p.adjust.method = "bonferroni", 
	    label = "p.adj",             
	    hide.ns = TRUE) +
      labs(title = "Prediction–Truth Discrepancy in Clonal CNV Proportions", x = "", y = "Absolute Difference") +
        theme_classic(base_size = 13) +
        theme(
          plot.title = element_text(margin = margin(b = 40),size=12),
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(color = "black",angle = 45, hjust = 1,vjust=1),
          axis.title = element_text(face = "bold"),
          axis.line = element_line(color = "black")
        )+
		scale_fill_manual(values =  col_method)+
		scale_color_manual(values=border_colors)+
  coord_cartesian(ylim = c(0,1), clip = "off")+
 theme(plot.margin = margin(40, 5,5, 5))
ggsave(paste0(outdir,"/eva_3method_",Plot_for,"_dist.pdf"),p,width=5,height=6)

method_set <- c("copyscAT","epiAneuFinder","TeaCNV")
for(m in method_set[1:2]){
	plt_dat <- res_all %>%dplyr::filter(method == m)

	plt_dat$prop_pred <- as.numeric(plt_dat$prop_pred)
	plt_dat$prop_truth <- as.numeric(plt_dat$prop_truth)
	plt_dat$Ncell_pred <- as.numeric(plt_dat$Ncell_pred)
	plt_dat$Ncell_truth <- as.numeric(plt_dat$Ncell_truth)
	plt_dat <- plt_dat %>%
	  mutate(
	    quartile_label = dplyr::case_when(
	      prop_truth < 0.3 ~ "<30%",
	      prop_truth <= 0.6 & prop_truth >=0.3 ~ "30%~60%",
	      TRUE ~ ">60%"
	    )
	  )
	  plt_dat$quartile_label <- factor(plt_dat$quartile_label,levels=c("<30%","30%~60%",">60%"))

	pear <- cor.test(plt_dat$prop_pred, plt_dat$prop_truth, method = "pearson")
	label_text <- sprintf("Pearson r = %.2f (p = %.2g)",
	                      pear$estimate, pear$p.value)
	group_colors <- c("#A6CEE3", "#6898C5", "#3A65B0", "#6A3D9A")
	p1 <- ggplot(plt_dat, aes(x = quartile_label, y = prop_pred,color=quartile_label)) +
		geom_half_violin(aes(fill = quartile_label),side = "r",bw=0.05,trim = FALSE, alpha = 0.7,
		 color = NA,scale = "width",position = position_nudge(x = 0.15)) + 
	  theme_few(base_size = 14) +
	   theme(legend.position = "none",                           
	        panel.grid.major.y = element_line(color = "grey90"), 
	        axis.title = element_text(face = "bold")) +  
	   labs(title = "Clonal CNV proportion",x = "Truth", y = "Prediction") +   
	    scale_fill_manual(values = group_colors) + 
	    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
	    annotate("text", x = 3.5, y = 0, label = label_text,
	           hjust = 0, vjust = 1, size = 4, fontface = "italic")+
	    coord_flip()

	ggsave(paste0(outdir,"/ClonalCNVproportion_",m,".pdf"),p1,width =5,height = 4)



	p2 = ggplot(plt_dat, aes(x=as.numeric(prop_pred), y=as.numeric(prop_truth))) +
	    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", size = 0.5) + 
	   geom_density_2d_filled(contour_var = "ndensity",alpha = 0.6)+
	   scale_fill_viridis_d( option = "C",na.value = "white")+
		scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
	  	scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
	    theme_few(base_size = 14)+
	    labs(title = "Clonal CNV proportion",x="Prediction",y="Truth",color = "Density") +
	   theme(
	    legend.position = "right",
	    plot.title = element_text(face = "bold", hjust = 0.5),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    axis.text = element_text(color = "black")
	  )+
	  annotate("text", x = 0.05, y = 0.95, label = label_text,
	           hjust = 0, vjust = 1, size = 4, fontface = "italic")

	ggsave(paste0(outdir,"/ClonalCNVfrequency_",m,"-2Ddensity.pdf"),p2,width =6,height = 5)


}










