#' @title plot_cnv_distribution()
#' @param df Dataframe containing genomic positions (chr, start, end) and values for plot
#' @param chromosome_length Dataframe containing chrom and length
#' @examples
#' # usage
#' plot_cnv_distribution(plot_data, chromosome_length)

plot_cnv_distribution_v0 <- function(df,
	chromosome_length=NULL,
	chromosome_levels=paste0("chr", c(1:22,"X")),
	y_end=1,
	color_by = NULL,
	color_palette = NULL,
	title = "",
	x_lab = "Chromosome",
	y_lab = "-log10(q-value)",
	rect_alpha = 0.2,
	segment_alpha = 0.8
	){

	# data check
	required_cols_chrom <- c("chrom", "length")
	if (!all(required_cols_chrom %in% colnames(chromosome_length))) {
	    stop("'chromosome_length' must contain the following columns: ", paste(required_cols_chrom, collapse = ", "))
	}
	chromosome_length <- chromosome_length %>%
		mutate(
			cumulative_start = lag(cumsum(as.numeric(length)), default = 0),
			line_pos = cumulative_start+length/2,
			cumulative_end=cumulative_start+length
		)


	required_cols <- c("chr", "start","end")
	if (!all(required_cols %in% colnames(df))) {
	    stop("'df' must contain the following columns: ", paste(required_cols, collapse = ", "))
	}
	df <- df %>%
	  mutate(chr = factor(chr, levels = chromosome_levels)) %>%
	  arrange(chr, start) %>%
	  left_join(chromosome_length, by = c("chr"="chrom"))	
	df <- df %>%
	  mutate(
	    Start_abs = start + cumulative_start,
	    End_abs = end + cumulative_start
	  ) %>%
	  dplyr::select(-cumulative_start,-length) 
	df$chr <- factor(df$chr, levels = chromosome_levels)

	plot_data <- df
	gg <- ggplot(plot_data)
	gg <- gg +
	    geom_rect(
	      data = chromosome_length,
	      aes(xmin = cumulative_start, 
	          xmax = cumulative_end,
	          ymin = -Inf, 
	          ymax = Inf),
	      fill = rep(c("gray80", "white"), length.out = nrow(chromosome_length)),
	      alpha = rect_alpha,
	      show.legend = FALSE
	    )
	# 动态处理颜色映射
	if (!is.null(color_by)) {
	    # 校验分组列存在性
	    if (!color_by %in% colnames(plot_data)) {
	      stop("color_by指定的列不存在于plot_data中")
	    }
	    
	    # 设置默认颜色映射
	    if (is.null(color_palette)) {
	      color_palette <- scales::hue_pal()(n_distinct(plot_data[[color_by]]))
	      names(color_palette) <- unique(plot_data[[color_by]])
	    } 

	    # 添加带颜色的线段层
    gg <- gg +
      lapply(c("Start_abs", "End_abs"), function(pos) {
        geom_segment(
          aes_string(x = pos, 
                    xend = pos,
                    y = 0, 
                    yend = y_end,
                    color = color_by),
          alpha = segment_alpha
        )
      }) +
      geom_segment(
        aes_string(x = "Start_abs",
                  xend = "End_abs",
                  y = y_end,
                  yend = y_end,
                  color = color_by),
        alpha = segment_alpha
      ) +
      scale_color_manual(values = color_palette)
    
  } else {
    # 无颜色分组模式
    default_color <- ifelse(is.null(color_palette), "red", color_palette[1])
    
    gg <- gg +
      lapply(c("Start_abs", "End_abs"), function(pos) {
        geom_segment(
          aes_string(x = pos, 
                    xend = pos,
                    y = 0, 
                    yend = y_end),
          color = default_color,
          alpha = segment_alpha
        )
      }) +
      geom_segment(
        aes_string(x = "Start_abs",
                  xend = "End_abs",
                  y = y_end,
                  yend = y_end),
        color = default_color,
        alpha = segment_alpha
      )
  }
  # 通用图形设置
  gg <- gg +
    scale_x_continuous(
      breaks = chromosome_length$line_pos,
      labels = gsub("chr|chr19|chr21", "", levels(plot_data$chr)) # 需确保plot_data有Chromosome因子列
    ) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  return(gg)

}


plot_cnv_distribution <- function(df,
	chromosome_length=NULL,
	chromosome_levels=paste0("chr", c(1:22,"X")),
	y_end=1,
	color_by = NULL,
	color_palette = NULL,
	title = "",
	x_lab = "Chromosome",
	y_lab = "-log10(q-value)",
	rect_alpha = 0.2,
	segment_alpha = 0.8,
	gene_annotations = NULL,
	annot_y = 1.2, 
  annot_color = "darkred",
  annot_line_type = "dashed",
  text_size = 3
	){



	# data check
	required_cols_chrom <- c("chrom", "length")
	if (!all(required_cols_chrom %in% colnames(chromosome_length))) {
	    stop("'chromosome_length' must contain: ", paste(required_cols_chrom, collapse = ", "))
	}
	chromosome_length <- chromosome_length %>%
		mutate(
			cumulative_start = lag(cumsum(as.numeric(length)), default = 0),
			line_pos = cumulative_start+length/2,
			cumulative_end=cumulative_start+length
		)

	if (!is.null(gene_annotations)) {
    # 校验基因注释数据
    required_gene_cols <- c("gene", "chr", "start", "end")
    if (!all(required_gene_cols %in% colnames(gene_annotations))) {
      stop("gene_annotations must contain: ", paste(required_gene_cols, collapse = ", "))
    }
    
    # 转换染色体因子并合并位置信息
    gene_annot <- gene_annotations %>%
      mutate(
        chr = factor(chr, levels = chromosome_levels),
        gene_mid = start + (end - start)/2
      ) %>%
      left_join(
        chromosome_length %>% select(chrom, cumulative_start),
        by = c("chr" = "chrom")
      ) %>%
      mutate(
        gene_mid_abs = gene_mid + cumulative_start,
        y_start = annot_y,
        y_end = y_end
      )
  }


	required_cols <- c("chr", "start","end")
	if (!all(required_cols %in% colnames(df))) {
	    stop("'df' must contain: ", paste(required_cols, collapse = ", "))
	}
	df <- df %>%
	  mutate(chr = factor(chr, levels = chromosome_levels)) %>%
	  arrange(chr, start) %>%
	  left_join(chromosome_length, by = c("chr"="chrom"))	
	df <- df %>%
	  mutate(
	    Start_abs = start + cumulative_start,
	    End_abs = end + cumulative_start
	  ) %>%
	  dplyr::select(-cumulative_start,-length) 
	df$chr <- factor(df$chr, levels = chromosome_levels)

	plot_data <- df
	gg <- ggplot(plot_data)
	gg <- gg +
	    geom_rect(
	      data = chromosome_length,
	      aes(xmin = cumulative_start, 
	          xmax = cumulative_end,
	          ymin = -Inf, 
	          ymax = Inf),
	      fill = rep(c("gray80", "white"), length.out = nrow(chromosome_length)),
	      alpha = rect_alpha,
	      show.legend = FALSE
	    )
	# 动态处理颜色映射
	if (!is.null(color_by)) {
	    # 校验分组列存在性
	    if (!color_by %in% colnames(plot_data)) {
	      stop("color_by指定的列不存在于plot_data中")
	    }
	    
	    # 设置默认颜色映射
	    if (is.null(color_palette)) {
	      color_palette <- scales::hue_pal()(n_distinct(plot_data[[color_by]]))
	      names(color_palette) <- unique(plot_data[[color_by]])
	    } 

	    # 添加带颜色的线段层
    gg <- gg +
      lapply(c("Start_abs", "End_abs"), function(pos) {
        geom_segment(
          aes_string(x = pos, 
                    xend = pos,
                    y = 0, 
                    yend = y_end,
                    color = color_by),
          alpha = segment_alpha
        )
      }) +
      geom_segment(
        aes_string(x = "Start_abs",
                  xend = "End_abs",
                  y = y_end,
                  yend = y_end,
                  color = color_by),
        alpha = segment_alpha
      ) +
      scale_color_manual(values = color_palette)
    
  } else {
    # 无颜色分组模式
    default_color <- ifelse(is.null(color_palette), "red", color_palette[1])
    
    gg <- gg +
      lapply(c("Start_abs", "End_abs"), function(pos) {
        geom_segment(
          aes_string(x = pos, 
                    xend = pos,
                    y = 0, 
                    yend = y_end),
          color = default_color,
          alpha = segment_alpha
        )
      }) +
      geom_segment(
        aes_string(x = "Start_abs",
                  xend = "End_abs",
                  y = y_end,
                  yend = y_end),
        color = default_color,
        alpha = segment_alpha
      )
  }

  if (!is.null(gene_annotations)) {
    gg <- gg + 
      # 添加引线
      geom_segment(
        data = gene_annot,
        aes(x = gene_mid_abs, xend = gene_mid_abs,
            y = y_start, yend = y_end),
        color = annot_color,
        linetype = annot_line_type,
        alpha = 0.6
      ) +
      # 添加文本标签
      geom_text(
        data = gene_annot,
        aes(x = gene_mid_abs, y = y_start + 0.1,
            label = gene),
        color = annot_color,
        size = text_size,
        angle = 90,
        hjust = 0,
        vjust = 0.5
      )
  }

  # 通用图形设置
  gg <- gg +
    scale_x_continuous(
      breaks = chromosome_length$line_pos,
      labels = gsub("chr|chr19|chr21", "", levels(plot_data$chr)) # 需确保plot_data有Chromosome因子列
    ) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  

  return(gg)

}


