suppressMessages(suppressWarnings({library(patchwork)}))
suppressMessages(suppressWarnings({library(ggrepel)}))
suppressMessages(suppressWarnings({library(dplyr)}))
suppressMessages(suppressWarnings({library(ggplot2)}))
suppressMessages(suppressWarnings({library(Signac)}))
suppressMessages(suppressWarnings({library(Seurat)}))
suppressMessages(suppressWarnings({library(glue)}))
suppressMessages(suppressWarnings({library(ggsci)}))
myFindRegion <- function (object, region, sep = c("-", "-"), assay = NULL, extend.upstream = 0, 
                         extend.downstream = 0) 
{
  if (!is(object = region, class2 = "GRanges")) {
    region <- tryCatch(expr = suppressWarnings(expr = StringToGRanges(regions = region, 
                                                                      sep = sep)), error = function(x) {
                                                                        region <- LookupGeneCoords(object = object, assay = assay, 
                                                                                                   gene = region)
                                                                        return(region)
                                                                      })
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(x = region, upstream = extend.upstream, 
                                           downstream = extend.downstream))
  return(region)
}

myCoveragePlotSingle <- function(obj_cells,obj_peaks, region, 
                                     extend.upstream = 0, 
                                     extend.downstream = 0,
                                     ymax=NULL,
                                     window=100,heights = NULL,color_fill=NULL){
  # plot
  rg <- myFindRegion(object = obj_cells, region = region, 
                     sep = c("-", "-"),
                     extend.upstream = extend.upstream,
                     extend.downstream = extend.downstream)
  a1=CoveragePlot(obj_cells, 
                  region = rg,
                  ymax=ymax,
                  window=window, 
                  extend.upstream = 0, 
                  extend.downstream = 0,
                  annotation = FALSE,
                  peaks = FALSE)
  if(!is.null(color_fill)){a1 <- a1+scale_fill_manual(values=color_fill)}

  a2=AnnotationPlot(obj_cells, rg, extend.upstream = 0, extend.downstream = 0)
  a3=PeakPlot(obj_peaks, region=rg)
  return(CombineTracks(plotlist = list(a1,a2,a3),heights = NULL))
}


myCoveragePlotMultiple <- function(obj_cells,obj_peaks, region, 
                               ymax=NULL,window=2000, heights = NULL,
                               extend.upstream = 0, 
                               extend.downstream = 0,color_fill=NULL,...){
  if (length(x = region) == 1) {
    region <- list(region)
  }
  plot.list = lapply(X = seq_along(region), FUN = function(x) {
    myCoveragePlotSingle(obj_cells,obj_peaks, region[[x]], 
                             extend.upstream = extend.upstream, 
                             extend.downstream = extend.downstream,
                             ymax=ymax,window=window,heights = heights,color_fill=color_fill)
  })
  return(wrap_plots(plot.list, ...))
}

myCoverageModify<-function(plot_res, 
                           region,
                           ncol=NULL,
                           group_color=NULL, 
                           genes_color='black', 
                           peaks_color='red',
                           line_color_with_group=TRUE,...){
  if(is.null(ncol)){
    row_first = c(1)
  }else{
    row_first = seq(1,length(region),ncol)
  }
  if(is.null(group_color)){
    fill_theme = scale_fill_igv()
    color_theme = scale_color_igv()
  }else{
    fill_theme = scale_fill_manual(values=group_color)
    color_theme = scale_color_manual(values=group_color)
  }
  for(gene_plot in 1:length(plot_res)){
    tmp_layer = plot_res[[gene_plot]][[2]]$layers
    gene_title=region[gene_plot]
    for(ix in 1:length(tmp_layer)){
      if('GeomText'%in%class(tmp_layer[[ix]]$geom)){
        y = tmp_layer[[ix]]$data
        gene_title = y[which.max(y$width), 'gene_name']
      }
    }
    # gene_names = plot_res[[gene_plot]][[2]]$layers[[4]]$data$gene_name
    # gene_names = intersect(region, gene_names)
    # if(length(gene_names)>=1){
    #   gene_title=gene_names
    # }else{
    #   gene_title=''
    # }
    #print(gene_title)
    if(gene_plot%in%row_first){
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          #theme_void()+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            strip.text = element_text(size=8, face="bold"),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size=8, face="bold"),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
        if('GeomText'%in%class(a2$layers[[ix]]$geom)){
          label_ix=ix
        }else{
          if(ix==1){
            a2$layers[[ix]]$aes_params$size=2
          }else{
            a2$layers[[ix]]$aes_params$size=0.1
          }
        }
        
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              axis.title.y = element_text(angle = 0, vjust = 0.5,
                                          size=8, face="bold",
                                          margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      
      plot_res[[gene_plot]][[2]] = a2
      
        
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      if(length(a3$layers)>0){
          a3$layers[[1]]$aes_params$size=3
      }
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5,
                                      size=8, face="bold",
                                      margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      plot_res[[gene_plot]][[3]] = a3
    }else{
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
        if('GeomText'%in%class(a2$layers[[ix]]$geom)){
          label_ix=ix
        }else{
          if(ix==1){
            a2$layers[[ix]]$aes_params$size=2
          }else{
            a2$layers[[ix]]$aes_params$size=0.1
          }
        }
        
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              #axis.title.y = element_text(angle = 0, vjust = 0.5,
              #                            margin=margin(l=0,r=-4,t=0,b=0, unit = 'line')),
              axis.title.y=element_blank()
        )
      
      plot_res[[gene_plot]][[2]] = a2
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      if(length(a3$layers)>0){
          a3$layers[[1]]$aes_params$size=3
      }
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.y=element_blank()
        )
      plot_res[[gene_plot]][[3]] = a3
    }
  }
  wp = list()
  for(i in 1:length(plot_res)){
    wp[[i]] = plot_res[[i]]
  }
  return(wrap_plots(wp, ncol = ncol))
}
myCoveragePlot <- function(obj, region, group.by, ymax='q75',
                           window=2000, heights = c(20,1,1),
                           ncol=NULL,
                           group_color=NULL, genes_color='black', peaks_color='red',line_color_with_group=TRUE,...
){
  plot_res =CoveragePlot(obj, region = region, 
                         group.by = group.by,
                         ymax=ymax, window=window, heights = heights,ncol=ncol,...)
  if(is.null(ncol)){
    row_first = c(1)
  }else{
    row_first = seq(1,length(region),ncol)
  }
  if(is.null(group_color)){
    fill_theme = scale_fill_igv()
    color_theme = scale_color_igv()
  }else{
    fill_theme = scale_fill_manual(group_color)
    color_theme = scale_color_manual(group_color)
  }
  for(gene_plot in 1:length(plot_res)){
    tmp_layer = plot_res[[gene_plot]][[2]]$layers
    gene_title=region[gene_plot]
    for(ix in 1:length(tmp_layer)){
        if('GeomText'%in%class(tmp_layer[[ix]]$geom)){
          y = tmp_layer[[ix]]$data
          gene_title = y[which.max(y$width), 'gene_name']
        }
    }
    # gene_names = plot_res[[gene_plot]][[2]]$layers[[4]]$data$gene_name
    # gene_names = intersect(region, gene_names)
    # if(length(gene_names)>=1){
    #   gene_title=gene_names
    # }else{
    #   gene_title=''
    # }
    #print(gene_title)
    if(gene_plot%in%row_first){
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          #theme_void()+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            strip.text = element_text(size=8, face="bold"),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(
            panel.spacing.y = unit(x = 0, units = "line"),
            strip.placement = 'outside',
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size=8, face="bold"),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
          if('GeomText'%in%class(a2$layers[[ix]]$geom)){
              label_ix=ix
          }else{
              if(ix==1){
                  a2$layers[[ix]]$aes_params$size=2
              }else{
                  a2$layers[[ix]]$aes_params$size=0.1
              }
          }
      
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              axis.title.y = element_text(angle = 0, vjust = 0.5,
                                          size=8, face="bold",
                                          margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      
      plot_res[[gene_plot]][[2]] = a2
      
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      a3$layers[[1]]$aes_params$size=3
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5,
                                      size=8, face="bold",
                                      margin=margin(l=0,r=-2,t=0,b=0, unit = 'line'))
        )
      plot_res[[gene_plot]][[3]] = a3
    }else{
      # coverage 1
      y_upper = max(plot_res[[gene_plot]][[1]]$data$coverage)
      if(line_color_with_group){
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(aes(yintercept=0,color=group), size=0.8)+
          color_theme+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }else{
        plot_res[[gene_plot]][[1]] = plot_res[[gene_plot]][[1]]+
          fill_theme+
          scale_y_continuous(breaks=0,sec.axis = dup_axis(breaks = 0),
                             expand=c(0,0))+
          geom_hline(yintercept = y_upper, size=0.5)+
          geom_hline(yintercept=0, size=0.5)+
          labs(title=gene_title)+
          theme(#panel.border = element_rect(fill=NA),
            #panel.spacing = unit(-1,'lines'),
            panel.spacing.y = unit(x = 0, units = "line"),
            #strip.placement = 'outside',
            strip.text.y.left = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(size = 0.2),
            plot.title = element_text(hjust=0.5,size=12, face="bold.italic")
          )
      }
      
      # gene
      a2 = plot_res[[gene_plot]][[2]]
      label_ix = 4
      for(ix in 1:length(a2$layers)){
        if('GeomText'%in%class(a2$layers[[ix]]$geom)){
          label_ix=ix
        }else{
          if(ix==1){
            a2$layers[[ix]]$aes_params$size=2
          }else{
            a2$layers[[ix]]$aes_params$size=0.1
          }
        }
        
      }
      # a2$layers[[1]]$aes_params$size=2
      # a2$layers[[2]]$aes_params$size=0.1
      # a2$layers[[3]]$aes_params$size=0.1
      a2$layers[[label_ix]] = NULL
      #gene_y_lims = layer_scales(a2)$y$range$range
      tmp_da = a2$layer[[1]]$data
      dodge_y = as.numeric(tmp_da[tmp_da$gene_name==gene_title, 'dodge'][1])
      gene_y_lims = c(dodge_y-0.1,dodge_y+0.1)
      a2=a2+
        scale_color_manual(values = c(genes_color))+
        scale_size_manual(values=0.1)+
        scale_y_continuous(limits = gene_y_lims)+
        theme(panel.border = element_rect(fill=NA, size=0.5),
              strip.placement = 'inside',
              axis.line = element_blank(),
              #axis.title.y = element_text(angle = 0, vjust = 0.5,
              #                            margin=margin(l=0,r=-4,t=0,b=0, unit = 'line')),
              axis.title.y=element_blank()
        )
      
      plot_res[[gene_plot]][[2]] = a2
      
      # peaks 3
      a3 = plot_res[[gene_plot]][[3]]
      a3$layers[[1]]$aes_params$size=3
      a3 =a3+
        scale_color_manual(values = c(peaks_color))+
        #scale_y_continuous(limits = c(-0,0))+
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_rect(fill=NA, size=0.5),
          strip.placement = 'inside',
          axis.line = element_blank(),
          axis.title.y=element_blank()
        )
      plot_res[[gene_plot]][[3]] = a3
    }
  }
  wp = list()
  for(i in 1:length(plot_res)){
    wp[[i]] = plot_res[[i]]
  }
  return(wrap_plots(wp, ncol = ncol))
}



#' current cell identities
#' @param idents A list of identities to include in the plot. If NULL, include
#' all identities
#' @param slot Which slot to pull expression data from
#'
#' @importFrom SeuratObject GetAssayData DefaultAssay
#' @importFrom ggplot2 ggplot geom_violin facet_wrap aes theme_classic theme
#' element_blank scale_y_discrete scale_x_continuous scale_fill_manual
#' @importFrom scales hue_pal
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges start end
#' @importFrom patchwork wrap_plots
#' @importFrom fastmatch fmatch
#'
#' @export
#' @concept visualization
#' @examples
#' \donttest{
#' ExpressionPlot(atac_small, features = "TSPAN6", assay = "RNA")
#' }


ExpressionPlot_gradient <- function(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  slot = "data",
  fill_by_group=FALSE,
  theme_gradient=NULL
) {

  # get data
  assay <- Signac:::SetIfNull(x = assay, y = DefaultAssay(object = object))
  data.plot <- GetAssayData(
    object = object,
    assay = assay,
    layer = slot
  )[features, ]
  obj.groups <- Signac:::GetGroups(
    object = object,
    group.by = group.by,
    idents = NULL
  )
  obj.groups <- obj.groups[colnames(object[[assay]])] 
  # if levels set, define colors based on all groups
  levels.use <- levels(x = obj.groups)
  if (!is.null(x = levels.use) &fill_by_group) {
    colors_all <- scales::hue_pal()(length(x = levels.use))
    names(x = colors_all) <- levels.use
  }
  if (!is.null(x = idents)) {
    cells.keep <- names(x = obj.groups)[
      fmatch(x = obj.groups, table = idents, nomatch = 0L) > 0
    ]
    if (length(x = features) > 1) {
      data.plot <- data.plot[, cells.keep]
    } else {
      data.plot <- data.plot[cells.keep]
    }
    obj.groups <- obj.groups[cells.keep]
  }
   # construct data frame
  if (length(x = features) == 1) {
    df <- data.frame(
      gene = features,
      expression = data.plot,
      group = obj.groups
    )
  } else {
    df <- data.frame()
    for (i in features) {
      df.1 <- data.frame(
        gene = i,
        expression = data.plot[i, ],
        group = obj.groups
      )
      df <- rbind(df, df.1)
    }
  }
  missing.levels <- setdiff(x = levels(x = df$group), y = unique(x = df$group))
  if (!is.null(x = idents)) {
    missing.levels <- intersect(x = missing.levels, y = idents)
  }
  if (length(x = missing.levels) > 0) {
    # fill missing idents with zero
    for (i in features) {
      df.1 <- data.frame(
        gene = i,
        expression = 0,
        group = missing.levels
      )
      df <- rbind(df, df.1)
    }
  }
  p.list <- list()
  lower.limit <- ifelse(test = slot == "scale.data", yes = NA, no = 0)
  for (i in seq_along(along.with = features)) {
    df.use <- df[df$gene == features[[i]], ]
    if(fill_by_group){
      p <- ggplot(data = df.use, aes(x = expression, y = gene, fill = group)) +
      geom_violin(size = 1/4) +
      facet_wrap(~group, ncol = 1, strip.position = "right") +
      theme_classic() +
      scale_y_discrete(position = "top") +
      scale_x_continuous(position = "bottom", limits = c(lower.limit, NA)) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none"
      )
      if (!is.null(x = levels.use)) {
        p <- p + scale_fill_manual(values = colors_all)
      }
    }else{
      df.use <- df.use %>%
        group_by(group)%>%
        dplyr::mutate(group_mean = mean(expression,na.rm=TRUE))
      p <- ggplot(data = df.use, aes(x = expression, y = gene)) +
        geom_violin(aes(fill=group_mean),scale = "width",size = 1/4) +
        facet_wrap(~group, ncol = 1, strip.position = "right") +
        theme_classic() +
        scale_y_discrete(position = "top") +
        scale_x_continuous(position = "bottom", limits = c(lower.limit, NA)) +
        theme(
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "none"
        )
      if (!is.null(theme_gradient) ) {

        p <- p + theme_gradient
      }else{
        p <- p + scale_fill_gradientn(colours = colorRampPalette(c('#F4CAB4','#DB967F', '#B75D55','#922935'))(10))
        #scale_fill_gradient(limits=c(min(df.use$group_mean),max(df.use$group_mean)),low = "yellow",high = "red")
      }

    }

    p.list[[i]] <- p
  }
  p <- wrap_plots(p.list, ncol = length(x = p.list))
  return(p)
}













