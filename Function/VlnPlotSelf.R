VlnPlotSelf <- function(obj,
                         features = NULL,lim = NULL,round_p = 0 ,
                         ncol = 1){
    require(Seurat);
    require(cowplot);
    require(ggplot2);
    
    # calculate median
    give.m = function(x) {
        la.df = data.frame(y = median(x) + max(x)/8,label = paste0(round(mean(x), round_p)));
        return(la.df)
    }
    P = lapply(features,FUN = function(x){
        
        if(is.null(lim)){
            a = min(obj@meta.data[[x]])
            b = max(obj@meta.data[[x]])
            lim = c(a,b)
        }
        
        Seurat::VlnPlot(obj,features =x,cols = cb_pattern,pt.size = 0,raster=FALSE) +
            stat_summary(fun = mean, geom = "point", col = "black") +
            stat_summary(fun.data = give.m,size = 4,position = 'identity',
                         geom = "text",
                         col = "black")+NoLegend()+ylab('') + ylim(lim[1],lim[2])
    })
    plot_grid_args <- c(P[c(1:length(features))],ncol = ncol)
    do.call(cowplot::plot_grid, plot_grid_args)
}
