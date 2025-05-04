
#' List of color palettes that can be used in plots
#' 
#' A collection of some original and some borrowed color palettes for plots
#' 
#' @export
colorPalettes <- list(
  
  #DISCLOSURE: This is a collection of palettes that includes some original palettes and some palettes originally
  #implemented by others in other packages.
  #They are included here for convenience because they help improve plot aesthetics.
  
  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------
 
   #40-colors
  colour40 = c("1"="#DC143C","2"="#0000FF","3"="#20B2AA","4"="#FFA500","5"="#9370DB","6"="#98FB98","7"="#F08080","8"="#1E90FF","9"="#7CFC00","10"="#FFFF00",
               "11"="#808000","12"="#FF00FF","13"="#FA8072","14"="#7B68EE","15"="#9400D3","16"="#800080","17"="#A0522D","18"="#D2B48C","19"="#D2691E","20"="#87CEEB",
               "21"="#40E0D0","22"="#5F9EA0","23"="#FF1493","24"="#0000CD","25"="#008B8B","26"="#FFE4B5","27"="#8A2BE2","28"="#228B22","29"="#E9967A","30"="#4682B4",
               "31"="#32CD32","32"="#F0E68C","33"="#FFFFE0","34"="#EE82EE","35"="#FF6347","36"="#6A5ACD","37"="#9932CC","38"="#8B008B","39"="#8B4513","40"="#DEB887"),
  #20-colors
  stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D"),
  
  stallion2 = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767"),
  
  calm = c("1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
           "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36", "15"="#D4477D", "16"="#403552", "17"="#76D73C", "18"="#96CED5", "19"="#CE54D1", "20"="#C48736"),
  
  kelly = c("1"="#FFB300", "2"="#803E75", "3"="#FF6800", "4"="#A6BDD7", "5"="#C10020", "6"="#CEA262", "7"="#817066", "8"="#007D34", "9"="#F6768E", "10"="#00538A",
            "11"="#FF7A5C", "12"="#53377A", "13"="#FF8E00", "14"="#B32851", "15"="#F4C800", "16"="#7F180D", "17"="#93AA00", "18"="#593315", "19"="#F13A13", "20"="#232C16"),
  
  #16-colors
  bear = c("1"="#faa818", "2"="#41a30d","3"="#fbdf72", "4"="#367d7d",  "5"="#d33502", "6"="#6ebcbc", "7"="#37526d",
           "8"="#916848", "9"="#f5b390", "10"="#342739", "11"="#bed678","12"="#a6d9ee", "13"="#0d74b6",
           "14"="#60824f","15"="#725ca5", "16"="#e0598b"),
  
  #15-colors
  ironMan = c("9"='#371377',"3"='#7700FF',"2"='#9E0142',"10"='#FF0080', "14"='#DC494C',"12"="#F88D51","1"="#FAD510","8"="#FFFF5F","4"='#88CFA4',
              "13"='#238B45',"5"="#02401B", "7"="#0AD7D3","11"="#046C9A", "6"="#A2A475", "15"='grey35'),
  
  circus = c("1"="#D52126", "2"="#88CCEE", "3"="#FEE52C", "4"="#117733", "5"="#CC61B0", "6"="#99C945", "7"="#2F8AC4", "8"="#332288",
             "9"="#E68316", "10"="#661101", "11"="#F97B72", "12"="#DDCC77", "13"="#11A579", "14"="#89288F", "15"="#E73F74"),
  
  #12-colors
  paired = c("9"="#A6CDE2","1"="#1E78B4","3"="#74C476","12"="#34A047","11"="#F59899","2"="#E11E26",
             "10"="#FCBF6E","4"="#F47E1F","5"="#CAB2D6","8"="#6A3E98","6"="#FAF39B","7"="#B15928"),
  
  #11-colors
  grove = c("11"="#1a1334","9"="#01545a","1"="#017351","6"="#03c383","8"="#aad962","2"="#fbbf45","10"="#ef6a32","3"="#ed0345","7"="#a12a5e","5"="#710162","4"="#3B9AB2"),
  
  #7-colors
  summerNight = c("1"="#2a7185", "2"="#a64027", "3"="#fbdf72","4"="#60824f","5"="#9cdff0","6"="#022336","7"="#725ca5"),
  
  #5-colors
  zissou = c("1"="#3B9AB2", "4"="#78B7C5", "3"="#EBCC2A", "5"="#E1AF00", "2"="#F21A00"), #wesanderson
  darjeeling = c("1"="#FF0000", "2"="#00A08A", "3"="#F2AD00", "4"="#F98400", "5"="#5BBCD6"), #wesanderson
  rushmore = c("1"="#E1BD6D", "5"="#EABE94", "2"="#0B775E", "4"="#35274A" , "3"="#F2300F"), #wesanderson
  captain = c("1"="grey","2"="#A1CDE1","3"="#12477C","4"="#EC9274","5"="#67001E"),
  
  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------
  
  #10-colors
  horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60'),
  
  #9-colors
  horizonExtra =c("1"="#000436","4"="#021EA9","6"="#1632FB","8"="#6E34FC","3"="#C732D5","9"="#FD619D","7"="#FF9965","5"="#FFD32B","2"="#FFFC5A"),
  blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D"),
  sambaNight = c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500'), #buencolors
  solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D'),  #buencolors
  whitePurple = c("9"='#f7fcfd',"6"='#e0ecf4',"8"='#bfd3e6',"5"='#9ebcda',"2"='#8c96c6',"4"='#8c6bb1',"7"='#88419d',"3"='#810f7c',"1"='#4d004b'),
  whiteBlue = c("9"='#fff7fb',"6"='#ece7f2',"8"='#d0d1e6',"5"='#a6bddb',"2"='#74a9cf',"4"='#3690c0',"7"='#0570b0',"3"='#045a8d',"1"='#023858'),
  whiteRed = c("1"="white", "2"="red"),
  comet = c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"),
  
  #7-colors
  greenBlue = c("4"='#e0f3db',"7"='#ccebc5',"2"='#a8ddb5',"5"='#4eb3d3',"3"='#2b8cbe',"6"='#0868ac',"1"='#084081'),
  
  #6-colors
  beach = c("4"="#87D2DB","1"="#5BB1CB","6"="#4F66AF","3"="#F15F30","5"="#F7962E","2"="#FCEE2B"),
  
  #5-colors
  coolwarm = c("1"="#4858A7", "4"="#788FC8", "5"="#D6DAE1", "3"="#F49B7C", "2"="#B51F29"),
  fireworks = c("5"="white","2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  greyMagma = c("2"="grey", "4"="#FB8861FF", "5"="#B63679FF", "3"="#51127CFF", "1"="#000004FF"),
  fireworks2 = c("5"="black", "2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  purpleOrange = c("5"="#581845", "2"="#900C3F", "4"="#C70039", "3"="#FF5744", "1"="#FFC30F")
)

#' Optimized discrete color palette generation
#'
#' This function assesses the number of inputs and returns a discrete color palette that is tailored to provide the most
#' possible color contrast from the designated color set.
#'
#' @param values A character vector containing the sample names that will be used. Each entry in this character vector will be
#' given a unique color from the designated palette set.
#' @param set The name of a color palette.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteDiscrete <- function(
  values = NULL,
  set = "stallion",  
  reverse = FALSE
){
  
  values <- gtools::mixedsort(values)
  n <- length(unique(values))
  pal <- colorPalettes[[set]]
  palOrdered <- pal[gtools::mixedsort(names(pal))] #mixed sort gets 1,2,3,4..10,11,12
  
  if(n > length(palOrdered)){
    message("Length of unique values greater than palette, interpolating..")
    palOut <- colorRampPalette(pal)(n)
  }else{
    palOut <- palOrdered[seq_len(n)]
  }
  
  if(reverse){
    palOut <- rev(palOut)
  }
  
  names(palOut) <- unique(values)
  
  return(palOut)
  
}

#' Continuous Color Palette
#'
#' @param set The name of a color palette provided in the list object.
#' @param n The number of unique colors to generate as part of this continuous color palette.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteContinuous <- function(
  set = "solarExtra", 
  n = 256, 
  reverse = FALSE
){

  pal <- colorPalettes[[set]]
  palOut <- colorRampPalette(pal)(n)
  
  if(reverse){
    palOut <- rev(palOut)
  }
  
  return(palOut)
  
}

library(Seurat)
library(RColorBrewer)

#color setting
VsampleColor = brewer.pal(n = 9, name = "Blues")
RedColor = brewer.pal(n = 9, name = "Reds")
#25-colors
myColor = c("red3","steelblue4","purple3","green3","brown","gold1","cyan3", "darkorange1","darkorange3","darkorange4","deeppink1","deeppink3","darkseagreen1","darkseagreen4","darkslategrey","cyan1","cyan4",VsampleColor[c(3,5,7,9)],RedColor[c(2,4,6,9)])
mySampleColor =  c("mediumseagreen","deepskyblue4","firebrick","orange2")
#48-colors
colour40=c("#DC143C","#0000CD","#20B2AA","#FFA500","#9370DB",'#aaffc3',"#F08080",'#4363d8', '#f58231', 
           '#911eb4', '#46f0f0', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',  '#ffd8b1', '#000075',"#1E90FF","#7CFC00","#FFFF00",
           "#808000","#FF00FF","#7B68EE","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
           "#FF1493","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
           "#FF6347",'#e6194b', '#3cb44b', '#ffe119')
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#68-colors
col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))[-c(6,32)]
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF","#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))


library(ggsci)
sciColor_ls <- list()
sciColor_ls[["cols7_1"]] <- pal_locuszoom( alpha = 0.8)(7)
sciColor_ls[["cols7_2"]] <- pal_jama()(7)
sciColor_ls[["cols7_3"]] = pal_startrek(alpha = 0.8)(7)
sciColor_ls[["cols10_1"]] <- pal_npg("nrc", alpha = 0.9)(10)
sciColor_ls[["cols10_2"]] <- pal_aaas(alpha = 0.7)(10)
sciColor_ls[["cols10_3"]] <- pal_jco(alpha = 0.8)(10)
sciColor_ls[["cols8"]] <- pal_nejm(alpha = 0.8)(8)
sciColor_ls[["cols9"]] <- pal_lancet(alpha = 0.8)(9)
sciColor_ls[["cols26"]] = pal_ucscgb(alpha = 0.8)(26)
sciColor_ls[["cols20"]] = pal_d3("category20")(20)
sciColor_ls[["cols12_1"]] = pal_futurama(alpha = 0.8)(12)
sciColor_ls[["cols12_2"]] = pal_rickandmorty(alpha = 0.8)(12)
sciColor_ls[["cols16"]] = pal_simpsons(alpha = 0.9)(16)

sciColor_ls[["comb_col49"]] = unique(c(sciColor_ls[["cols16"]],sciColor_ls[["cols7_2"]],sciColor_ls[["cols26"]]))
sciColor_ls[["comb_col36"]] = unique(c(sciColor_ls[["cols7_1"]],sciColor_ls[["cols10_1"]][-c(2,6)],sciColor_ls[["cols10_3"]][1:5],sciColor_ls[["cols8"]][-c(1:2)],pal_futurama(alpha = 0.9)(11)[-6]))
sciColor_ls[["comb_col34"]] = unique(c(sciColor_ls[["cols9"]],sciColor_ls[["cols20"]],sciColor_ls[["cols10_3"]][6:10]))


#For cell type
#sciColor_ls[["cols10_1"]]
cols_types <- c("T"="#4DBBD5E5","NK"="#00A087E5","B"="#8491B4E5","PC"="#91D1C2E5","MC"="#367d7d",
                "pDC"="#688EC1","Myeloid"="#96CEE0", #cold colors
                "EC"="#D39200","FB"="#F39B7FE5","Epi"="#E64B35E5",
                "Stroma"="#7E6148E5","Nephron"="#F08080")
#show_col(cols_types)
cols_types_13 <- c("#2A8D8F","#74B980","#6EC4B3","#E9995E","#E9999F",
                   "#F1D99F", "#A6C256","#E6C14D","#8A9FE4","#5F76B8",
                   "#B986D2","#82539A","#B98A28")

###Discrete colors

col20 <- c("#5A8F64","#87AB8D","#8FBDC9","#C9C2DB","#A69BC1",
  "#7B2E09","#CF6B77","#DEA9BD","#EEE1C3","#C6B88B",
  "#A86421","#D59742","#97B0DD","#7696BA","#467CC6",
  "#9A476B","#8D7301","#E7C098","#B8C6A5","#86AC91")

##
warmCold_5 <- c("#983E31","#C66656","#88AFBF","#1B7B96","#0C546E")

nat_colors <- c("#C9553A", "#F2B134", "#5B9C6D", "#3C6B97", "#855FA8", "#D1AAE1", "#AF9F91", "#998A71", "#7B7678", "#3B3A36")
jacs_colors <- c("#FFD700", "#FFA500", "#FF8C00", "#FF6347", "#FF69B4", "#EE82EE", "#1E90FF", "#00CED1", "#32CD32", "#008080")
nc_colors <- c("#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3", "#03A9F4", "#00BCD4", "#009688", "#4CAF50")
ce_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#000000")
pnas_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
sci_colors <- c("#008EB0", "#95C623", "#E87E04", "#D91E36", "#2F2F2F", "#A3A3A3", "#EC008C", "#FBCB09", "#00AEB3", "#E40046")
nm_colors <- c("#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695", "#CF1C90", "#f97f27")
sci_adv_colors <- c("#1E90FF", "#FF4500", "#228B22", "#800080", "#00BFFF", "#FFD700", "#2E8B57", "#8B008B", "#FF6347", "#7B68EE")
#
nb_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")

cols13 <- c("#A9011B","#E4A826","#D2DCAD","#DCD6B2","#94BEA7","#4E7989","#75A0AE","#80944E","#547DB1","#9055A2","#D43B36","#F4C28F","#BAAFD1")
color_clone <- c("#FFD47F","#D9B9D4","#92A5D1","#C5DFF4","#C9DCC4","#DAA87C","#7C9895","#AEB2D1","#C25759","#D69D98","#599CB4","#92B5CA")



