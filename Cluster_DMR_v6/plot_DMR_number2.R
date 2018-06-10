##################################################################################################################
## Suffixes of all self-defined global variables must be "_g".
## Example:  
## Rscript  plot_DMR_number2.R         
          
outDir_g = "z-figures2"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }
##################################################################################################################


 
 



cat("\n\n\n\n\n Start part 1: \n\n")
##################################################################################################################
library(ggplot2) 


MyTheme_1_g <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
  theme(  
    line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                                                        ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                                                 ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all text elements.           "serif" for a serif font. 所有文本相关属性.
    title = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    ## aspect.ratio = 1,   ##高宽比
    
    axis.title    = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis label (element_text; inherits from axis.title)
    axis.title.y  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=90,      lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis label (element_text; inherits from axis.title)
    axis.text     = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1,  lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    
    axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                      ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	              ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                  ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                  ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.spacing      = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.spacing.x    = grid::unit(1, "mm", data=NULL) ,
    panel.spacing.y    = grid::unit(1, "mm", data=NULL) ,
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=0.5, vjust=0.5,   angle=NULL, lineheight=NULL),     ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL,    size=NULL, linetype=NULL, fill=NULL ), 	                                                      ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face="bold", colour=NULL, size=rel(1.2), hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 

MySaveGgplot2_1_g <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200,   device=cairo_ps)         
}

## df contains two columns, the first column (cond_col=1) is sample type, the second column (val_col=2) is value. (must be).
whisk_1_g <- function(df, cond_col=1, val_col=2) {  
  require(reshape2)
  condname <- names(df)[cond_col]  ## save the name of the first column.
  names(df)[cond_col] <- "cond" 
  names(df)[val_col]  <- "value"
  b   <- boxplot(value~cond, data=df, plot=FALSE)   
  df2 <- cbind(as.data.frame(b$stats), c("min","lq","m","uq","max"))
  names(df2) <- c(levels(df$cond), "pos")
  df2 <- melt(df2, id="pos", variable.name="cond")
  df2 <- dcast(df2, cond~pos)   
  names(df2)[1] <- condname 
  print(df2)
  df2
}



MyHistogram_1 <- function(x1, y1,  mytype1, mytype2,   mytype3,  path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4 ) {
  if( ! file.exists(path2)   ) { dir.create(path2, recursive = TRUE) }
  dataframeB  <- data.frame(list( xAxis = x1,  yAxis =y1,  Type1=mytype1, Type2=mytype2 , Type3=mytype3  )  ) 
  
  FigureTemp1 <- ggplot( data=dataframeB,  aes(x=xAxis, y=yAxis, fill=Type1 ) )  +  facet_wrap(Type3~Type2, scales = "free_x") +
    geom_bar(stat="identity",  position=position_dodge()  ) +  
    xlab(xLab2) + ylab("Number of DMRs") +  ggtitle(title2)  +  
    MyTheme_1_g( hjust1=1, vjust1=1,  angle1=60,   textSize=14 )   
  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_number",    sep="",  collapse=NULL),  
                    height1=height2,  width1=width2)  
}



MyLines_1 <- function(x1, y1,  mytype1, mytype2,  mytype3,  path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4 ) {
  if( ! file.exists(path2)   ) { dir.create(path2, recursive = TRUE) }
  dataframeB  <- data.frame(list( xAxis = x1,  yAxis =y1,  Type1=mytype1, Type2=mytype2 , Type3=mytype3 )  ) 
  
  FigureTemp1 <- ggplot( data=dataframeB,  aes(x=xAxis, y=yAxis, colour=Type1 , group=Type1 ) )  +  
    facet_wrap(Type3~Type2, scales = "free_x") +
    geom_line() +  geom_point() + 
    xlab(xLab2) + ylab("Number of DMRs") +  ggtitle(title2)  +  
    MyTheme_1_g( hjust1=1, vjust1=1,  angle1=60,   textSize=14 )   
  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_number",    sep="",  collapse=NULL),  
                    height1=height2,  width1=width2)  
}



MyLines_2 <- function(x1, y1,  mytype1, mytype2,   path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4 ) {
  if( ! file.exists(path2)   ) { dir.create(path2, recursive = TRUE) }
  dataframeB  <- data.frame(list( xAxis = x1,  yAxis =y1,  Type1=mytype1, Type2=mytype2  )  ) 
  
  FigureTemp1 <- ggplot( data=dataframeB, 
                         aes(x=xAxis, y=yAxis, 
                             linetype=interaction(Type2, Type1), color=interaction(Type1, Type2) , group=interaction(Type1, Type2) ) )  +  
    geom_line( ) + geom_point() + 
    xlab(xLab2) + ylab("Number of DMRs") +  ggtitle(title2)  +  
    MyTheme_1_g( hjust1=1, vjust1=1,  angle1=60,   textSize=14 )   
  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_number",    sep="",  collapse=NULL),  
                    height1=height2,  width1=width2)  
}



MyLines_3 <- function(x1, y1,  mytype1, mytype2,   mytype3, path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4 ) {
  if( ! file.exists(path2)   ) { dir.create(path2, recursive = TRUE) }
  dataframeB  <- data.frame(list( xAxis = x1,  yAxis =y1,  Type1=mytype1, Type2=mytype2 , Type3=mytype3  )  ) 
  
  FigureTemp1 <- ggplot( data=dataframeB, 
                         aes(x=xAxis, y=yAxis, 
                             group=interaction(Type1, Type2),  color=interaction(Type1, Type2) , 
                             linetype=interaction(Type1, Type2)    ) )  +  
    geom_line(size=0.7 ) + geom_point( size=2, alpha=1, aes(shape=interaction(Type1, Type2)  )   ) + 
    xlab(xLab2) + ylab("Number of DMRs") +  ggtitle(title2)  +  facet_wrap(~Type3, scales = "free_x", nrow=2) +
    scale_linetype_manual(values= rep(c("solid", "dotted"), times=length(levels(mytype2)))  )  +  
    scale_shape_manual(values = rep(c(19, 2), times=length(levels(mytype2))) )+
    scale_color_manual(values = rep( c("red3",  "cyan3",  "black"), each=length(levels(mytype1)))  ) +
    MyTheme_1_g( hjust1=1, vjust1=1,  angle1=60,   textSize=14 )   
  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_dotted",    sep="",  collapse=NULL),  
                    height1=height2,  width1=width2)  
  
  FigureTemp2 <- ggplot( data=dataframeB, 
                         aes(x=xAxis, y=yAxis, 
                             group=interaction(Type1, Type2),  color=interaction(Type1, Type2) , 
                             linetype=interaction(Type1, Type2)     ) )  +  
    geom_line(  ) + geom_point(aes(shape=interaction(Type1, Type2)),  size=2, alpha=1  ) + 
    xlab(xLab2) + ylab("Number of DMRs") +  ggtitle(title2)  +  facet_wrap(~Type3, scales = "free_x", nrow=2) +
    scale_linetype_manual(values= rep(c("solid", "dashed"), times=length(levels(mytype2)))  )  +  
    scale_shape_manual(values = rep(c(19, 2), times=length(levels(mytype2))) )+
    scale_color_manual(values = rep( c("red3",  "cyan3",  "black"), each=length(levels(mytype1)))  ) +
    MyTheme_1_g( hjust1=1, vjust1=1,  angle1=60,   textSize=14 )   
  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_dashed",    sep="",  collapse=NULL),  
                    height1=height2,  width1=width2)  
}






MyLines_4 <- function(x1, y1,  mytype1, mytype2,   mytype3, path2,   fileName2,  title2,  xLab2,  height2=4,  width2=4 ) {
  if( ! file.exists(path2)   ) { dir.create(path2, recursive = TRUE) }
  dataframeB  <- data.frame(list( xAxis = x1,  yAxis =y1,  Type1=mytype1, Type2=mytype2 , Type3=mytype3  )  ) 
  
  FigureTemp1 <- ggplot( data=dataframeB, 
                         aes(x=xAxis, y=yAxis, 
                             group=interaction(Type1, Type2),  color=interaction(Type1, Type2) , 
                             linetype=interaction(Type1, Type2)    ) )  +  
    geom_line(size=0.7 ) + geom_point( size=2, alpha=1, aes(shape=interaction(Type1, Type2)  )   ) + 
    xlab(xLab2) + ylab("Number of DMRs") +  ggtitle(title2)  +  facet_wrap(~Type3, scales = "free_x" ) +
    scale_linetype_manual(values= rep(c("solid", "dotted"), times=length(levels(mytype2)))  )  +  
    scale_shape_manual(values = rep(c(19, 2), times=length(levels(mytype2))) )+
    scale_color_manual(values = rep( c("red3",  "cyan3",  "black"), each=length(levels(mytype1)))  ) +
    MyTheme_1_g( hjust1=1, vjust1=1,  angle1=60,   textSize=14 )   
  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_dotted",    sep="",  collapse=NULL),  
                    height1=height2,  width1=width2)  
  
  FigureTemp2 <- ggplot( data=dataframeB, 
                         aes(x=xAxis, y=yAxis, 
                             group=interaction(Type1, Type2),  color=interaction(Type1, Type2) , 
                             linetype=interaction(Type1, Type2)     ) )  +  
    geom_line(  ) + geom_point(aes(shape=interaction(Type1, Type2)),  size=2, alpha=1  ) + 
    xlab(xLab2) + ylab("Number of DMRs") +  ggtitle(title2)  +  facet_wrap(~Type3, scales = "free_x" ) +
    scale_linetype_manual(values= rep(c("solid", "dashed"), times=length(levels(mytype2)))  )  +  
    scale_shape_manual(values = rep(c(19, 2), times=length(levels(mytype2))) )+
    scale_color_manual(values = rep( c("red3",  "cyan3",  "black"), each=length(levels(mytype1)))  ) +
    MyTheme_1_g( hjust1=1, vjust1=1,  angle1=60,   textSize=14 )   
  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_dashed",    sep="",  collapse=NULL),  
                    height1=height2,  width1=width2)  
}











###################
rawMatrix_ICSI <- read.table("ICSI-number.txt", header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_ICSI)
rawMatrix_ICSI2 = cbind(rawMatrix_ICSI, typeN=rep("ICSI", nrow(rawMatrix_ICSI) ) )



rawMatrix_IVF <- read.table("IVF-number.txt", header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_IVF)
rawMatrix_IVF2 = cbind(rawMatrix_IVF, typeN=rep("IVF", nrow(rawMatrix_IVF) ) )



mergeMatrix = rbind(rawMatrix_ICSI2, rawMatrix_IVF2)


MyHistogram_1(x1 = mergeMatrix$Comparisons, y1 = mergeMatrix$Number,  
              mytype1 = mergeMatrix$typeN, mytype2 = mergeMatrix$Analysis,  mytype3 = mergeMatrix$DMRs ,
              path2 =  outDir_g ,   fileName2 = "Histogram",  title2 = "",  
              xLab2 = "Comparisons",  height2=12,  width2=12 )



MyLines_1(x1 = mergeMatrix$Comparisons, y1 = mergeMatrix$Number,  
          mytype1 = mergeMatrix$typeN, mytype2 = mergeMatrix$Analysis,  mytype3 = mergeMatrix$DMRs ,
              path2 =  outDir_g ,   fileName2 = "Lines",  title2 = "",  
              xLab2 = "Comparisons",  height2=12,  width2=12 )




rawMatrix_2 = mergeMatrix
dim(rawMatrix_2)
rawMatrix_2[25:36, 1] = rawMatrix_2[1:12, 1]
rawMatrix_2[61:72, 1] = rawMatrix_2[1:12, 1]


MyLines_3(x1 = as.character( rawMatrix_2$Comparisons ), y1 = rawMatrix_2$Number,  
          mytype1 = mergeMatrix$typeN, mytype2 = mergeMatrix$Analysis,  mytype3 = mergeMatrix$DMRs ,
          path2 = outDir_g ,   fileName2 = "Lines2",  title2 = "",  
          xLab2 = "Comparisons",  height2=10,  width2=8 )




MyLines_4(x1 = as.character( rawMatrix_2$Comparisons ), y1 = rawMatrix_2$Number,  
          mytype1 = mergeMatrix$typeN, mytype2 = mergeMatrix$Analysis,  mytype3 = mergeMatrix$DMRs ,
          path2 = outDir_g ,   fileName2 = "Lines3",  title2 = "",  
          xLab2 = "Comparisons",  height2=6,  width2=12 )


























