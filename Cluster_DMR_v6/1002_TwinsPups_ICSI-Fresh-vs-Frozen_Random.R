##################################################################################################################
## Suffixes of all self-defined global variables must be "_g".
## Example:  
## Rscript  1002_TwinsPups_ICSI-Fresh-vs-Frozen_Random.R     12_Twins/A-rmXY     1002_TwinsPups_ICSI-Fresh-vs-Frozen_Random/A-rmXY    >     1002_TwinsPups_ICSI-Fresh-vs-Frozen_Random_A-rmXY.runLog  2>&1     
         
args_g <- commandArgs(TRUE)
print("##########################")
print("args: ")
print(args_g[1])   
print(args_g[2])     
print("##########################")

inputDir_g = args_g[1];     ## the path of input files
outDir_g   = args_g[2];     ## the path of output files

#   inputDir_g =  "12_Twins/C-onlyX"
#   outDir_g   =  "1002_TwinsPups_ICSI-Fresh-vs-Frozen_Random_test/C-onlyX"


print(inputDir_g)   
print(outDir_g)

if( ! file.exists(inputDir_g) ) { print("##### Error-1 #####") }
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }


lowestCoverage_g = 10   ## The lowest coverage of kept CpG sites in  CovSamples_g sampples of each group.
CovSamples_g = 5        ## The number of covered samples of each group for kept CpG sites.
CovSamples_g = as.integer(CovSamples_g)
##################################################################################################################


 


cat("\n\n\n\n\n Start part 1: \n\n")
##################################################################################################################
library(methylKit)
library(genomation)
library(ggplot2) 
library(ggfortify)
library(cluster)
library(lfda)
library(MASS)
library(factoextra)
library(magrittr)  
library(dplyr)  
library(rgl)
library(gdata)
library(ggrepel)
library(scatterplot3d)
library(car) 
library(plotly)
library(plot3D)
library(FactoMineR)
library(fpc)   
library(dendextend)
library(ggpubr)




myOutDir_sub1_g = paste(outDir_g, "/1_DefineVariables",  sep="") 
if( ! file.exists(myOutDir_sub1_g) ) { dir.create(myOutDir_sub1_g, recursive = TRUE) }


Files_one_1_g   = c()
Files_two_2_g   = c()
Files_three_3_g = c()
Files_four_4_g  = c()
Files_five_5_g  = c()

length( Files_one_1_g )
length( Files_two_2_g )
length( Files_three_3_g )
length( Files_four_4_g )
length( Files_five_5_g )

ShortName_one_1_g   = c()
ShortName_two_2_g   = c()
ShortName_three_3_g = c()
ShortName_four_4_g  = c()
ShortName_five_5_g  = c()

length( ShortName_one_1_g )
length( ShortName_two_2_g )
length( ShortName_three_3_g )
length( ShortName_four_4_g )
length( ShortName_five_5_g )



# Label of each group:
Label_one_1_g   = "random1_1"
Label_two_2_g   = "random1_2"
Label_three_3_g = "random2_1"
Label_four_4_g  = "random2_2"
Label_five_5_g  =  c()

Files_one_1_g1   = list.files( path = paste(inputDir_g, "/", "Pups_Fresh/ICSI-fresh_1", sep="") ,      pattern = "" ,    full.names = TRUE )
Files_two_2_g1   = list.files( path = paste(inputDir_g, "/", "Pups_Fresh/ICSI-fresh_2", sep="") ,      pattern = "" ,    full.names = TRUE )
Files_three_3_g1 = list.files( path = paste(inputDir_g, "/", "Pups_Frozen/ICSI-frozen_1", sep="") ,    pattern = "" ,    full.names = TRUE )
Files_four_4_g1  = list.files( path = paste(inputDir_g, "/", "Pups_Frozen/ICSI-frozen_2", sep="") ,    pattern = "" ,    full.names = TRUE )
# Files_five_5_g1  = list.files( path = paste(inputDir_g, "/", "Pups_Frozen/ICSI-frozen_2", sep="") ,    pattern = "" ,    full.names = TRUE )

ShortName_one_1_g1   = list.files( path = paste(inputDir_g, "/", "Pups_Fresh/ICSI-fresh_1", sep="") ,      pattern = "" ,    full.names = FALSE )
ShortName_two_2_g1   = list.files( path = paste(inputDir_g, "/", "Pups_Fresh/ICSI-fresh_2", sep="") ,      pattern = "" ,    full.names = FALSE )
ShortName_three_3_g1 = list.files( path = paste(inputDir_g, "/", "Pups_Frozen/ICSI-frozen_1", sep="") ,    pattern = "" ,    full.names = FALSE )
ShortName_four_4_g1  = list.files( path = paste(inputDir_g, "/", "Pups_Frozen/ICSI-frozen_2", sep="") ,    pattern = "" ,    full.names = FALSE )
# Files_five_5_g1  = list.files( path = paste(inputDir_g, "/", "Pups_Frozen/ICSI-frozen_2", sep="") ,    pattern = "" ,    full.names = FALSE )



myNumFiles1 = length(Files_one_1_g1)
myNumFiles2 = length(Files_three_3_g1)
myNumRandom1 = ceiling(myNumFiles1/2)
myNumRandom2 = ceiling(myNumFiles2/2)

myRandomIndex1A = sample( x=c(1:myNumFiles1), size=myNumRandom1 )
myRandomIndex2A = sample( x=c(1:myNumFiles2), size=myNumRandom2 )

Files_one_1_g   = c(  Files_three_3_g1[myRandomIndex2A], Files_two_2_g1[-myRandomIndex1A]  )
Files_two_2_g   = c(  Files_four_4_g1[myRandomIndex2A] , Files_one_1_g1[-myRandomIndex1A]  )    
Files_three_3_g = c( Files_two_2_g1[myRandomIndex1A],  Files_three_3_g1[-myRandomIndex2A]  )
Files_four_4_g  = c( Files_one_1_g1[myRandomIndex1A],  Files_four_4_g1[-myRandomIndex2A]   )

ShortName_one_1_g   = c(  ShortName_three_3_g1[myRandomIndex2A], ShortName_two_2_g1[-myRandomIndex1A]  )
ShortName_two_2_g   = c(  ShortName_four_4_g1[myRandomIndex2A] , ShortName_one_1_g1[-myRandomIndex1A]  )    
ShortName_three_3_g = c(  ShortName_two_2_g1[myRandomIndex1A],  ShortName_three_3_g1[-myRandomIndex2A]  )
ShortName_four_4_g  = c(  ShortName_one_1_g1[myRandomIndex1A],  ShortName_four_4_g1[-myRandomIndex2A]   )





sink( paste(myOutDir_sub1_g, "/1_All-Files.txt",  sep="")  )
cat("\n\n Files_one_1_g:\n")
Files_one_1_g 
cat("\n\n Files_two_2_g:\n")
Files_two_2_g 
cat("\n\n Files_three_3_g:\n")
Files_three_3_g 
cat("\n\n Files_four_4_g:\n")
Files_four_4_g 
cat("\n\n Files_five_5_g:\n")
Files_five_5_g 
sink()  


sink( paste(myOutDir_sub1_g, "/2_All-Files-shortName.txt",  sep="")  )
cat("\n\n ShortName_one_1_g:\n")
ShortName_one_1_g 
cat("\n\n ShortName_two_2_g:\n")
ShortName_two_2_g 
cat("\n\n ShortName_three_3_g:\n")
ShortName_three_3_g 
cat("\n\n ShortName_four_4_g:\n")
ShortName_four_4_g 
cat("\n\n ShortName_five_5_g:\n")
ShortName_five_5_g 
sink()  


sink( paste(myOutDir_sub1_g, "/3_Number-of-files-of-each-category.txt",  sep="")  )
length( Files_one_1_g )
length( Files_two_2_g  )
length( Files_three_3_g )
length( Files_four_4_g  )
length( Files_five_5_g  )
cat("\n\n\n")
length( ShortName_one_1_g )
length( ShortName_two_2_g  )
length( ShortName_three_3_g )
length( ShortName_four_4_g  )
length( ShortName_five_5_g  )
sink() 






Sex_one_1_g   = Files_one_1_g 
Sex_two_2_g   = Files_two_2_g
Sex_three_3_g = Files_three_3_g 
Sex_four_4_g  = Files_four_4_g
Sex_five_5_g  = Files_five_5_g

Type_one_1_g   = rep( x=Label_one_1_g,     times=length( Files_one_1_g ) ) 
Type_two_2_g   = rep( x=Label_two_2_g,     times=length( Files_two_2_g ) ) 
Type_three_3_g = rep( x=Label_three_3_g,   times=length( Files_three_3_g ) ) 
Type_four_4_g  = rep( x=Label_four_4_g,    times=length( Files_four_4_g ) ) 
Type_five_5_g  = rep( x=Label_five_5_g,    times=length( Files_five_5_g ) ) 


Treatment_one_1_g   = rep( x=1,  times=length( Files_one_1_g ) ) 
Treatment_two_2_g   = rep( x=2,  times=length( Files_two_2_g ) ) 
Treatment_three_3_g = rep( x=3,  times=length( Files_three_3_g ) ) 
Treatment_four_4_g  = rep( x=4,  times=length( Files_four_4_g ) ) 
Treatment_five_5_g  = rep( x=5,  times=length( Files_five_5_g ) ) 


SampleID_one_1_g   = Type_one_1_g 
SampleID_two_2_g   = Type_two_2_g
SampleID_three_3_g = Type_three_3_g 
SampleID_four_4_g  = Type_four_4_g
SampleID_five_5_g  = Type_five_5_g


for(i in c(1 : length(SampleID_one_1_g) ) ) {
  if( length(SampleID_one_1_g) )  {
    if(  grepl(pattern="boy",    x=Sex_one_1_g[i]) ) {Sex_one_1_g[i] = "boy"}
    if(  grepl(pattern="girl",   x=Sex_one_1_g[i]) ) {Sex_one_1_g[i] = "girl"}
    if(  grepl(pattern="father", x=Sex_one_1_g[i]) ) {Sex_one_1_g[i] = "father"}
    if(  grepl(pattern="mother", x=Sex_one_1_g[i]) ) {Sex_one_1_g[i] = "mother"}
    SampleID_one_1_g[i] = paste(SampleID_one_1_g[i], i, sep="_")
  }
}
for(i in c(1 : length(SampleID_two_2_g) ) ) {
  if( length(SampleID_two_2_g) )  {
    if(  grepl(pattern="boy",    x=Sex_two_2_g[i]) ) {Sex_two_2_g[i] = "boy"}
    if(  grepl(pattern="girl",   x=Sex_two_2_g[i]) ) {Sex_two_2_g[i] = "girl"}
    if(  grepl(pattern="father", x=Sex_two_2_g[i]) ) {Sex_two_2_g[i] = "father"}
    if(  grepl(pattern="mother", x=Sex_two_2_g[i]) ) {Sex_two_2_g[i] = "mother"}
    SampleID_two_2_g[i] = paste(SampleID_two_2_g[i], i, sep="_")
  }
}
for(i in c(1 : length(SampleID_three_3_g) ) ) {
  if( length(SampleID_three_3_g) )  {
    if(  grepl(pattern="boy",    x=Sex_three_3_g[i]) ) {Sex_three_3_g[i] = "boy"}
    if(  grepl(pattern="girl",   x=Sex_three_3_g[i]) ) {Sex_three_3_g[i] = "girl"}
    if(  grepl(pattern="father", x=Sex_three_3_g[i]) ) {Sex_three_3_g[i] = "father"}
    if(  grepl(pattern="mother", x=Sex_three_3_g[i]) ) {Sex_three_3_g[i] = "mother"}
    SampleID_three_3_g[i] = paste(SampleID_three_3_g[i], i, sep="_")
  }
}
for(i in c(1 : length(SampleID_four_4_g) ) ) {
  if( length(SampleID_four_4_g) >= 1 )  {
    if(  grepl(pattern="boy",    x=Sex_four_4_g[i]) ) {Sex_four_4_g[i] = "boy"}
    if(  grepl(pattern="girl",   x=Sex_four_4_g[i]) ) {Sex_four_4_g[i] = "girl"}
    if(  grepl(pattern="father", x=Sex_four_4_g[i]) ) {Sex_four_4_g[i] = "father"}
    if(  grepl(pattern="mother", x=Sex_four_4_g[i]) ) {Sex_four_4_g[i] = "mother"}
    SampleID_four_4_g[i] = paste(SampleID_four_4_g[i], i, sep="_")
  }
}
for(i in c(1 : length(SampleID_five_5_g) ) ) {
  if( length(SampleID_five_5_g) >= 1 )  {
    if(  grepl(pattern="boy",    x=Sex_five_5_g[i]) ) {Sex_five_5_g[i] = "boy"}
    if(  grepl(pattern="girl",   x=Sex_five_5_g[i]) ) {Sex_five_5_g[i] = "girl"}
    if(  grepl(pattern="father", x=Sex_five_5_g[i]) ) {Sex_five_5_g[i] = "father"}
    if(  grepl(pattern="mother", x=Sex_five_5_g[i]) ) {Sex_five_5_g[i] = "mother"}
    SampleID_five_5_g[i] = paste(SampleID_five_5_g[i], i, sep="_")
  }
}




ShortName_All_vector_g <- c(
  ShortName_one_1_g,
  ShortName_two_2_g,
  ShortName_three_3_g,
  ShortName_four_4_g,
  ShortName_five_5_g
)
ShortName_All_list_g <- as.list( ShortName_All_vector_g )


Files_All_vector_g <- c(
  Files_one_1_g,
  Files_two_2_g,
  Files_three_3_g,
  Files_four_4_g,
  Files_five_5_g
)
Files_All_list_g <- as.list( Files_All_vector_g )


Sex_All_vector_g <- c(
  Sex_one_1_g,
  Sex_two_2_g,
  Sex_three_3_g,
  Sex_four_4_g,
  Sex_five_5_g
)
Sex_All_list_g <- as.list( Sex_All_vector_g )


Type_All_vector_g <- c(
  Type_one_1_g,
  Type_two_2_g,
  Type_three_3_g,
  Type_four_4_g,
  Type_five_5_g
)
Type_All_list_g <- as.list( Type_All_vector_g )


Treatment_All_vector_g <- c(
  Treatment_one_1_g,
  Treatment_two_2_g,
  Treatment_three_3_g,
  Treatment_four_4_g,
  Treatment_five_5_g
)
Treatment_All_list_g <- as.list( Treatment_All_vector_g )



SampleID_All_vector_g <- c(
  SampleID_one_1_g,
  SampleID_two_2_g,
  SampleID_three_3_g,
  SampleID_four_4_g,
  SampleID_five_5_g
)
SampleID_All_list_g <- as.list( SampleID_All_vector_g )



if( is.null( Label_one_1_g   ) ) {Label_one_1_g   = "NA1"}
if( is.null( Label_two_2_g   ) ) {Label_two_2_g   = "NA2"}
if( is.null( Label_three_3_g ) ) {Label_three_3_g = "NA3"}
if( is.null( Label_four_4_g  ) ) {Label_four_4_g  = "NA4"}
if( is.null( Label_five_5_g  ) ) {Label_five_5_g  = "NA5"}
Label_one_1_g
Label_two_2_g
Label_three_3_g
Label_four_4_g
Label_five_5_g





## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=cyan3,   IVF-frozen=blue3,  ICSI-fresh=purple3, ICSI-frozen=red3
mySex_All_shape_g  = c( "boy"=16, "girl"=17, "father"=1, "mother"=2 ) 
myTech_All_color_g = c( "black",    "blue",    "cyan3" ,  "red" ,   "purple"   )
names(myTech_All_color_g) = c(Label_one_1_g, Label_two_2_g, Label_three_3_g, Label_four_4_g, Label_five_5_g)

cat("\n\n Shape and color: \n")
mySex_All_shape_g
myTech_All_color_g


cat("\n\n Length of variables: \n")
length( ShortName_All_vector_g )
length( ShortName_All_list_g )
length( Files_All_vector_g )
length( Files_All_list_g )
length( Sex_All_vector_g )
length( Sex_All_list_g )
length( Type_All_vector_g )
length( Type_All_list_g )
length( Treatment_All_vector_g )
length( Treatment_All_list_g )
length( SampleID_All_vector_g )
length( SampleID_All_list_g )


allSamples_info_g = cbind(ShortName_All_vector_g, Files_All_vector_g,  Sex_All_vector_g,  
                          Type_All_vector_g,  Treatment_All_vector_g,  SampleID_All_vector_g)

write.table(allSamples_info_g , 
            file = paste(myOutDir_sub1_g,   "4_allSamples_info.txt",  sep="/"), 
            append =FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



MySex_Shape_g <- function(  mySex_vector   ) {
  mySex_shape2  = mySex_vector
  for(i in c(1:length(mySex_shape2)) ) {
    if(mySex_shape2[i] == "boy")    { mySex_shape2[i] = c("boy"=16)   }
    if(mySex_shape2[i] == "girl")   { mySex_shape2[i] = c("girl"=17)  }
    if(mySex_shape2[i] == "father") { mySex_shape2[i] = c("father"=1) }
    if(mySex_shape2[i] == "mother") { mySex_shape2[i] = c("mother"=2) }
  }
  mySex_shape2 = as.numeric(mySex_shape2)
  names(mySex_shape2) = mySex_vector
  return(mySex_shape2)
}


MyTech_color_g <- function(  myTech_vector  ) {
  myTech_color  = myTech_vector
  for(i in c(1:length(myTech_color)) ) {
    if(myTech_color[i] == Label_one_1_g)   { myTech_color[i] = "black"  }
    if(myTech_color[i] == Label_two_2_g)   { myTech_color[i] = "blue"   }
    if(myTech_color[i] == Label_three_3_g) { myTech_color[i] = "cyan3"  }
    if(myTech_color[i] == Label_four_4_g)  { myTech_color[i] = "red"    }
    if(myTech_color[i] == Label_five_5_g)  { myTech_color[i] = "purple" }
  }
  names(myTech_color) = myTech_vector
  return(myTech_color)
}


MyTech_number_g <- function(  myTech_vector  ) {
  myTech_color  = myTech_vector
  for(i in c(1:length(myTech_color)) ) {
    if(myTech_color[i] == Label_one_1_g)   { myTech_color[i] = 1 }
    if(myTech_color[i] == Label_two_2_g)   { myTech_color[i] = 2 }
    if(myTech_color[i] == Label_three_3_g) { myTech_color[i] = 3 }
    if(myTech_color[i] == Label_four_4_g)  { myTech_color[i] = 4 }
    if(myTech_color[i] == Label_five_5_g)  { myTech_color[i] = 5 }
  }
  return(myTech_color)
}




continue_on_error <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())' " )
}
# This option is very important.
options(error=continue_on_error) 
##################################################################################################################
cat("\n\n End part 1: \n\n\n\n\n")





cat("\n\n\n\n\n Start part 2: \n\n")
##################################################################################################################
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
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
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





MyCluster_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) { 
  if( nrow(mymeth1)*(1-sdThres1) > 10 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",         sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",         sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  dev.off()
  }
}

MyCluster_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {  
  myTempValue = tryCatch(
    clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE, plot=FALSE  ), 
    error = function(err){"000"}
  )
  if( length( as.character(myTempValue) ) > 1) {
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",         sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",         sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  dev.off()
  }
}

MyCluster_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
  myTempFunction <- function() {
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99)                    
  
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    )                     
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 )   
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  )
}





MyPCA_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) { 
  if( nrow(mymeth1)*(1-sdThres1) > 10 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  dev.off()
  pdf( file=paste(path1, "/", file1, ".screeplot.pdf", sep="") , width=width1/2, height=height1/2  )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  dev.off()
  }
}

MyPCA_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {
  myTempFunction <- function() {  
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  dev.off()
  pdf( file=paste(path1, "/", file1, ".screeplot.pdf", sep="") , width=width1, height=height1  )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  dev.off()
  }

  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  )
}

MyPCA_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
  myTempFunction <- function() {
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99)                    

  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    )                     
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 )                    
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  ) 
}





#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
MyPrcompObj_1_g <- function(  prcompObj2,   path2,   file2,  dataFrame_temp2   ) {
  
  sink( file = paste(path2, "/1A_PCA_results_",  file2,  ".txt",  sep="") )
  print( prcompObj2 )
  print( names(prcompObj2) )
  sink()
  
  sink( file = paste(path2, "/1B_PCA_summary_",  file2,  ".txt",   sep="") )
  print( summary(prcompObj2) )
  sink()
  
  sink( file = paste(path2, "/2_PCA_all_",  file2,  ".txt",   sep="") )
  print("####################### prcompObj2$sdev #########################")
  print(prcompObj2$sdev)
  print("####################### prcompObj2$rotation #########################")
  print(prcompObj2$rotation)
  print("####################### prcompObj2$center #########################")
  print(prcompObj2$center)
  print("####################### prcompObj2$scale #########################")
  print(prcompObj2$scale)
  print("####################### prcompObj2$x #########################")
  print(prcompObj2$x)
  sink()
  
  pdf( file = paste(path2, "/3_PCA_info_",  file2,  ".pdf",   sep="")  )
  plot(prcompObj2, type="lines")
  print( fviz_eig(prcompObj2) )
  dev.off() 
  
  my_fviz_pca_ind1 <- fviz_pca_ind(prcompObj2,
                                   col.ind = "cos2", # Color by the quality of representation
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                   repel = TRUE     # Avoid text overlapping
  )
  my_fviz_pca_ind2 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE
  )
  my_fviz_pca_ind3 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none"
  )
  my_fviz_pca_ind4 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   ggtheme = MyTheme_1_g()
  )
  my_fviz_pca_ind5 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none",
                                   ggtheme = MyTheme_1_g()
  )
  
  pdf( file=paste(path2, "/4_PCA-2D_", file2, ".pdf",  sep="") )
  print(my_fviz_pca_ind1)
  print(my_fviz_pca_ind2)
  print(my_fviz_pca_ind3)
  print(my_fviz_pca_ind4)
  print(my_fviz_pca_ind5)
  dev.off() 
  
  #############################
  prcompObj2_matrix <- prcompObj2$x
  prcompObj2_Contri  <- (prcompObj2$sdev)^2
  prcompObj2_Contri  <- prcompObj2_Contri/sum(prcompObj2_Contri)
  prcompObj2_Contri  <- prcompObj2_Contri * 100
  prcompObj2_Contri  <- round(prcompObj2_Contri, 2)
  
  label1_2two <-   paste( "PC1 ",  "(", prcompObj2_Contri[1], "%)", sep="" )
  label2_2two <-   paste( "PC2 ",  "(", prcompObj2_Contri[2], "%)", sep="" )
  label3_2two <-   paste( "PC3 ",  "(", prcompObj2_Contri[3], "%)", sep="" ) 
  
  dataframeA_2two  <- data.frame( as.data.frame(prcompObj2_matrix), mySex= as.vector(dataFrame_temp2$mysex), 
                                  myTech=as.vector(dataFrame_temp2$mytech),    myLabel=as.vector(dataFrame_temp2$mysampleID)   ) 
  
  FigureTemp1_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=3, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1_2two ,  path1=path2, fileName1="5A_PCA-PC1-PC2",  height1=3.5,  width1=6)
  
  FigureTemp2_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=5, alpha=0.5  )+ xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2_2two ,  path1=path2, fileName1="5B_PCA-PC1-PC2-alpha",   height1=3.5,  width1=6)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot",   height1=3.5,  width1=6)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot-alpha",   height1=3.5,  width1=6)
  
  FigureTemp4 <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=6, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +     
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4_2two ,  path1=path2, fileName1="5D_PCA-PC1-PC2-big",   height1=3.5,  width1=6)
  
  FigureTemp5_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=2, alpha=0.7  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp5_2two ,  path1=path2, fileName1="5E_PCA-PC1-PC2-Labels",   height1=3.5,  width1=6)
  
  FigureTemp6_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=5, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp6_2two ,  path1=path2, fileName1="5F_PCA-PC1-PC2-Labels",   height1=3.5,  width1=6)
  
  
  ## PC1, 2, and 3     dev.off() 
  pdf( file = paste(path2, "6_PCA-3d-by-scatter3D.pdf",  sep="/") )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),     
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),          
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ##############
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ######
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  dev.off()
}

MyPCA_1A_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1,  dataFrame_temp1    ) {   
  if( nrow(mymeth1)*(1-sdThres1) > 10 ) {   
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  #myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  #myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, file1, "filterByQuantile_1", sep="/")
  path1_2 = paste(path1, file1, "filterByQuantile_2", sep="/")
  path1_3 = paste(path1, file1, "filterByQuantile_3", sep="/")
  path1_4 = paste(path1, file1, "filterByQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  #if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  #if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  #MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  #MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  }
}

MyPCA_2A_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1,  dataFrame_temp1   ) { 
  myTempFunction <- function() {
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }   
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  #myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  #myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, file1, "filterNoQuantile_1", sep="/")
  path1_2 = paste(path1, file1, "filterNoQuantile_2", sep="/")
  path1_3 = paste(path1, file1, "filterNoQuantile_3", sep="/")
  path1_4 = paste(path1, file1, "filterNoQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  #MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  #MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )  
  }
  
  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  )
}

MyPCA_3A_g <- function(  mymeth2 ,  path2,   file2, width2, height2,  dataFrame_temp2   )  {
  myTempFunction <- function() {
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99 , dataFrame_temp1=dataFrame_temp2   )                     

  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 ,  dataFrame_temp1=dataFrame_temp2   )                     
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"000111"}
  ) 
}





My_for_MultidimensionalScaling_1_g <- function(  path3,   dataFrame_temp3 ) {
  FigureTemp1_2two  <- ggplot( data = dataFrame_temp3, aes(x = myX, y = myY, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=3, alpha=1  ) + xlab("Dimension 1") +   ylab("Dimension 2") +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1_2two ,  path1=path3, fileName1="A_MDS",  height1=3.5,  width1=6)
  
  FigureTemp2_2two  <- ggplot( data = dataFrame_temp3, aes(x = myX, y = myY, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=5, alpha=0.5  )+ xlab("Dimension 1") +   ylab("Dimension 2") +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2_2two ,  path1=path3, fileName1="B_MDS-alpha",   height1=3.5,  width1=6)
  
  FigureTemp3_2two  <- ggplot( data = dataFrame_temp3, aes(x = myX, y = myY, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=1  ) + xlab("Dimension 1") +   ylab("Dimension 2") +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path3, fileName1="C_MDS-smallDot",   height1=3.5,  width1=6)
  
  FigureTemp3_2two  <- ggplot( data = dataFrame_temp3, aes(x = myX, y = myY, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=0.5  ) + xlab("Dimension 1") +   ylab("Dimension 2") +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path3, fileName1="D_MDS-smallDot-alpha",   height1=3.5,  width1=6)
  
  FigureTemp4 <- ggplot( data = dataFrame_temp3, aes(x = myX, y = myY, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=6, alpha=0.5  ) + xlab("Dimension 1") +   ylab("Dimension 2") +     
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4_2two ,  path1=path3, fileName1="E_MDS-big",   height1=3.5,  width1=6)
  
  FigureTemp5_2two  <- ggplot( data = dataFrame_temp3, aes(x = myX, y = myY, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=2, alpha=0.7  ) + xlab("Dimension 1") +   ylab("Dimension 2") +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp5_2two ,  path1=path3, fileName1="F_MDS-Labels",   height1=3.5,  width1=6)
  
  FigureTemp6_2two  <- ggplot( data = dataFrame_temp3, aes(x = myX, y = myY, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=5, alpha=0.5  ) + xlab("Dimension 1") +   ylab("Dimension 2") +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp6_2two ,  path1=path3, fileName1="G_MDS-Labels",   height1=3.5,  width1=6)
  
}

MyMultidimensionalScaling_1_g <- function(  meLevelMatrix2,   path2,   dataFrame_temp2 ) {
  if( ! file.exists(path2) ) { dir.create(path2, recursive = TRUE) }
  path2_temp1 = paste(path2, "/Classical-Metric-Multidimensional-Scaling", sep="")
  if( ! file.exists(path2_temp1) ) { dir.create(path2_temp1, recursive = TRUE) }
  res.dist_temp1 <- get_dist( t(meLevelMatrix2) ,   method = "euclidean")
  mds_temp1 = cmdscale(res.dist_temp1,   k = 2)
  dataframeA_temp1  <- data.frame( myX=mds_temp1[,1],  myY=mds_temp1[,2],  myLabel=rownames(mds_temp1),
                                   mySex= as.vector(dataFrame_temp2$mysex),  myTech=as.vector(dataFrame_temp2$mytech)  ) 
  My_for_MultidimensionalScaling_1_g(  path3=path2_temp1,   dataFrame_temp3=dataframeA_temp1 )
  
  path2_temp2 = paste(path2, "/Nonmetric-Multidimensional-Scaling", sep="")
  if( ! file.exists(path2_temp2) ) { dir.create(path2_temp2, recursive = TRUE) }
  res.dist_temp2 <- res.dist_temp1
  mds_temp2 = isoMDS(res.dist_temp2,   k = 2)
  dataframeA_temp2  <- data.frame( myX=mds_temp2$points[,1],  myY=mds_temp2$points[,2],  myLabel=rownames(mds_temp2$points),
                                   mySex= as.vector(dataFrame_temp2$mysex),  myTech=as.vector(dataFrame_temp2$mytech)  ) 
  My_for_MultidimensionalScaling_1_g(  path3=path2_temp2,   dataFrame_temp3=dataframeA_temp2 )
  
  cat("\n\n\n\n\nMDS:")
  cat(mds_temp1)
  cat("\n\n")
  print(mds_temp2)
  cat("\n\n\n\n\n")
  
  sink( file=paste(path2, "MDS_info.txt", sep="/") )
  cat("\n\n\n\n\nMDS:")
  cat(mds_temp1)
  cat("\n\n\n\n\n\n\n")
  print(mds_temp2)
  cat("\n\n\n\n\n")
  sink() 
}

myDMRs_annotation_g <- function(  myDiffDMR_5,   path2_5  ) {
  if( ! file.exists(path2_5) ) { dir.create(path2_5, recursive = TRUE) }
  
  ########################## all genes.
  Ann_gene = annotateWithGeneParts( as(myDiffDMR_5,  "GRanges"),  gene.obj_2000bp_g)
  sink( file=paste(path2_5, "1000A-distribution-onGenes.txt", sep="/")   )
  getFeatsWithTargetsStats( Ann_gene,  percentage=TRUE)
  print(Ann_gene)
  sink()
  pdf( file=paste(path2_5, "1000B-distribution-onGenes.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_gene,  precedence=TRUE, main="Annotation on all genes") )
  dev.off()
  

  ########################## ImprintedRegions 1
  diffImprinted1Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint1.obj_g$ImprintedRegions,  flank=imprint1.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002A-distribution-onImprintedRegions1.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted1Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002B-distribution-onImprintedRegions1.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted1Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 2
  diffImprinted2Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint2.obj_g$ImprintedRegions,  flank=imprint2.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002C-distribution-onImprintedRegions2.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted2Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002D-distribution-onImprintedRegions2.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted2Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 3
  diffImprinted3Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint3.obj_g$ImprintedRegions,  flank=imprint3.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002E-distribution-onImprintedRegions3.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted3Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002F-distribution-onImprintedRegions3.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted3Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 4
  diffImprinted4Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint4.obj_g$ImprintedRegions,  flank=imprint4.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002G-distribution-onImprintedRegions4.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted4Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002H-distribution-onImprintedRegions4.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted4Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 5
  diffImprinted5Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint5.obj_g$ImprintedRegions,  flank=imprint5.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002I-distribution-onImprintedRegions5.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted5Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002J-distribution-onImprintedRegions5.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted5Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  
  
  ########################## H3K4me1 peaks
  Ann_H3K4me1 = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                         feature = region_1_H3K4me1.obj_g$H3K4me1,  flank = region_1_H3K4me1.obj_g$shores,
                                         feature.name = "H3K4me1",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "1A-distribution-on-H3K4me1.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K4me1,  percentage=TRUE)
  print(Ann_H3K4me1)
  sink()
  pdf( file=paste(path2_5, "1B-distribution-on-H3K4me1.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K4me1, precedence=TRUE, main="Annotation on H3K4me1 peaks"))
  dev.off()
  
  
  ########################## H3K4me3 peaks
  Ann_H3K4me3 = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                         feature = region_2_H3K4me3.obj_g$H3K4me3,  flank = region_2_H3K4me3.obj_g$shores,
                                         feature.name = "H3K4me3",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "2A-distribution-on-H3K4me3.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K4me3,  percentage=TRUE)
  print(Ann_H3K4me3)
  sink()
  pdf( file=paste(path2_5, "2B-distribution-on-H3K4me3.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K4me3, precedence=TRUE, main="Annotation on H3K4me3 peaks"))
  dev.off()
  
  
  
  ########################## H3K27ac peaks
  Ann_H3K27ac = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                         feature = region_3_H3K27ac.obj_g$H3K27ac,  flank = region_3_H3K27ac.obj_g$shores,
                                         feature.name = "H3K27ac",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "3A-distribution-on-H3K27ac.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K27ac,  percentage=TRUE)
  print(Ann_H3K27ac)
  sink()
  pdf( file=paste(path2_5, "3B-distribution-on-H3K27ac.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K27ac, precedence=TRUE, main="Annotation on H3K27ac peaks"))
  dev.off()
  
  
  ########################## H3K27me3 peaks
  Ann_H3K27me3 = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                          feature = region_4_H3K27me3.obj_g$H3K27me3,  flank = region_4_H3K27me3.obj_g$shores,
                                          feature.name = "H3K27me3",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "4A-distribution-on-H3K27me3.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K27me3,  percentage=TRUE)
  print(Ann_H3K27me3)
  sink()
  pdf( file=paste(path2_5, "4B-distribution-on-H3K27me3.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K27me3, precedence=TRUE, main="Annotation on H3K27me3 peaks"))
  dev.off()
  
  
  
  ########################## 380imprint regions
  Ann_380imprint = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_5_380imprint.obj_g$imprint,  flank = region_5_380imprint.obj_g$shores,
                                            feature.name = "imprint",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "5A-distribution-on-380imprint.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_380imprint,  percentage=TRUE)
  print(Ann_380imprint)
  sink()
  pdf( file=paste(path2_5, "5B-distribution-on-380imprint.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_380imprint, precedence=TRUE, main="Annotation on 380imprint regions"))
  dev.off()
  
  
  
  ########################## centromere regions
  Ann_centromere = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_6_centromeres.obj_g$centromere,  flank = region_6_centromeres.obj_g$shores,
                                            feature.name = "centromere",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "6A-distribution-on-centromere.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_centromere,  percentage=TRUE)
  print(Ann_centromere)
  sink()
  pdf( file=paste(path2_5, "6B-distribution-on-centromere.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_centromere, precedence=TRUE, main="Annotation on centromere regions"))
  dev.off()
  
  
  
  ########################## CpGislands regions
  Ann_CpGislands = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_7_CpGislands.obj_g$CpGislands,  flank = region_7_CpGislands.obj_g$shores,
                                            feature.name = "CpGislands",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "7A-distribution-on-CpGislands.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_CpGislands,  percentage=TRUE)
  print(Ann_CpGislands)
  sink()
  pdf( file=paste(path2_5, "7B-distribution-on-CpGislands.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_CpGislands, precedence=TRUE, main="Annotation on CpGislands regions"))
  dev.off()
  
 
  
  ########################## HouseKeeping genes
  Ann_HouseKeeping = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                              feature = region_9_HouseKeeping.obj_g$HouseKeeping,  flank = region_9_HouseKeeping.obj_g$shores,
                                              feature.name = "HouseKeeping",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "9A-distribution-on-HouseKeeping.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_HouseKeeping,  percentage=TRUE)
  print(Ann_HouseKeeping)
  sink()
  pdf( file=paste(path2_5, "9B-distribution-on-HouseKeeping.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_HouseKeeping, precedence=TRUE, main="Annotation on HouseKeeping genes"))
  dev.off()
  
  
  
  ########################## rRNA_Genes genes
  Ann_rRNA_Genes = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_10_rRNA.obj_g$rRNA_Genes,  flank = region_10_rRNA.obj_g$shores,
                                            feature.name = "rRNA_Genes",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "10A-distribution-on-rRNA_Genes.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_rRNA_Genes,  percentage=TRUE)
  print(Ann_rRNA_Genes)
  sink()
  pdf( file=paste(path2_5, "10B-distribution-on-rRNA_Genes.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_rRNA_Genes, precedence=TRUE, main="Annotation on rRNA_Genes genes"))
  dev.off()
  


  ########################## activeEnhancers 
  Ann_activeEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_11A_activeEn.obj_g$activeEnhancers,  flank = region_11A_activeEn.obj_g$shores,
                                                 feature.name = "activeEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11A-distribution-on-activeEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_activeEnhancers,  percentage=TRUE)
  print(Ann_activeEnhancers)
  sink()
  pdf( file=paste(path2_5, "11A-distribution-on-activeEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_activeEnhancers, precedence=TRUE, main="Annotation on activeEnhancers"))
  dev.off()
  
  
  
  ########################## otherEnhancers 
  Ann_otherEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                feature = region_11B_otherEn.obj_g$otherEnhancers,  flank = region_11B_otherEn.obj_g$shores,
                                                feature.name = "otherEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11B-distribution-on-otherEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_otherEnhancers,  percentage=TRUE)
  print(Ann_otherEnhancers)
  sink()
  pdf( file=paste(path2_5, "11B-distribution-on-otherEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_otherEnhancers, precedence=TRUE, main="Annotation on otherEnhancers"))
  dev.off()
  
  ########################## posiedEnhancers 
  Ann_posiedEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_11C_posiedEn.obj_g$posiedEnhancers,  flank = region_11C_posiedEn.obj_g$shores,
                                                 feature.name = "posiedEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11C-distribution-on-posiedEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_posiedEnhancers,  percentage=TRUE)
  print(Ann_posiedEnhancers)
  sink()
  pdf( file=paste(path2_5, "11C-distribution-on-posiedEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_posiedEnhancers, precedence=TRUE, main="Annotation on posiedEnhancers"))
  dev.off()
  
  
  
  ########################## primedEnhancers 
  Ann_primedEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_11D_primedEn.obj_g$primedEnhancers,  flank = region_11D_primedEn.obj_g$shores,
                                                 feature.name = "primedEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11D-distribution-on-primedEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_primedEnhancers,  percentage=TRUE)
  print(Ann_primedEnhancers)
  sink()
  pdf( file=paste(path2_5, "11D-distribution-on-primedEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_primedEnhancers, precedence=TRUE, main="Annotation on primedEnhancers"))
  dev.off()
  
 
}  

myDiff_DMC_DMR_g <- function(  methobj2,   path2  ) {
  sink(paste(path2,  "0_methobj2-for-diffMe.txt",  sep="/"))
  print(methobj2)
  sink() 
  
  ## Depending on the sample size per each set it will either use Fisher’s exact or logistic 
  ##  regression to calculate P-values.
  ##  If you have replicates, the function will automatically use logistic regression.
  print("################## start calculateDiffMeth:")
  myDiff_2two_sub2 = calculateDiffMeth(methobj2,  mc.cores=8 )   
  print("################## End calculateDiffMeth") 
  dim(myDiff_2two_sub2)
  names(myDiff_2two_sub2)
  head(myDiff_2two_sub2)
  
  
  
  myOutDir_sub2_g3 = paste(path2, "/Segmentation",  sep="") 
  if( ! file.exists(myOutDir_sub2_g3) ) { dir.create(myOutDir_sub2_g3, recursive = TRUE) }
  
  pdf( file=paste(path2, "Segmentation.pdf", sep="/")  )
  res_temp_1 = methSeg( myDiff_2two_sub2 , diagnostic.plot=TRUE)
  methSeg2bed(res_temp_1,   filename= paste(myOutDir_sub2_g3, "Segmentation.bed", sep="/") )
  dev.off()
  
  
  
  myQvalue_2two_sub2 = myDiff_2two_sub2$qvalue
  myMethDi_2two_sub2 = myDiff_2two_sub2$meth.diff
  
  pdf(paste(path2,  "1A_qvalue_distribution.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=FALSE) )   
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=TRUE ) )   
  dev.off() 

  myQvalue_2two_sub2_log10 <- -log10(myQvalue_2two_sub2)
  pdf(paste(path2,  "1B_qvalue_distribution_log10.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=TRUE ) )
  dev.off() 
  
  pdf(paste(path2,  "1C_methDiff_distribution.pdf",  sep="/"))
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=TRUE ) )
  dev.off() 
             
  qvalue_cutoff = 0.05
  methDiff_cutoff = 5
  
  myColor1_2two_sub2 <- rep( "no",   times= length(myQvalue_2two_sub2) )
  myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
  number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.01
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.05
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.1
    methDiff_cutoff = 5
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.5
    methDiff_cutoff = 1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  
  print("##############################")
  print("##############################")
  print("The final parameters:")
  print(number_yes)
  print(qvalue_cutoff)
  print(methDiff_cutoff)
  print("##############################")
  print("##############################")

  
  DataFrame2_2two_sub2 <- data.frame(myx1 = myMethDi_2two_sub2,   
                                     myy1 =  -log10(myQvalue_2two_sub2),  
                                     mycolor1 = myColor1_2two_sub2 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + ylim(0, 200)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  + ylim(0, 100)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)   + ylim(0, 30)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )

  
  
  myDiff25p.hypo_2two_sub2  = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hypo" )  ## less enrich in ART
  myDiff25p.hyper_2two_sub2 = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hyper")  ## more enrich in ART
  myDiff25p_2two_sub2       = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff)
  myDiffTemp_2two_sub2      = getMethylDiff(myDiff_2two_sub2, difference=0,  qvalue=0.05)
  
  write.table(myDiff_2two_sub2 , file = paste(path2,"2A_diffMe-allsites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hypo_2two_sub2 , file = paste(path2,"2B_diffMe-hypo.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hyper_2two_sub2 , file = paste(path2,"2C_diffMe-hyper.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p_2two_sub2 , file = paste(path2,"2D_AlldiffMesites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiffTemp_2two_sub2 , file = paste(path2,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  sink( file=paste(path2, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=FALSE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  sink()
  
  pdf( file=paste(path2, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=TRUE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  dev.off()
  

  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hypo_2two_sub2,    path2_5 = paste(path2, "Annotation_Hypo",  sep="/") ),
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hyper_2two_sub2,   path2_5 = paste(path2, "Annotation_Hyper", sep="/") ), 
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p_2two_sub2,         path2_5 = paste(path2, "Annotation_All",   sep="/") ),
        error = function(err){"000"}
  )
 
}





#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myMainFunction_1_g  <- function(  myobj_temp1,   path_temp1,   binSize_temp1, binBases_temp1, dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }

  meth_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  write.table(dataFrame_temp1 , 
              file = paste(path_temp1_sub1,   "0_dataFrame.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(meth_2two) )
  print( getSampleID(meth_2two) )
  print( getTreatment(meth_2two) )
  sink()
  
  

  
  meth_2two = tileMethylCounts( meth_2two,   win.size=binSize_temp1,   step.size=binSize_temp1,   cov.bases = binBases_temp1  )   
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2_dimensions-tiles-merged.txt", sep="/")  )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
  
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(getData(meth_2two)[,1:4] , 
              file = paste(path_temp1_sub1,   "3C_regions-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  
  meth_sub4_matrix_g <- getData(meth_2two)
  mySamples_firstRow = colnames( mat_2two )
  myCoverageMatrix = mat_2two 
  myMeLevelMatrix  = mat_2two 
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myCoverageMatrix[,i] = meth_sub4_matrix_g[,3*i+2]
    myMeLevelMatrix[,i]  = 100*meth_sub4_matrix_g[,3*i+3]/meth_sub4_matrix_g[,3*i+2] 
  }
  
  
  pdf( file=paste(path_temp1_sub1, "4_Coverage_MeLevel_Cor.pdf", sep="/")  )
  par(mfrow=c(3,1))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    ## i =1 
    myTempVec1 = myCoverageMatrix[,i]
    myTempBool1 = (myTempVec1<=70) 
    myTempDataframe = data.frame(x1=myTempVec1[myTempBool1], y1=myMeLevelMatrix[myTempBool1,i] )
    boxplot(y1 ~ x1, data = myTempDataframe , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i] )  
    
    myTempBool2 = ( (myTempVec1 >= 70) & (myTempVec1 <= 140)  )
    myTempDataframe2 = data.frame(x2=myTempVec1[myTempBool2], y2=myMeLevelMatrix[myTempBool2,i] )
    boxplot(y2 ~ x2, data = myTempDataframe2 , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i]  ) 
    
    myTempBool2 =  (myTempVec1 >= 140)  
    myTempDataframe2 = data.frame(x2=myTempVec1[myTempBool2], y2=myMeLevelMatrix[myTempBool2,i] )
    boxplot(y2 ~ x2, data = myTempDataframe2 , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i]  ) 
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4A_Coverage.eachSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  for( i in c(1:length(as.character(dataFrame_temp1$myfiles))) ) {
    myTempCov1 = myCoverageMatrix[,i]
    file_name = as.character(dataFrame_temp1$myfiles)[i]
    cat(   min(myTempCov1, na.rm = TRUE), "\t", mean(myTempCov1, na.rm = TRUE), "\t", 
           median(myTempCov1, na.rm = TRUE), "\t",  max(myTempCov1, na.rm = TRUE), "\t", file_name, "\n"  )    
  }
  sink()
  
  sink( file=paste(path_temp1_sub1, "4B_MeLevel.eachSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  for( i in c(1:length(as.character(dataFrame_temp1$myfiles))) ) {
    myTempCov1 = myMeLevelMatrix[,i]
    file_name = as.character(dataFrame_temp1$myfiles)[i]
    cat(   min(myTempCov1, na.rm = TRUE), "\t", mean(myTempCov1, na.rm = TRUE), "\t",  
           median(myTempCov1, na.rm = TRUE), "\t",  max(myTempCov1, na.rm = TRUE), "\t", file_name, "\n"  )    
  }
  sink()
  
  
  myCoverageMatrix_vector = as.vector( myCoverageMatrix ) 
  myMeLevelMatrix_vector  = as.vector( myMeLevelMatrix ) 
  
  
  pdf( file=paste(path_temp1_sub1, "5_Coverage_MeLevel_Cor_poolAllSamples.pdf", sep="/")  )
  par(mfrow=c(3,1))
  myTempBool3 =  (myCoverageMatrix_vector<=50) 
  myTempDataframe3 = data.frame(x3=myCoverageMatrix_vector[myTempBool3], y3=myMeLevelMatrix_vector[myTempBool3])
  boxplot(y3 ~ x3, data = myTempDataframe3  , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples") 
  myTempBool4 =  ( (myCoverageMatrix_vector>50)  & (myCoverageMatrix_vector<=100) )
  myTempDataframe4 = data.frame(x4=myCoverageMatrix_vector[myTempBool4], y4=myMeLevelMatrix_vector[myTempBool4] )
  boxplot(y4 ~ x4, data = myTempDataframe4  , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples") 
  myTempBool5 =  ( (myCoverageMatrix_vector>100)  & (myCoverageMatrix_vector<=150) )
  myTempFunction100 <- function() {
    myTempDataframe5 = data.frame(x5=myCoverageMatrix_vector[myTempBool5], y5=myMeLevelMatrix_vector[myTempBool5] )
    boxplot(y5 ~ x5, data = myTempDataframe5 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
    myTempBool6 =  ( (myCoverageMatrix_vector>150)  & (myCoverageMatrix_vector<=200) )
    myTempDataframe6 = data.frame(x6=myCoverageMatrix_vector[myTempBool6], y6=myMeLevelMatrix_vector[myTempBool6] )
    boxplot(y6 ~ x6, data = myTempDataframe6 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
    myTempBool7 = (myCoverageMatrix_vector>200) 
    myTempDataframe7 = data.frame(x7=myCoverageMatrix_vector[myTempBool7], y7=myMeLevelMatrix_vector[myTempBool7] )
    boxplot(y7 ~ x7, data = myTempDataframe7 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
  }
  tryCatch(
    myTempFunction100(),
    error = function(err){"000111222"}
  )
  dev.off()
  
  
  
  sink( file=paste(path_temp1_sub1, "5A_Coverage.allSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  cat(   min(myCoverageMatrix_vector, na.rm = TRUE), "\t", mean(myCoverageMatrix_vector, na.rm = TRUE), "\t", 
         median(myCoverageMatrix_vector, na.rm = TRUE), "\t",  max(myCoverageMatrix_vector, na.rm = TRUE), "\t", "all", "\n"  )    
  sink()
  
  sink( file=paste(path_temp1_sub1, "5B_MeLevel.allSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  cat(   min(myMeLevelMatrix_vector, na.rm = TRUE), "\t", mean(myMeLevelMatrix_vector, na.rm = TRUE), "\t",  
         median(myMeLevelMatrix_vector, na.rm = TRUE), "\t",  max(myMeLevelMatrix_vector, na.rm = TRUE), "\t", "all", "\n"  )    
  sink()
  
  
  
  
  
  pdf( file=paste(path_temp1_sub1, "6A_Coverage.allSample.pdf", sep="/")  )
  print( hist(myCoverageMatrix_vector) )
  print( qplot(myCoverageMatrix_vector, binwidth=10)  )
  print( qplot(myCoverageMatrix_vector, binwidth=5) ) 
  print( qplot(myCoverageMatrix_vector, binwidth=1) )  
  dev.off() 
  
  pdf( file=paste(path_temp1_sub1, "6B_MeLevel.allSample.pdf", sep="/")  )
  print( hist(myMeLevelMatrix_vector) )
  print( qplot(myMeLevelMatrix_vector, binwidth=10)  )
  print( qplot(myMeLevelMatrix_vector, binwidth=5) ) 
  print( qplot(myMeLevelMatrix_vector, binwidth=1)  ) 
  dev.off() 
  
  pdf( file=paste(path_temp1_sub1, "7A_Coverage_eachSample.pdf", sep="/")  )
  ## par(mfrow=c(2,2))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myTempVec1 = myCoverageMatrix[,i]
    myTitle1 = as.vector( dataFrame_temp1$mysampleID )[i]
    print( hist(myTempVec1, main=myTitle1) )
    print( qplot(myTempVec1, binwidth=10, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=5, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=1, main=myTitle1)  )
  }
  dev.off()
  
  
  
  pdf( file=paste(path_temp1_sub1, "7B_MeLevel_eachSample.pdf", sep="/")  )
  ## par(mfrow=c(2,2))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myTempVec1 = myMeLevelMatrix[,i]
    myTitle1 = as.vector( dataFrame_temp1$mysampleID )[i]
    print( hist(myTempVec1, main=myTitle1) )
    print( qplot(myTempVec1, binwidth=10, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=5, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=1, main=myTitle1)  )
  }
  dev.off()
  
  
  
  
  
  
  

  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  tryCatch(
      MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=15,   height2=5 ),
      error = function(err){"000"}
  )
  
  ## MDS（multidimensional scaling）多维尺度分析
  ## Classical (Metric) Multidimensional Scaling
  path_temp3_sub6_g = paste(path_temp1, "3_MultidimensionalScaling", sep="/")
  if( ! file.exists(path_temp3_sub6_g) ) { dir.create(path_temp3_sub6_g, recursive = TRUE) }
  MyMultidimensionalScaling_1_g(  meLevelMatrix2=mat_2two,   path2=path_temp3_sub6_g,   dataFrame_temp2=dataFrame_temp1  )
  
  path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  tryCatch(
      MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=7,   height2=5,  dataFrame_temp2=dataFrame_temp1   ),
      error = function(err){"000"}
  )
                                                                 
  path_temp1_sub8 = paste(path_temp1, "5_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  tryCatch(
        myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  ) ,
        error = function(err){"000"}
  )
} 

myMainFunction1bp_1_g  <- function(  myobj_temp1,   path_temp1,   dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }
  
  meth_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  
                           treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  write.table(dataFrame_temp1 , 
              file = paste(path_temp1_sub1,   "0_dataFrame.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(meth_2two) )
  print( getSampleID(meth_2two) )
  print( getTreatment(meth_2two) )
  sink()
  
   
  mat_2two   = percMethylation( meth_2two )

  sink( file=paste(path_temp1_sub1 , "2_dimensions-tiles-merged.txt", sep="/")  )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  

  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(getData(meth_2two)[,1:4] , 
              file = paste(path_temp1_sub1,   "3C_regions-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  meth_sub4_matrix_g <- getData(meth_2two)
  mySamples_firstRow = colnames( mat_2two )
  myCoverageMatrix = mat_2two 
  myMeLevelMatrix  = mat_2two 
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myCoverageMatrix[,i] = meth_sub4_matrix_g[,3*i+2]
    myMeLevelMatrix[,i]  = 100*meth_sub4_matrix_g[,3*i+3]/meth_sub4_matrix_g[,3*i+2] 
  }
  
  
  pdf( file=paste(path_temp1_sub1, "4_Coverage_MeLevel_Cor.pdf", sep="/")  )
  par(mfrow=c(3,1))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    ## i =1 
    myTempVec1 = myCoverageMatrix[,i]
    myTempBool1 = (myTempVec1<=50) 
    myTempDataframe = data.frame(x1=myTempVec1[myTempBool1], y1=myMeLevelMatrix[myTempBool1,i] )
    boxplot(y1 ~ x1, data = myTempDataframe , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i] )  
    
    myTempBool2 = ( (myTempVec1 >= 50) & (myTempVec1 <= 90)  )
    myTempDataframe2 = data.frame(x2=myTempVec1[myTempBool2], y2=myMeLevelMatrix[myTempBool2,i] )
    boxplot(y2 ~ x2, data = myTempDataframe2 , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i]  ) 
    
    myTempBool2 =  (myTempVec1 >= 90)  
    myTempDataframe2 = data.frame(x2=myTempVec1[myTempBool2], y2=myMeLevelMatrix[myTempBool2,i] )
    boxplot(y2 ~ x2, data = myTempDataframe2 , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i]  ) 
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4A_Coverage.eachSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  for( i in c(1:length(as.character(dataFrame_temp1$myfiles))) ) {
    myTempCov1 = myCoverageMatrix[,i]
    file_name = as.character(dataFrame_temp1$myfiles)[i]
    cat(   min(myTempCov1, na.rm = TRUE), "\t", mean(myTempCov1, na.rm = TRUE), "\t", 
           median(myTempCov1, na.rm = TRUE), "\t",  max(myTempCov1, na.rm = TRUE), "\t", file_name, "\n"  )    
  }
  sink()
  
  sink( file=paste(path_temp1_sub1, "4B_MeLevel.eachSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  for( i in c(1:length(as.character(dataFrame_temp1$myfiles))) ) {
    myTempCov1 = myMeLevelMatrix[,i]
    file_name = as.character(dataFrame_temp1$myfiles)[i]
    cat(   min(myTempCov1, na.rm = TRUE), "\t", mean(myTempCov1, na.rm = TRUE), "\t",  
           median(myTempCov1, na.rm = TRUE), "\t",  max(myTempCov1, na.rm = TRUE), "\t", file_name, "\n"  )    
  }
  sink()
  
  
  myCoverageMatrix_vector = as.vector( myCoverageMatrix ) 
  myMeLevelMatrix_vector  = as.vector( myMeLevelMatrix ) 
  
  
  pdf( file=paste(path_temp1_sub1, "5_Coverage_MeLevel_Cor_poolAllSamples.pdf", sep="/")  )
  par(mfrow=c(3,1))
  myTempBool3 =  (myCoverageMatrix_vector<=40) 
  myTempDataframe3 = data.frame(x3=myCoverageMatrix_vector[myTempBool3], y3=myMeLevelMatrix_vector[myTempBool3])
  boxplot(y3 ~ x3, data = myTempDataframe3  , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples") 
  myTempBool4 =  ( (myCoverageMatrix_vector>40)  & (myCoverageMatrix_vector<=70) )
  myTempDataframe4 = data.frame(x4=myCoverageMatrix_vector[myTempBool4], y4=myMeLevelMatrix_vector[myTempBool4] )
  boxplot(y4 ~ x4, data = myTempDataframe4  , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples") 
  myTempBool5 =  ( (myCoverageMatrix_vector>70)  & (myCoverageMatrix_vector<=100) )
  myTempFunction100 <- function() {
  myTempDataframe5 = data.frame(x5=myCoverageMatrix_vector[myTempBool5], y5=myMeLevelMatrix_vector[myTempBool5] )
  boxplot(y5 ~ x5, data = myTempDataframe5 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
  myTempBool6 =  ( (myCoverageMatrix_vector>100)  & (myCoverageMatrix_vector<=130) )
  myTempDataframe6 = data.frame(x6=myCoverageMatrix_vector[myTempBool6], y6=myMeLevelMatrix_vector[myTempBool6] )
  boxplot(y6 ~ x6, data = myTempDataframe6 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
  myTempBool7 = (myCoverageMatrix_vector>130) 
  myTempDataframe7 = data.frame(x7=myCoverageMatrix_vector[myTempBool7], y7=myMeLevelMatrix_vector[myTempBool7] )
  boxplot(y7 ~ x7, data = myTempDataframe7 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
  }
  tryCatch(
    myTempFunction100(),
    error = function(err){"000111222"}
  )
  dev.off()

  
  
  sink( file=paste(path_temp1_sub1, "5A_Coverage.allSample.txt", sep="/")  )
    cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
    cat(   min(myCoverageMatrix_vector, na.rm = TRUE), "\t", mean(myCoverageMatrix_vector, na.rm = TRUE), "\t", 
           median(myCoverageMatrix_vector, na.rm = TRUE), "\t",  max(myCoverageMatrix_vector, na.rm = TRUE), "\t", "all", "\n"  )    
  sink()
  
  sink( file=paste(path_temp1_sub1, "5B_MeLevel.allSample.txt", sep="/")  )
    cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
    cat(   min(myMeLevelMatrix_vector, na.rm = TRUE), "\t", mean(myMeLevelMatrix_vector, na.rm = TRUE), "\t",  
           median(myMeLevelMatrix_vector, na.rm = TRUE), "\t",  max(myMeLevelMatrix_vector, na.rm = TRUE), "\t", "all", "\n"  )    
  sink()
  
  
  pdf( file=paste(path_temp1_sub1, "6A_Coverage.allSample.pdf", sep="/")  )
  print( hist(myCoverageMatrix_vector) )
  print( qplot(myCoverageMatrix_vector, binwidth=10)  )
  print( qplot(myCoverageMatrix_vector, binwidth=5) ) 
  print( qplot(myCoverageMatrix_vector, binwidth=1) )  
  dev.off() 
  
  pdf( file=paste(path_temp1_sub1, "6B_MeLevel.allSample.pdf", sep="/")  )
  print( hist(myMeLevelMatrix_vector) )
  print( qplot(myMeLevelMatrix_vector, binwidth=10)  )
  print( qplot(myMeLevelMatrix_vector, binwidth=5) ) 
  print( qplot(myMeLevelMatrix_vector, binwidth=1)  ) 
  dev.off() 
  
  pdf( file=paste(path_temp1_sub1, "7A_Coverage_eachSample.pdf", sep="/")  )
  ## par(mfrow=c(2,2))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myTempVec1 = myCoverageMatrix[,i]
    myTitle1 = as.vector( dataFrame_temp1$mysampleID )[i]
    print( hist(myTempVec1, main=myTitle1) )
    print( qplot(myTempVec1, binwidth=10, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=5, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=1, main=myTitle1)  )
  }
  dev.off()
  
  
  
  pdf( file=paste(path_temp1_sub1, "7B_MeLevel_eachSample.pdf", sep="/")  )
  ## par(mfrow=c(2,2))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myTempVec1 = myMeLevelMatrix[,i]
    myTitle1 = as.vector( dataFrame_temp1$mysampleID )[i]
    print( hist(myTempVec1, main=myTitle1) )
    print( qplot(myTempVec1, binwidth=10, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=5, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=1, main=myTitle1)  )
  }
  dev.off()
  
  
  
  
  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=15,   height2=5 )

  ## MDS（multidimensional scaling）多维尺度分析
  ## Classical (Metric) Multidimensional Scaling
  path_temp3_sub6_g = paste(path_temp1, "3_MultidimensionalScaling", sep="/")
  if( ! file.exists(path_temp3_sub6_g) ) { dir.create(path_temp3_sub6_g, recursive = TRUE) }
  MyMultidimensionalScaling_1_g(  meLevelMatrix2=mat_2two,   path2=path_temp3_sub6_g,   dataFrame_temp2=dataFrame_temp1  )
  
  path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=7,   height2=5,  dataFrame_temp2=dataFrame_temp1   )

  path_temp1_sub8 = paste(path_temp1, "5_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )

} 

myMainFunction_speRegions_g  <- function(  myobj_temp1,   path_temp1,  speRegions_temp1,   binBases_temp1,  dataFrame_temp1  ) {
  ## speRegions_temp1 = gene.obj_2000bp_g$promoters
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }
  
  meth_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  write.table(dataFrame_temp1 , 
              file = paste(path_temp1_sub1,   "0_dataFrame.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(meth_2two) )
  print( getSampleID(meth_2two) )
  print( getTreatment(meth_2two) )
  sink()  
  

  meth_2two = regionCounts(meth_2two,  speRegions_temp1,   cov.bases=binBases_temp1,  strand.aware=FALSE)  
  mat_2two  = percMethylation( meth_2two )
  

  sink( file=paste(path_temp1_sub1 , "2_dimensions-tiles-merged.txt", sep="/")  )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
 
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(getData(meth_2two)[,1:4] , 
              file = paste(path_temp1_sub1,   "3C_regions-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  

  
  
  meth_sub4_matrix_g <- getData(meth_2two)
  mySamples_firstRow = colnames( mat_2two )
  myCoverageMatrix = mat_2two 
  myMeLevelMatrix  = mat_2two 
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myCoverageMatrix[,i] = meth_sub4_matrix_g[,3*i+2]
    myMeLevelMatrix[,i]  = 100*meth_sub4_matrix_g[,3*i+3]/meth_sub4_matrix_g[,3*i+2] 
  }
  
  
  pdf( file=paste(path_temp1_sub1, "4_Coverage_MeLevel_Cor.pdf", sep="/")  )
  par(mfrow=c(3,1))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    ## i =1 
    myTempVec1 = myCoverageMatrix[,i]
    myTempBool1 = (myTempVec1<=70) 
    myTempDataframe = data.frame(x1=myTempVec1[myTempBool1], y1=myMeLevelMatrix[myTempBool1,i] )
    boxplot(y1 ~ x1, data = myTempDataframe , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i] )  
    
    myTempBool2 = ( (myTempVec1 >= 70) & (myTempVec1 <= 140)  )
    myTempDataframe2 = data.frame(x2=myTempVec1[myTempBool2], y2=myMeLevelMatrix[myTempBool2,i] )
    boxplot(y2 ~ x2, data = myTempDataframe2 , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i]  ) 
    
    myTempBool2 =  (myTempVec1 >= 140)  
    myTempDataframe2 = data.frame(x2=myTempVec1[myTempBool2], y2=myMeLevelMatrix[myTempBool2,i] )
    boxplot(y2 ~ x2, data = myTempDataframe2 , xlab="Coverage", ylab="Methylation Level (%)", main=mySamples_firstRow[i]  ) 
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4A_Coverage.eachSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  for( i in c(1:length(as.character(dataFrame_temp1$myfiles))) ) {
    myTempCov1 = myCoverageMatrix[,i]
    file_name = as.character(dataFrame_temp1$myfiles)[i]
    cat(   min(myTempCov1, na.rm = TRUE), "\t", mean(myTempCov1, na.rm = TRUE), "\t", 
           median(myTempCov1, na.rm = TRUE), "\t",  max(myTempCov1, na.rm = TRUE), "\t", file_name, "\n"  )    
  }
  sink()
  
  sink( file=paste(path_temp1_sub1, "4B_MeLevel.eachSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  for( i in c(1:length(as.character(dataFrame_temp1$myfiles))) ) {
    myTempCov1 = myMeLevelMatrix[,i]
    file_name = as.character(dataFrame_temp1$myfiles)[i]
    cat(   min(myTempCov1, na.rm = TRUE), "\t", mean(myTempCov1, na.rm = TRUE), "\t",  
           median(myTempCov1, na.rm = TRUE), "\t",  max(myTempCov1, na.rm = TRUE), "\t", file_name, "\n"  )    
  }
  sink()
  
  
  myCoverageMatrix_vector = as.vector( myCoverageMatrix ) 
  myMeLevelMatrix_vector  = as.vector( myMeLevelMatrix ) 
  
  
  pdf( file=paste(path_temp1_sub1, "5_Coverage_MeLevel_Cor_poolAllSamples.pdf", sep="/")  )
  par(mfrow=c(3,1))
  myTempBool3 =  (myCoverageMatrix_vector<=50) 
  myTempDataframe3 = data.frame(x3=myCoverageMatrix_vector[myTempBool3], y3=myMeLevelMatrix_vector[myTempBool3])
  boxplot(y3 ~ x3, data = myTempDataframe3  , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples") 
  myTempBool4 =  ( (myCoverageMatrix_vector>50)  & (myCoverageMatrix_vector<=100) )
  myTempDataframe4 = data.frame(x4=myCoverageMatrix_vector[myTempBool4], y4=myMeLevelMatrix_vector[myTempBool4] )
  boxplot(y4 ~ x4, data = myTempDataframe4  , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples") 
  myTempBool5 =  ( (myCoverageMatrix_vector>100)  & (myCoverageMatrix_vector<=150) )
  myTempFunction100 <- function() {
    myTempDataframe5 = data.frame(x5=myCoverageMatrix_vector[myTempBool5], y5=myMeLevelMatrix_vector[myTempBool5] )
    boxplot(y5 ~ x5, data = myTempDataframe5 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
    myTempBool6 =  ( (myCoverageMatrix_vector>150)  & (myCoverageMatrix_vector<=200) )
    myTempDataframe6 = data.frame(x6=myCoverageMatrix_vector[myTempBool6], y6=myMeLevelMatrix_vector[myTempBool6] )
    boxplot(y6 ~ x6, data = myTempDataframe6 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
    myTempBool7 = (myCoverageMatrix_vector>200) 
    myTempDataframe7 = data.frame(x7=myCoverageMatrix_vector[myTempBool7], y7=myMeLevelMatrix_vector[myTempBool7] )
    boxplot(y7 ~ x7, data = myTempDataframe7 , xlab="Coverage", ylab="Methylation Level (%)", main="All Samples" ) 
  }
  tryCatch(
    myTempFunction100(),
    error = function(err){"000111222"}
  )
  dev.off()
  
  
  
  sink( file=paste(path_temp1_sub1, "5A_Coverage.allSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  cat(   min(myCoverageMatrix_vector, na.rm = TRUE), "\t", mean(myCoverageMatrix_vector, na.rm = TRUE), "\t", 
         median(myCoverageMatrix_vector, na.rm = TRUE), "\t",  max(myCoverageMatrix_vector, na.rm = TRUE), "\t", "all", "\n"  )    
  sink()
  
  sink( file=paste(path_temp1_sub1, "5B_MeLevel.allSample.txt", sep="/")  )
  cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
  cat(   min(myMeLevelMatrix_vector, na.rm = TRUE), "\t", mean(myMeLevelMatrix_vector, na.rm = TRUE), "\t",  
         median(myMeLevelMatrix_vector, na.rm = TRUE), "\t",  max(myMeLevelMatrix_vector, na.rm = TRUE), "\t", "all", "\n"  )    
  sink()
  
  
  
  
  pdf( file=paste(path_temp1_sub1, "6A_Coverage.allSample.pdf", sep="/")  )
  print( hist(myCoverageMatrix_vector) )
  print( qplot(myCoverageMatrix_vector, binwidth=10)  )
  print( qplot(myCoverageMatrix_vector, binwidth=5) ) 
  print( qplot(myCoverageMatrix_vector, binwidth=1) )  
  dev.off() 
  
  pdf( file=paste(path_temp1_sub1, "6B_MeLevel.allSample.pdf", sep="/")  )
  print( hist(myMeLevelMatrix_vector) )
  print( qplot(myMeLevelMatrix_vector, binwidth=10)  )
  print( qplot(myMeLevelMatrix_vector, binwidth=5) ) 
  print( qplot(myMeLevelMatrix_vector, binwidth=1)  ) 
  dev.off() 
  
  pdf( file=paste(path_temp1_sub1, "7A_Coverage_eachSample.pdf", sep="/")  )
  ## par(mfrow=c(2,2))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myTempVec1 = myCoverageMatrix[,i]
    myTitle1 = as.vector( dataFrame_temp1$mysampleID )[i]
    print( hist(myTempVec1, main=myTitle1) )
    print( qplot(myTempVec1, binwidth=10, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=5, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=1, main=myTitle1)  )
  }
  dev.off()
  
  
  
  pdf( file=paste(path_temp1_sub1, "7B_MeLevel_eachSample.pdf", sep="/")  )
  ## par(mfrow=c(2,2))
  for(i  in  c(1: ncol(myCoverageMatrix)) )  {
    myTempVec1 = myMeLevelMatrix[,i]
    myTitle1 = as.vector( dataFrame_temp1$mysampleID )[i]
    print( hist(myTempVec1, main=myTitle1) )
    print( qplot(myTempVec1, binwidth=10, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=5, main=myTitle1)  )
    print( qplot(myTempVec1, binwidth=1, main=myTitle1)  )
  }
  dev.off()
  
  
  
  
  

  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=15,   height2=5 )
  
  ## MDS（multidimensional scaling）多维尺度分析
  ## Classical (Metric) Multidimensional Scaling
  path_temp3_sub6_g = paste(path_temp1, "3_MultidimensionalScaling", sep="/")
  if( ! file.exists(path_temp3_sub6_g) ) { dir.create(path_temp3_sub6_g, recursive = TRUE) }
  MyMultidimensionalScaling_1_g(  meLevelMatrix2=mat_2two,   path2=path_temp3_sub6_g,   dataFrame_temp2=dataFrame_temp1  )
  
  path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=7,   height2=5,  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub8 = paste(path_temp1, "5_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )
  
} 

myMainFunction_speRegions_All_g  <- function(  myobj_temp2,   path_temp2,  binBases_temp2,  dataFrame_temp2  )  {
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/1_promoters_2000bp",  sep="") , 
                                speRegions_temp1= gene.obj_2000bp_g$promoters,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/2_promoters_1500bp",  sep="") , 
                                speRegions_temp1= gene.obj_1500bp_g$promoters,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ) , 
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/3_promoters_1000bp",  sep="") , 
                                speRegions_temp1= gene.obj_1000bp_g$promoters,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ), 
  error = function(err){"000"}
  )
  
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/4_exons",  sep="") , 
  #                              speRegions_temp1= gene.obj_2000bp_g$exons ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  )
  
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/5_introns",  sep="") , 
  #                              speRegions_temp1= gene.obj_2000bp_g$introns ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-1A_imprintRegions",  sep="") , 
                                speRegions_temp1= imprint1.obj_g$ImprintedRegions ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-1B_imprintRegions-shores",  sep="") , 
                                speRegions_temp1= imprint1.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-2A_imprintRegions",  sep="") , 
                                speRegions_temp1= imprint2.obj_g$ImprintedRegions ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-2B_imprintRegions-shores",  sep="") , 
                                speRegions_temp1= imprint2.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-3A_imprintRegions",  sep="") , 
                                speRegions_temp1= imprint3.obj_g$ImprintedRegions ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-3B_imprintRegions-shores",  sep="") , 
                                speRegions_temp1= imprint3.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-4A_imprintRegions",  sep="") , 
                                speRegions_temp1= imprint4.obj_g$ImprintedRegions ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-4B_imprintRegions-shores",  sep="") , 
                                speRegions_temp1= imprint4.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-5A_imprintRegions",  sep="") , 
                                speRegions_temp1= imprint5.obj_g$ImprintedRegions ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/6-5B_imprintRegions-shores",  sep="") , 
                                speRegions_temp1= imprint5.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/7-A_H3K4me1Peaks",  sep="") , 
                                speRegions_temp1= region_1_H3K4me1.obj_g$H3K4me1 ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/7-B_H3K4me1-shores",  sep="") , 
                                speRegions_temp1= region_1_H3K4me1.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/8-A_H3K4me3Peaks",  sep="") , 
                                speRegions_temp1= region_2_H3K4me3.obj_g$H3K4me3 ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/8-B_H3K4me3-shores",  sep="") , 
                                speRegions_temp1= region_2_H3K4me3.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/9-A_H3K27acPeaks",  sep="") , 
                                speRegions_temp1= region_3_H3K27ac.obj_g$H3K27ac ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/9-B_H3K27ac-shores",  sep="") , 
                                speRegions_temp1= region_3_H3K27ac.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/10-A_H3K27me3Peaks",  sep="") , 
                                speRegions_temp1= region_4_H3K27me3.obj_g$H3K27me3 ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/10-B_H3K27me3-shores",  sep="") , 
                                speRegions_temp1= region_4_H3K27me3.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/11-A_centromeres",  sep="") , 
                                #speRegions_temp1= region_6_centromeres.obj_g$centromere ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  )
  
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/11-B_centromeres-shores",  sep="") , 
                                #speRegions_temp1= region_6_centromeres.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  )
  
  
  #tryCatch( 
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/12-A_HouseKeeping",  sep="") , 
  #                              speRegions_temp1= region_9_HouseKeeping.obj_g$HouseKeeping ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  #error = function(err){"000"}
  #)
  
  #tryCatch(
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/12-B_HouseKeeping-shores",  sep="") , 
  #                              speRegions_temp1= region_9_HouseKeeping.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  #error = function(err){"000"}
  #)
  
  
  #tryCatch(
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/13-A_rRNAGenes",  sep="") , 
  #                              speRegions_temp1= region_10_rRNA.obj_g$rRNA_Genes ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  #error = function(err){"000"}
  #)
  
  #tryCatch(
  #myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/13-B_rRNA-shores",  sep="") , 
  #                              speRegions_temp1= region_10_rRNA.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  #error = function(err){"000"}
  #)
  
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/14-A_activeEnhancers",  sep="") , 
                                speRegions_temp1= region_11A_activeEn.obj_g$activeEnhancers ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/14-B_activeEnhancers-shores",  sep="") , 
                                speRegions_temp1= region_11A_activeEn.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/15-A_otherEnhancers",  sep="") , 
                                speRegions_temp1= region_11B_otherEn.obj_g$otherEnhancers,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/15-B_otherEnhancers-shores",  sep="") , 
                                speRegions_temp1= region_11B_otherEn.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/16-A_posiedEnhancers",  sep="") , 
                                speRegions_temp1= region_11C_posiedEn.obj_g$posiedEnhancers ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/16-B_posiedEnhancers-shores",  sep="") , 
                                speRegions_temp1= region_11C_posiedEn.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/17-A_primedEnhancers",  sep="") , 
                                speRegions_temp1= region_11D_primedEn.obj_g$primedEnhancers ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  tryCatch(
  myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/17-B_primedEnhancers-shores",  sep="") , 
                                speRegions_temp1= region_11D_primedEn.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
  error = function(err){"000"}
  )
  
  
  
  
  
  tryCatch(
    myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/18-A_CpGislands",  sep="") , 
                                  speRegions_temp1= region_7_CpGislands.obj_g$CpGislands ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
    error = function(err){"000"}
  )
  
  tryCatch(
    myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/18-B_CpGislands-shores",  sep="") , 
                                  speRegions_temp1= region_7_CpGislands.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
    error = function(err){"000"}
  )
  
  
    
  
  
  
  tryCatch(
    myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/19-A_Repeats",  sep="") , 
                                  speRegions_temp1= myrepeat.obj_g$Repeats ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
    error = function(err){"000"}
  )
  
  tryCatch(
    myMainFunction_speRegions_g(  myobj_temp1 = myobj_temp2,   path_temp1 = paste(path_temp2, "/19-B_Repeats-shores",  sep="") , 
                                  speRegions_temp1= myrepeat.obj_g$shores ,   binBases_temp1=binBases_temp2,  dataFrame_temp1 = dataFrame_temp2  ),
    error = function(err){"000"}
  )
  
  

}
##################################################################################################################
cat("\n\n End part 2: \n\n\n\n\n")





cat("\n\n\n\n\n Start part 3: \n\n")
##################################################################################################################
## Annotating differentially methylated bases or regions
myRefSeqGenes_g       = "/home/yp/AnnotationBED/hg38/otherRegions/RefSeq_Genes.bed"
myRepeats_g           = "/home/yp/AnnotationBED/hg38/otherRegions/Repeats_rmsk.bed"
myImprintedRegions1_g = "/home/yp/AnnotationBED/hg38/imprintedRegions/67.Regions.PlosGenetics.ImprintedGenes.hg38.bed"
myImprintedRegions2_g = "/home/yp/AnnotationBED/hg38/imprintedRegions/75Regions.GR.hg38.bed"
myImprintedRegions3_g = "/home/yp/AnnotationBED/hg38/imprintedRegions/369Regions.GR.hg38.bed"
myImprintedRegions4_g = "/home/yp/AnnotationBED/hg38/imprintedRegions/merge1.imprintedRegions.hg38.bed"
myImprintedRegions5_g = "/home/yp/AnnotationBED/hg38/imprintedRegions/merge2.imprintedRegions.hg38.bed"
region_1_H3K4me1_g        = "/home/yp/AnnotationBED/hg38/1-H3K4me1.bed"
region_2_H3K4me3_g        = "/home/yp/AnnotationBED/hg38/2-H3K4me3.bed"
region_3_H3K27ac_g        = "/home/yp/AnnotationBED/hg38/3-H3K27ac.bed"
region_4_H3K27me3_g       = "/home/yp/AnnotationBED/hg38/4-H3K27me3.bed"
region_5_380imprint_g     = "/home/yp/AnnotationBED/hg38/5-380imprint.bed"
region_6_centromeres_g    = "/home/yp/AnnotationBED/hg38/6-centromeres.bed"
region_7_CpGislands_g     = "/home/yp/AnnotationBED/hg38/7-CpGislands.bed"
region_8A_activeEn_g      = "/home/yp/AnnotationBED/hg38/8A-enhancers-active-H3K4me1-H3K27ac.bed"
region_8B_otherEn_g       = "/home/yp/AnnotationBED/hg38/8B-enhancers-H3K4me1-H3K27ac-H3K27me3.bed"
region_8C_posiedEn_g      = "/home/yp/AnnotationBED/hg38/8C-enhancers-posied-H3K4me1-H3K27me3.bed"
region_8D_primedEn_g      = "/home/yp/AnnotationBED/hg38/8D-enhancers-primed-H3K4me1.bed"
region_9_HouseKeeping_g   = "/home/yp/AnnotationBED/hg38/9-HouseKeepingGenes.bed"
region_10_rRNA_g          = "/home/yp/AnnotationBED/hg38/10-rRNA-Genes.bed"
region_11A_activeEn_g     = "/home/yp/AnnotationBED/hg38/11A-active-H3K4me1-H3K27ac.bed"
region_11B_otherEn_g      = "/home/yp/AnnotationBED/hg38/11B-H3K4me1-H3K27ac-H3K27me3.bed"
region_11C_posiedEn_g     = "/home/yp/AnnotationBED/hg38/11C-posied-H3K4me1-H3K27me3.bed"
region_11D_primedEn_g     = "/home/yp/AnnotationBED/hg38/11D-primed-H3K4me1.bed"

gene.obj_2000bp_g   = readTranscriptFeatures(myRefSeqGenes_g, remove.unusual=TRUE, up.flank=2000,  down.flank=0,      unique.prom=TRUE)
gene.obj_1500bp_g   = readTranscriptFeatures(myRefSeqGenes_g, remove.unusual=TRUE, up.flank=1000,  down.flank=500,    unique.prom=TRUE)
gene.obj_1000bp_g   = readTranscriptFeatures(myRefSeqGenes_g, remove.unusual=TRUE, up.flank=1000,  down.flank=0,      unique.prom=TRUE)

myrepeat.obj_g = readFeatureFlank(myRepeats_g,           flank=2000, feature.flank.name=c("Repeats",          "shores"))
imprint1.obj_g = readFeatureFlank(myImprintedRegions1_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj_g = readFeatureFlank(myImprintedRegions2_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj_g = readFeatureFlank(myImprintedRegions3_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj_g = readFeatureFlank(myImprintedRegions4_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj_g = readFeatureFlank(myImprintedRegions5_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
region_1_H3K4me1.obj_g        = readFeatureFlank(region_1_H3K4me1_g,          flank=2000, feature.flank.name=c("H3K4me1",           "shores"))
region_2_H3K4me3.obj_g        = readFeatureFlank(region_2_H3K4me3_g,          flank=2000, feature.flank.name=c("H3K4me3",           "shores"))
region_3_H3K27ac.obj_g        = readFeatureFlank(region_3_H3K27ac_g,          flank=2000, feature.flank.name=c("H3K27ac",           "shores"))
region_4_H3K27me3.obj_g       = readFeatureFlank(region_4_H3K27me3_g,         flank=2000, feature.flank.name=c("H3K27me3",          "shores"))
region_5_380imprint.obj_g     = readFeatureFlank(region_5_380imprint_g,       flank=2000, feature.flank.name=c("imprint",           "shores"))
region_6_centromeres.obj_g    = readFeatureFlank(region_6_centromeres_g,      flank=2000, feature.flank.name=c("centromere",        "shores"))
region_7_CpGislands.obj_g     = readFeatureFlank(region_7_CpGislands_g,       flank=2000, feature.flank.name=c("CpGislands",        "shores"))
region_8A_activeEn.obj_g      = readFeatureFlank(region_8A_activeEn_g,        flank=2000, feature.flank.name=c("activeEnhancers",   "shores"))
region_8B_otherEn.obj_g       = readFeatureFlank(region_8B_otherEn_g,         flank=2000, feature.flank.name=c("otherEnhancers",    "shores"))
region_8C_posiedEn.obj_g      = readFeatureFlank(region_8C_posiedEn_g,        flank=2000, feature.flank.name=c("posiedEnhancers",   "shores"))
region_8D_primedEn.obj_g      = readFeatureFlank(region_8D_primedEn_g,        flank=2000, feature.flank.name=c("primedEnhancers",   "shores"))
region_9_HouseKeeping.obj_g   = readFeatureFlank(region_9_HouseKeeping_g,     flank=2000, feature.flank.name=c("HouseKeeping",      "shores"))
region_10_rRNA.obj_g          = readFeatureFlank(region_10_rRNA_g,            flank=2000, feature.flank.name=c("rRNA_Genes",        "shores"))
region_11A_activeEn.obj_g     = readFeatureFlank(region_11A_activeEn_g,        flank=2000, feature.flank.name=c("activeEnhancers",   "shores"))
region_11B_otherEn.obj_g      = readFeatureFlank(region_11B_otherEn_g,         flank=2000, feature.flank.name=c("otherEnhancers",    "shores"))
region_11C_posiedEn.obj_g     = readFeatureFlank(region_11C_posiedEn_g,        flank=2000, feature.flank.name=c("posiedEnhancers",   "shores"))
region_11D_primedEn.obj_g     = readFeatureFlank(region_11D_primedEn_g,        flank=2000, feature.flank.name=c("primedEnhancers",   "shores"))
##################################################################################################################
cat("\n\n End part 3: \n\n\n\n\n")





cat("\n\n\n\n\n Start part 4: \n\n")
##################################################################################################################
myOutDir_sub1_g = paste(outDir_g, "/2_ReadRawFiles",  sep="") 
if( ! file.exists(myOutDir_sub1_g) ) { dir.create(myOutDir_sub1_g, recursive = TRUE) }

sink( file=paste(myOutDir_sub1_g, "1_length-variables.txt", sep="/") )
length( Files_All_vector_g )
length( Files_All_list_g )
length( SampleID_All_vector_g )
length( SampleID_All_list_g )
length( Treatment_All_vector_g )
length( Treatment_All_list_g )
length( Sex_All_vector_g )
length( Sex_All_list_g )
length( Type_All_vector_g )
length( Type_All_list_g )
print( "#################### Files_All_vector_g ####################" )
print( Files_All_vector_g )
print( "#################### Files_All_list_g ####################" )
print( Files_All_list_g )
print( "#################### SampleID_All_vector_g ####################" )
print( SampleID_All_vector_g )
print( "#################### SampleID_All_list_g ####################" )
print( SampleID_All_list_g )
print( "#################### Treatment_All_vector_g ####################" )
print( Treatment_All_vector_g )
print( "#################### Treatment_All_list_g ####################" )
print( Treatment_All_list_g )
print( "#################### Sex_All_vector_g ####################" )
print( Sex_All_vector_g )
print( "#################### Sex_All_list_g ####################" )
print( Sex_All_list_g )
print( "#################### Type_All_vector_g ####################" )
print( Type_All_vector_g )
print( "#################### Type_All_list_g ####################" )
print( Type_All_list_g )
sink()



sink( file=paste(myOutDir_sub1_g, "2_theLog-of-read-inputFiles.txt", sep="/") )
myobj_g = methRead(location = Files_All_list_g, 
                   sample.id  = SampleID_All_list_g, 
                   assembly   = "hg38", 
                   pipeline   = "bismarkCoverage",
                   header     = FALSE, 
                   sep        = "\t", 
                   context    = "CpG",
                   resolution = "base",  ## allowed values 'base' or 'region'.
                   treatment  = Treatment_All_vector_g, 
                   mincov     = lowestCoverage_g
)
sink()


sink( file=paste(myOutDir_sub1_g, "3_all-rawFiles.txt", sep="/") )
    print(myobj_g)
sink()

 

sink( file=paste(myOutDir_sub1_g, "4_dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  print( "######################" )
  print(   Files_All_vector_g[i]  )
  print(   dim(myobj_g[[i]])  )
}
sink()


sink( file=paste(myOutDir_sub1_g, "5_dimensions-of-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(myobj_g[[i]])[1], "\t", dim(myobj_g[[i]])[2], "\t", Files_All_vector_g[i], "\n"  )    
}
sink()




sink( file=paste(myOutDir_sub1_g, "6_Check-each-file.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  ## i = 9
  temp_matrix1 <- read.table(file = Files_All_vector_g[i], header = FALSE, sep = "\t" )
  temp_matrix1_cov = rowSums(temp_matrix1[,c(5,6)])
  temp_matrix1 = temp_matrix1[temp_matrix1_cov>=lowestCoverage_g, ]
  temp_matrix2 <- getData(myobj_g[[i]])
  # length(temp_matrix1_cov)
  # length(temp_matrix1_cov[temp_matrix1_cov>=5])
  # dim(temp_matrix1)
  # dim(temp_matrix2)
  myNumCs_1 = temp_matrix1[,5]   
  myNumTs_1 = temp_matrix1[,6]  
  myNumCs_2 = temp_matrix2[,6]  
  myNumTs_2 = temp_matrix2[,7]  
  if( identical(myNumCs_1,myNumCs_2) & identical(myNumTs_1,myNumTs_2) ) {
    cat("Yes\n")
  }else{
    cat(" ##No## ",  ShortName_All_vector_g[i], "\n" )
  }
}
sink()


continue_on_error_g <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error_g) 



##########################1
filtered.myobj_g1 = filterByCoverage( myobj_g,  lo.count = lowestCoverage_g,  lo.perc = NULL,  hi.count = NULL,  
                                     hi.perc = 99,  chunk.size = 1e+06,  save.db = FALSE )        
filtered.myobj_g = filterByCoverage( filtered.myobj_g1,  lo.count = lowestCoverage_g,  lo.perc = NULL,  hi.count = 151,  
                                     hi.perc = NULL,  chunk.size = 1e+06,  save.db = FALSE )        

sink( file=paste(myOutDir_sub1_g, "7A_dimensions-of-eachCov.filtered.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(filtered.myobj_g[[i]])[1], "\t", dim(filtered.myobj_g[[i]])[2], "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()

sink( file=paste(myOutDir_sub1_g, "7B_coverageStat-of-eachCov.filtered.txt", sep="/")  )
cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
for( i in c(1:length(Files_All_vector_g)) ) {
  myTempCov1 = getData( filtered.myobj_g[[i]] )[,5]
  cat(   min(myTempCov1), "\t", mean(myTempCov1), "\t",  median(myTempCov1), "\t",  max(myTempCov1), "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()
##########################1


 
##########################1
removed.myobj_g1 = filterByCoverage( myobj_g,  lo.count = NULL,  lo.perc = 99,  hi.count = NULL,  
                                    hi.perc = NULL,  chunk.size = 1e+06,  save.db = FALSE )

myOutDir_sub2_g2 = paste(myOutDir_sub1_g, "/removed_1",  sep="") 
if( ! file.exists(myOutDir_sub2_g2) ) { dir.create(myOutDir_sub2_g2, recursive = TRUE) }

for( i in c(1:length(Files_All_vector_g)) ) {
  file_name = Files_All_vector_g[i]
  file_name = gsub(inputDir_g, i, file_name)
  file_name = gsub("/", "__", file_name)
  write.table(getData(removed.myobj_g1[[i]])  , 
              file = paste(myOutDir_sub2_g2,   file_name,  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
}

sink( file=paste(myOutDir_sub1_g, "8A_dimensions-removed1.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(removed.myobj_g1[[i]])[1], "\t", dim(removed.myobj_g1[[i]])[2], "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()

sink( file=paste(myOutDir_sub1_g, "8B_coverageStat-of-eachCov.removed1.txt", sep="/")  )
cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
for( i in c(1:length(Files_All_vector_g)) ) {
  myTempCov1 = getData( removed.myobj_g1[[i]] )[,5]
  cat(   min(myTempCov1), "\t", mean(myTempCov1), "\t",  median(myTempCov1), "\t",  max(myTempCov1), "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()
##########################1




##########################1
removed.myobj_g = filterByCoverage( filtered.myobj_g1,  lo.count = 151,  lo.perc = NULL,  hi.count = NULL,  
                                    hi.perc = NULL,  chunk.size = 1e+06,  save.db = FALSE )

myOutDir_sub2_g2 = paste(myOutDir_sub1_g, "/removed_2",  sep="") 
if( ! file.exists(myOutDir_sub2_g2) ) { dir.create(myOutDir_sub2_g2, recursive = TRUE) }

for( i in c(1:length(Files_All_vector_g)) ) {
  file_name = Files_All_vector_g[i]
  file_name = gsub(inputDir_g, i, file_name)
  file_name = gsub("/", "__", file_name)
  write.table(getData(removed.myobj_g[[i]])  , 
              file = paste(myOutDir_sub2_g2,   file_name,  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
}

sink( file=paste(myOutDir_sub1_g, "9A_dimensions-removed.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(removed.myobj_g[[i]])[1], "\t", dim(removed.myobj_g[[i]])[2], "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()

sink( file=paste(myOutDir_sub1_g, "9B_coverageStat-of-eachCov.removed.txt", sep="/")  )
cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
for( i in c(1:length(Files_All_vector_g)) ) {
  myTempCov1 = getData( removed.myobj_g[[i]] )[,5]
  cat(   min(myTempCov1), "\t", mean(myTempCov1), "\t",  median(myTempCov1), "\t",  max(myTempCov1), "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()
##########################1





##########################1
lowestCov.myobj_g = filterByCoverage( myobj_g,  lo.count = lowestCoverage_g,  lo.perc = NULL,  hi.count = NULL,  
                                      hi.perc = 25,  chunk.size = 1e+06,  save.db = FALSE )        

sink( file=paste(myOutDir_sub1_g, "10_dimensions-of-eachCov.lowestCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(lowestCov.myobj_g[[i]])[1], "\t", dim(lowestCov.myobj_g[[i]])[2], "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()
##########################1




sink( file=paste(myOutDir_sub1_g, "11_dimensions-all.txt", sep="/")  )
cat("Raw",  "Kept1", "Kept_final", "Removed1", "Removed2","Files", "\n", sep="\t")
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(    dim(myobj_g[[i]])[1], dim(filtered.myobj_g1[[i]])[1],  dim(filtered.myobj_g[[i]])[1],    
          dim(removed.myobj_g1[[i]])[1] , dim(removed.myobj_g[[i]])[1] , ShortName_All_vector_g[i], "\n" , sep="\t" )    
}
sink()



pdf( file=paste(myOutDir_sub1_g, "12_Methylation_Level_distribution.pdf", sep="/")  )
par(mfrow=c(2,3))
for( i in c(1:length(myobj_g)) ) {
  getMethylationStats(myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  getMethylationStats(filtered.myobj_g1[[i]], plot=TRUE, both.strands=FALSE )
  getMethylationStats(filtered.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  getMethylationStats(lowestCov.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  getMethylationStats(removed.myobj_g1[[i]], plot=TRUE, both.strands=FALSE )
  mytemp1 <- function() {
  getMethylationStats(removed.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  }
  tryCatch(
    mytemp1(), 
    error = function(err){ getMethylationStats(removed.myobj_g1[[i]], plot=TRUE, both.strands=FALSE ) }
  )
}
dev.off()




pdf( file=paste(myOutDir_sub1_g, "13_Coverage_distribution.pdf", sep="/")  )
par(mfrow=c(2,3))
for( i in c(1:length(myobj_g)) ) {
  getCoverageStats(myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  getCoverageStats(filtered.myobj_g1[[i]], plot=TRUE, both.strands=FALSE )
  getCoverageStats(filtered.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  getCoverageStats(lowestCov.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  getCoverageStats(removed.myobj_g1[[i]], plot=TRUE, both.strands=FALSE )
  mytemp1 <- function() {
  getCoverageStats(removed.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
  }
  tryCatch(
    mytemp1(), 
    error = function(err){ getCoverageStats(removed.myobj_g1[[i]], plot=TRUE, both.strands=FALSE ) }
  )
}
dev.off()





# myobj_nor_g <- filtered.myobj_g
myobj_nor_g <- normalizeCoverage(filtered.myobj_g,  method = "median")


sink( file=paste(myOutDir_sub1_g, "14A_dimensions-of-eachCov.normalized.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(myobj_nor_g[[i]])[1], "\t", dim(myobj_nor_g[[i]])[2], "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()

sink( file=paste(myOutDir_sub1_g, "14B_coverageStat-of-eachCov.normalized.txt", sep="/")  )
cat(   "min", "\t", "mean", "\t",  "median", "\t",  "max", "\t", "sample", "\n"  )    
for( i in c(1:length(Files_All_vector_g)) ) {
  myTempCov1 = getData( myobj_nor_g[[i]] )[,5]
  cat(   min(myTempCov1), "\t", mean(myTempCov1), "\t",  median(myTempCov1), "\t",  max(myTempCov1), "\t", ShortName_All_vector_g[i], "\n"  )    
}
sink()



myOutDir_sub2_g3 = paste(myOutDir_sub1_g, "/Segmentation",  sep="") 
if( ! file.exists(myOutDir_sub2_g3) ) { dir.create(myOutDir_sub2_g3, recursive = TRUE) }

pdf( file=paste(myOutDir_sub1_g, "15_segmentation.normalized.pdf", sep="/")  )
for( i in c(1:length(myobj_g)) ) {
   res_temp_1 = methSeg( myobj_nor_g[[i]] , diagnostic.plot=TRUE)
   file_name = Files_All_vector_g[i]
   file_name = gsub(inputDir_g, "", file_name)
   path_temp = paste(myOutDir_sub2_g3,   file_name,  sep="/")
   path_temp = gsub(".CpG.bismark.cov", "", path_temp)
   if( ! file.exists(path_temp) ) { dir.create(path_temp, recursive = TRUE) }
   methSeg2bed(res_temp_1,   filename= paste(myOutDir_sub2_g3, file_name, sep="/") )
}
dev.off()





my_methylBase_object  = unite( myobj_nor_g,  destrand=FALSE,   mc.cores=16 ,  min.per.group = CovSamples_g   )   

sink( file=paste(myOutDir_sub1_g, "16A_my_methylBase_object.merged.dimension.txt", sep="/")  )
print( dim(my_methylBase_object) )
cat("\n\n\n\n\n")
my_methylBase_object 
sink()

write.table(x=my_methylBase_object, 
            file=paste(myOutDir_sub1_g, "16B_my_methylBase_object.merged.txt", sep="/") , 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

sink( file=paste(myOutDir_sub1_g, "16C_my_methylBase_object.merged.info.txt", sep="/")  )
getSampleID(my_methylBase_object)
cat("\n\n")
getTreatment(my_methylBase_object)
sink()


rm(myobj_g)
rm(filtered.myobj_g1)
rm(filtered.myobj_g)
rm(removed.myobj_g1)
rm(removed.myobj_g)
rm(lowestCov.myobj_g)
rm(myobj_nor_g)
##################################################################################################################
cat("\n\n End part 4: \n\n\n\n\n")
 





cat("\n\n\n\n\n Start part 5: \n\n")
##################################################################################################################
dataFrame_temp1_A <- data.frame(
  mysampleID  = c(SampleID_one_1_g,   SampleID_two_2_g),
  mytreatment = c(Treatment_one_1_g,  Treatment_two_2_g),
  mysex       = c(Sex_one_1_g,        Sex_two_2_g),
  mytech      = c(Type_one_1_g,       Type_two_2_g), 
  myfiles     = c(Files_one_1_g,      Files_two_2_g)
)

myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
                        path_temp1      = paste(outDir_g, "/3A_",  Label_one_1_g, "_vs_", Label_two_2_g,  "_1bp",  sep="") ,      
                        dataFrame_temp1 = dataFrame_temp1_A  )

##################################################################################################################
cat("\n\n End part 5: \n\n\n\n\n")




 


cat("\n\n\n\n\n Start part 6: \n\n")
##################################################################################################################
dataFrame_temp2_A <- data.frame(
  mysampleID  = c(SampleID_one_1_g,   SampleID_three_3_g),
  mytreatment = c(Treatment_one_1_g,  Treatment_three_3_g),
  mysex       = c(Sex_one_1_g,        Sex_three_3_g),
  mytech      = c(Type_one_1_g,       Type_three_3_g), 
  myfiles     = c(Files_one_1_g,      Files_three_3_g)
)

myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
                        path_temp1      = paste(outDir_g, "/4A_",  Label_one_1_g, "_vs_", Label_three_3_g,  "_1bp",  sep="") ,      
                        dataFrame_temp1 = dataFrame_temp2_A  )
##################################################################################################################
cat("\n\n End part 6: \n\n\n\n\n")





 


cat("\n\n\n\n\n Start part 7: \n\n")
##################################################################################################################
dataFrame_temp3_A <- data.frame(
  mysampleID  = c(SampleID_one_1_g,   SampleID_four_4_g),
  mytreatment = c(Treatment_one_1_g,  Treatment_four_4_g),
  mysex       = c(Sex_one_1_g,        Sex_four_4_g),
  mytech      = c(Type_one_1_g,       Type_four_4_g), 
  myfiles     = c(Files_one_1_g,      Files_four_4_g)
)

myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
                        path_temp1      = paste(outDir_g, "/5A_",  Label_one_1_g, "_vs_", Label_four_4_g,  "_1bp",  sep="") ,      
                        dataFrame_temp1 = dataFrame_temp3_A  )
##################################################################################################################
cat("\n\n End part 7: \n\n\n\n\n")





cat("\n\n\n\n\n Start part 8: \n\n")
##################################################################################################################
dataFrame_temp4_A <- data.frame(
  mysampleID  = c(SampleID_two_2_g,   SampleID_three_3_g),
  mytreatment = c(Treatment_two_2_g,  Treatment_three_3_g),
  mysex       = c(Sex_two_2_g,        Sex_three_3_g),
  mytech      = c(Type_two_2_g,       Type_three_3_g), 
  myfiles     = c(Files_two_2_g,      Files_three_3_g)
)

myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
                        path_temp1      = paste(outDir_g, "/6A_",  Label_two_2_g, "_vs_", Label_three_3_g,  "_1bp",  sep="") ,      
                        dataFrame_temp1 = dataFrame_temp4_A  ) 
##################################################################################################################
cat("\n\n End part 8: \n\n\n\n\n")





 


cat("\n\n\n\n\n Start part 9: \n\n")
##################################################################################################################
dataFrame_temp5_A <- data.frame(
  mysampleID  = c(SampleID_two_2_g,   SampleID_four_4_g),
  mytreatment = c(Treatment_two_2_g,  Treatment_four_4_g),
  mysex       = c(Sex_two_2_g,        Sex_four_4_g),
  mytech      = c(Type_two_2_g,       Type_four_4_g), 
  myfiles     = c(Files_two_2_g,      Files_four_4_g)
)

myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
                        path_temp1      = paste(outDir_g, "/7A_",  Label_two_2_g, "_vs_", Label_four_4_g,  "_1bp",  sep="") ,      
                        dataFrame_temp1 = dataFrame_temp5_A  )
##################################################################################################################
cat("\n\n End part 9: \n\n\n\n\n")





 


cat("\n\n\n\n\n Start part 10: \n\n")
##################################################################################################################
dataFrame_temp6_A <- data.frame(
  mysampleID  = c(SampleID_three_3_g,   SampleID_four_4_g),
  mytreatment = c(Treatment_three_3_g,  Treatment_four_4_g),
  mysex       = c(Sex_three_3_g,        Sex_four_4_g),
  mytech      = c(Type_three_3_g,       Type_four_4_g), 
  myfiles     = c(Files_three_3_g,      Files_four_4_g)
)

myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
                        path_temp1      = paste(outDir_g, "/8A_",  Label_three_3_g, "_vs_", Label_four_4_g,  "_1bp",  sep="") ,      
                        dataFrame_temp1 = dataFrame_temp6_A  )
##################################################################################################################
cat("\n\n End part 10: \n\n\n\n\n")































##################################################################################################################
dataFrame_temp1_A <- data.frame(
  mysampleID  = c(SampleID_one_1_g,   SampleID_two_2_g),
  mytreatment = c(Treatment_one_1_g,  Treatment_two_2_g),
  mysex       = c(Sex_one_1_g,        Sex_two_2_g),
  mytech      = c(Type_one_1_g,       Type_two_2_g), 
  myfiles     = c(Files_one_1_g,      Files_two_2_g)
)

#myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
#                        path_temp1      = paste(outDir_g, "/3A_",  Label_one_1_g, "_vs_", Label_two_2_g,  "_1bp",  sep="") ,      
#                        dataFrame_temp1 = dataFrame_temp1_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/3B_",  Label_one_1_g, "_vs_", Label_two_2_g,  "_50bp",  sep="")  ,   
                     binSize_temp1 = 50,   binBases_temp1 = 1,   dataFrame_temp1 = dataFrame_temp1_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/3C_",  Label_one_1_g, "_vs_", Label_two_2_g,  "_100bp",  sep="")  ,   
                     binSize_temp1 = 100,   binBases_temp1 = 3,   dataFrame_temp1 = dataFrame_temp1_A  )

myMainFunction_speRegions_All_g(myobj_temp2 = my_methylBase_object,   
                                path_temp2 = paste(outDir_g, "/3D_",  Label_one_1_g, "_vs_", Label_two_2_g,  "_speRegions",  sep="") ,  
                                binBases_temp2=2,  dataFrame_temp2 = dataFrame_temp1_A  )  
##################################################################################################################





 



##################################################################################################################
dataFrame_temp2_A <- data.frame(
  mysampleID  = c(SampleID_one_1_g,   SampleID_three_3_g),
  mytreatment = c(Treatment_one_1_g,  Treatment_three_3_g),
  mysex       = c(Sex_one_1_g,        Sex_three_3_g),
  mytech      = c(Type_one_1_g,       Type_three_3_g), 
  myfiles     = c(Files_one_1_g,      Files_three_3_g)
)

#myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
#                        path_temp1      = paste(outDir_g, "/4A_",  Label_one_1_g, "_vs_", Label_three_3_g,  "_1bp",  sep="") ,      
#                        dataFrame_temp1 = dataFrame_temp2_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/4B_",  Label_one_1_g, "_vs_", Label_three_3_g,  "_50bp",  sep="")  ,   
                     binSize_temp1 = 50,   binBases_temp1 = 1,   dataFrame_temp1 = dataFrame_temp2_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/4C_",  Label_one_1_g, "_vs_", Label_three_3_g,  "_100bp",  sep="")  ,   
                     binSize_temp1 = 100,   binBases_temp1 = 3,   dataFrame_temp1 = dataFrame_temp2_A  )

myMainFunction_speRegions_All_g(myobj_temp2 = my_methylBase_object,   
                                path_temp2 = paste(outDir_g, "/4D_",  Label_one_1_g, "_vs_", Label_three_3_g,  "_speRegions",  sep="") ,  
                                binBases_temp2=2,  dataFrame_temp2 = dataFrame_temp2_A  )  
##################################################################################################################






 



##################################################################################################################
dataFrame_temp3_A <- data.frame(
  mysampleID  = c(SampleID_one_1_g,   SampleID_four_4_g),
  mytreatment = c(Treatment_one_1_g,  Treatment_four_4_g),
  mysex       = c(Sex_one_1_g,        Sex_four_4_g),
  mytech      = c(Type_one_1_g,       Type_four_4_g), 
  myfiles     = c(Files_one_1_g,      Files_four_4_g)
)

#myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
#                        path_temp1      = paste(outDir_g, "/5A_",  Label_one_1_g, "_vs_", Label_four_4_g,  "_1bp",  sep="") ,      
#                        dataFrame_temp1 = dataFrame_temp3_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/5B_",  Label_one_1_g, "_vs_", Label_four_4_g,  "_50bp",  sep="")  ,   
                     binSize_temp1 = 50,   binBases_temp1 = 1,   dataFrame_temp1 = dataFrame_temp3_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/5C_",  Label_one_1_g, "_vs_", Label_four_4_g,  "_100bp",  sep="")  ,   
                     binSize_temp1 = 100,   binBases_temp1 = 3,   dataFrame_temp1 = dataFrame_temp3_A  )

myMainFunction_speRegions_All_g(myobj_temp2 = my_methylBase_object,   
                                path_temp2 = paste(outDir_g, "/5D_",  Label_one_1_g, "_vs_", Label_four_4_g,  "_speRegions",  sep="") ,  
                                binBases_temp2=2,  dataFrame_temp2 = dataFrame_temp3_A  )  
##################################################################################################################







##################################################################################################################
dataFrame_temp4_A <- data.frame(
  mysampleID  = c(SampleID_two_2_g,   SampleID_three_3_g),
  mytreatment = c(Treatment_two_2_g,  Treatment_three_3_g),
  mysex       = c(Sex_two_2_g,        Sex_three_3_g),
  mytech      = c(Type_two_2_g,       Type_three_3_g), 
  myfiles     = c(Files_two_2_g,      Files_three_3_g)
)

#myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
#                        path_temp1      = paste(outDir_g, "/6A_",  Label_two_2_g, "_vs_", Label_three_3_g,  "_1bp",  sep="") ,      
#                        dataFrame_temp1 = dataFrame_temp4_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/6B_",  Label_two_2_g, "_vs_", Label_three_3_g,  "_50bp",  sep="")  ,   
                     binSize_temp1 = 50,   binBases_temp1 = 1,   dataFrame_temp1 = dataFrame_temp4_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/6C_",  Label_two_2_g, "_vs_", Label_three_3_g,  "_100bp",  sep="")  ,   
                     binSize_temp1 = 100,   binBases_temp1 = 3,   dataFrame_temp1 = dataFrame_temp4_A  )

myMainFunction_speRegions_All_g(myobj_temp2 = my_methylBase_object,   
                                path_temp2 = paste(outDir_g, "/6D_",  Label_two_2_g, "_vs_", Label_three_3_g,  "_speRegions",  sep="") ,  
                                binBases_temp2=2,  dataFrame_temp2 = dataFrame_temp4_A  )  
##################################################################################################################






 



##################################################################################################################
dataFrame_temp5_A <- data.frame(
  mysampleID  = c(SampleID_two_2_g,   SampleID_four_4_g),
  mytreatment = c(Treatment_two_2_g,  Treatment_four_4_g),
  mysex       = c(Sex_two_2_g,        Sex_four_4_g),
  mytech      = c(Type_two_2_g,       Type_four_4_g), 
  myfiles     = c(Files_two_2_g,      Files_four_4_g)
)

#myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
#                        path_temp1      = paste(outDir_g, "/7A_",  Label_two_2_g, "_vs_", Label_four_4_g,  "_1bp",  sep="") ,      
#                        dataFrame_temp1 = dataFrame_temp5_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/7B_",  Label_two_2_g, "_vs_", Label_four_4_g,  "_50bp",  sep="")  ,   
                     binSize_temp1 = 50,   binBases_temp1 = 1,   dataFrame_temp1 = dataFrame_temp5_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/7C_",  Label_two_2_g, "_vs_", Label_four_4_g,  "_100bp",  sep="")  ,   
                     binSize_temp1 = 100,   binBases_temp1 = 3,   dataFrame_temp1 = dataFrame_temp5_A  )

myMainFunction_speRegions_All_g(myobj_temp2 = my_methylBase_object,   
                                path_temp2 = paste(outDir_g, "/7D_",  Label_two_2_g, "_vs_", Label_four_4_g,  "_speRegions",  sep="") ,  
                                binBases_temp2=2,  dataFrame_temp2 = dataFrame_temp5_A  )  
##################################################################################################################






 



##################################################################################################################
dataFrame_temp6_A <- data.frame(
  mysampleID  = c(SampleID_three_3_g,   SampleID_four_4_g),
  mytreatment = c(Treatment_three_3_g,  Treatment_four_4_g),
  mysex       = c(Sex_three_3_g,        Sex_four_4_g),
  mytech      = c(Type_three_3_g,       Type_four_4_g), 
  myfiles     = c(Files_three_3_g,      Files_four_4_g)
)

#myMainFunction1bp_1_g(  myobj_temp1     = my_methylBase_object,   
#                        path_temp1      = paste(outDir_g, "/8A_",  Label_three_3_g, "_vs_", Label_four_4_g,  "_1bp",  sep="") ,      
#                        dataFrame_temp1 = dataFrame_temp6_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/8B_",  Label_three_3_g, "_vs_", Label_four_4_g,  "_50bp",  sep="")  ,   
                     binSize_temp1 = 50,   binBases_temp1 = 1,   dataFrame_temp1 = dataFrame_temp6_A  )

myMainFunction_1_g(  myobj_temp1 = my_methylBase_object,     
                     path_temp1 = paste(outDir_g, "/8C_",  Label_three_3_g, "_vs_", Label_four_4_g,  "_100bp",  sep="")  ,   
                     binSize_temp1 = 100,   binBases_temp1 = 3,   dataFrame_temp1 = dataFrame_temp6_A  )

myMainFunction_speRegions_All_g(myobj_temp2 = my_methylBase_object,   
                                path_temp2 = paste(outDir_g, "/8D_",  Label_three_3_g, "_vs_", Label_four_4_g,  "_speRegions",  sep="") ,  
                                binBases_temp2=2,  dataFrame_temp2 = dataFrame_temp6_A  )  
##################################################################################################################





