#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.16;

## Perl5 version >= 5.16
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################



open(INPUT_FH,   "<",   "venn.txt"  )      or   die "$!"; 
open(OUTPUT_FH,  ">",   "VennDiagram_3sets_byPerl.R" )      or   die "$!"; 


my @lines = <INPUT_FH>; 
print   OUTPUT_FH   "library(VennDiagram)\n";


say($#lines);
for ( my $i=1; $i<=$#lines; $i++ ) {
        my $temp = $lines[$i];  
        $temp =~ m/^[\sX]+\s(\d+)\s/  or  die;
        my $number = $1;
        if($i==1)  {print   OUTPUT_FH   "myC = $number \n";}
        if($i==2)  {print   OUTPUT_FH   "myB = $number \n";}
        if($i==3)  {print   OUTPUT_FH   "myBC = $number \n";}
        if($i==4)  {print   OUTPUT_FH   "myA = $number \n";}
        if($i==5)  {print   OUTPUT_FH   "myAB = $number \n";}
        if($i==6)  {print   OUTPUT_FH   "myABC = $number \n";}
}








my $string_1 = <<END;   



x_A <- paste('A',  1:myA, sep='')
x_B <- paste('B',  1:myB, sep='')
x_C <- paste('C',  1:myC, sep='')

myBC = paste('BC',  1:myBC, sep='')
myAB = paste('AB',  1:myAB, sep='')
myABC = paste('ABC',  1:myABC, sep='')



x1 <- c(x_A, myAB, myABC  )
x2 <- c(x_B, myBC, myAB, myABC  )
x3 <- c(x_C, myBC, myABC  )


myFile <- "A_vs_B_vs_C.svg"
     
venn.diagram(x=list( "Analysis_3_random"=x3,  "Analysis_1"=x1,    "Analysis_2"=x2 ), 
             filename=myFile, 
             height = 3,   width = 5,  
             resolution =500,   imagetype = "svg",   units = "px", 
             compression ="lzw",  na = "stop",  
             main = "A_vs_B_vs_C", sub = NULL, 
             main.pos= c(0.5, 1),  main.fontface = "plain",
             main.fontfamily = "serif",  main.col ="black",
             main.cex = 1,  main.just = c(0.5, 1), 
             sub.pos = c(0.5, 1),  sub.fontface = "plain", 
             sub.fontfamily ="serif", sub.col = "black", sub.cex = 1,
             sub.just =c(0.5, 1), #category.names = names(x), 
             force.unique =TRUE, print.mode = "raw", sigdigs = 3, 
             direct.area =FALSE, area.vector = 0, hyper.test = TRUE,
             fill=c( "red3" , "skyblue", "blue3"), 
             cat.col=c( "red3" , "skyblue", "blue3" ), lwd=0, margin=0.2 )

 



END


print  OUTPUT_FH   $string_1;


close(OUTPUT_FH);

system(" Rscript    VennDiagram_3sets_byPerl.R ");  





