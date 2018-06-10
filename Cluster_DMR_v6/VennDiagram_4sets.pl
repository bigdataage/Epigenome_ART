#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.16;

## Perl5 version >= 5.16
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################



open(INPUT_FH,   "<",   "venn.txt"  )      or   die "$!"; 
open(OUTPUT_FH,  ">",   "VennDiagram_4sets_byPerl.R" )      or   die "$!"; 


my @lines = <INPUT_FH>; 
print   OUTPUT_FH   "library(VennDiagram)\n";


say($#lines);
for ( my $i=1; $i<=$#lines; $i++ ) {
        my $temp = $lines[$i];  
        $temp =~ m/^[\sX]+\s(\d+)\s/  or  die;
        my $number = $1;
        if($i==1)  {print   OUTPUT_FH   "myD = $number \n";}
        if($i==2)  {print   OUTPUT_FH   "myC = $number \n";}
        if($i==3)  {print   OUTPUT_FH   "myCD = $number \n";}
        if($i==4)  {print   OUTPUT_FH   "myB = $number \n";}
        if($i==5)  {print   OUTPUT_FH   "myBD = $number \n";}
        if($i==6)  {print   OUTPUT_FH   "myBCD = $number \n";}
        if($i==7)  {print   OUTPUT_FH   "myA = $number \n";}
        if($i==8)  {print   OUTPUT_FH   "myAC = $number \n";}
        if($i==9)  {print   OUTPUT_FH   "myACD = $number \n";}
        if($i==10) {print   OUTPUT_FH   "myAB = $number \n";}
        if($i==11) {print   OUTPUT_FH   "myABD = $number \n";}
        if($i==12) {print   OUTPUT_FH   "myABC = $number \n";}
        if($i==13) {print   OUTPUT_FH   "myABCD = $number \n";}
}








my $string_1 = <<END;   



x_A <- paste('A',  1:myA, sep='')
x_B <- paste('B',  1:myB, sep='')
x_C <- paste('C',  1:myC, sep='')
x_D <- paste('D',  1:myD, sep='')
myCD = paste('CD',  1:myCD, sep='')
myBD = paste('BD',  1:myBD, sep='')
myBCD = paste('BCD',  1:myBCD, sep='')
myAC = paste('AC',  1:myAC, sep='')
myACD = paste('ACD',  1:myACD, sep='')
myAB = paste('AB',  1:myAB, sep='')
myABD = paste('ABD',  1:myABD, sep='')
myABC = paste('ABC',  1:myABC, sep='')
myABCD = paste('ABCD',  1:myABCD, sep='')

x1 <- c(x_A, myAC, myACD, myAB,  myABD, myABC, myABCD   )
x2 <- c(x_B, myBD, myBCD, myAB,  myABD, myABC, myABCD  )
x3 <- c(x_C, myCD, myBCD, myAC,  myACD, myABC, myABCD  )
x4 <- c(x_D, myCD, myBD,  myBCD, myACD, myABD, myABCD   ) 

myFile <- "A_vs_B_vs_C_vs_D.svg"

venn.diagram(x=list("hyper1"=x1,   "hyper2"=x2,  "hyper3"=x3,  "hyper4"=x4), 
             filename=myFile, 
             height = 3,   width = 7,  
             resolution =500,   imagetype = "svg",   units = "px", 
             compression ="lzw",  na = "stop",  
             main = "A_vs_B_vs_C", sub = NULL, 
             main.pos= c(0.5, 1.05),  main.fontface = "plain",
             main.fontfamily = "serif",  main.col ="black",
             main.cex = 1,  main.just = c(0.5, 1), 
             sub.pos = c(0.5, 1.05),  sub.fontface = "plain", 
             sub.fontfamily ="serif", sub.col = "black", sub.cex = 1,
             sub.just =c(0.5, 1), #category.names = names(x), 
             force.unique =TRUE, print.mode = "raw", sigdigs = 3, 
             direct.area =FALSE, area.vector = 0, hyper.test = TRUE,
fill=c("skyblue", "blue3", "pink", "red3"), cat.col=c("skyblue", "blue3", "pink", "red3"), lwd=0, margin=0.1 )

 



END


print  OUTPUT_FH   $string_1;


close(OUTPUT_FH);

system(" Rscript    VennDiagram_4sets_byPerl.R ");  





