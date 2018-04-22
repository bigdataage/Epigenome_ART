library(VennDiagram)


myA = 3698   # Parents_hyper
myB = 2996   # Parents_hypo
myC = 7924   # Children_hyper
myD = 5913    # Children_hypo

myBD = 2110
myBC = 99
myAD = 84
myAC = 2091


x_A <- paste('A',  1:myA, sep='')
x_B <- paste('B',  1:myB, sep='')
x_C <- paste('C',  1:myC, sep='')
x_D <- paste('D',  1:myD, sep='')
x_BD  <- paste('BD',  1:myBD, sep='')
x_BC  <- paste('BC',  1:myBC, sep='')
x_AD  <- paste('AD',  1:myAD, sep='')
x_AC  <- paste('AC',  1:myAC, sep='')

x1 <- c(x_A,   x_AC, x_AD )
x2 <- c(x_B,   x_BD, x_BC )
x3 <- c(x_C,   x_AC, x_BC )
x4 <- c(x_D,   x_BD, x_AD  ) 

myFile <- "A_vs_B_vs_C_vs_D.svg"

venn.diagram(x=list("Parents_hyper"=x1,   "Parents_hypo"=x2,  "Children_hyper"=x3,  "Children_hypo"=x4), 
             filename=myFile, 
             height = 3,   width = 4.5,  
             resolution =500,   imagetype = "svg",   units = "px", 
             compression ="lzw",  na = "stop",  
             main = "A_vs_B_vs_C_vs_D", sub = NULL, 
             main.pos= c(0.5, 1.05),  main.fontface = "plain",
             main.fontfamily = "serif",  main.col ="black",
             main.cex = 1,  main.just = c(0.5, 1), 
             sub.pos = c(0.5, 1.05),  sub.fontface = "plain", 
             sub.fontfamily ="serif", sub.col = "black", sub.cex = 1,
             sub.just =c(0.5, 1), #category.names = names(x), 
             force.unique =TRUE, print.mode = "raw", sigdigs = 3, 
             direct.area =FALSE, area.vector = 0, hyper.test = TRUE,
             fill=c("skyblue", "blue3", "pink", "red"),   cat.col=c("skyblue", "blue3", "pink", "red"),   lwd=0,   margin=0.1 )


