## Rscript  diffMe_select.R  3A_ICSI-fresh1_vs_ICSI-fresh2_1bp

args_g <- commandArgs(TRUE)
print("##########################")
print("args: ")
print(args_g[1])   
print("##########################")

nameOfRegion = args_g[1];     ## name of genomic region

# nameOfRegion =  "3A_ICSI-fresh1_vs_ICSI-fresh2_1bp"



matrix1 = read.table(file="2E_AlldiffMesites_q0.05_diff0.txt", header = TRUE, sep = "\t", quote = "\"'", dec = ".")
dim(matrix1)
colnames(matrix1)

myQvalue = matrix1$qvalue
myDiff   = matrix1$meth.diff   




print("#####################################################")
print("## Number of DMCs for six conditions: ")
nrow( matrix1[( (myQvalue<0.001) & (abs(myDiff)>10) ),]  )
nrow( matrix1[( (myQvalue<0.01)  & (abs(myDiff)>10) ),]  )
nrow( matrix1[( (myQvalue<0.05)  & (abs(myDiff)>10) ),]  )
nrow( matrix1[( (myQvalue<0.001) & (abs(myDiff)>5) ),]  )
nrow( matrix1[( (myQvalue<0.01)  & (abs(myDiff)>5) ),]  )
nrow( matrix1[( (myQvalue<0.05)  & (abs(myDiff)>5) ),]  )
print("#####################################################")


## Select the DMCs.
diff1_all = matrix1[( (myQvalue<0.05)  & (abs(myDiff)>5) ),]
print("#####################################################")
print("## Dimension of diff1_all : ")
dim(diff1_all) 
print("#####################################################")


diff1_name = c()
for(i in c(1:nrow(diff1_all)) ) {
  diff1_name[i] = paste(nameOfRegion, i, sep="_"  )     
}




############### for ChIPseeker
diff1_ChIPseeker = cbind(diff1_all[,1:4],  diff1_name, rep(1, nrow(diff1_all)),   diff1_all[,5:7] )   
colnames(diff1_ChIPseeker) = c("chr", "start",   "end",  "strand", "name",  "strength",  "pvalue",    "qvalue",  "meth.diff")
diff1_ChIPseeker[,3] = diff1_ChIPseeker[,3] + 1

diff1_ChIPseeker_hyper = diff1_ChIPseeker[ diff1_ChIPseeker$meth.diff > 0, ]
diff1_ChIPseeker_hypo  = diff1_ChIPseeker[ diff1_ChIPseeker$meth.diff < 0, ]

print("#####################################################")
print("diff1_ChIPseeker:")
dim(diff1_ChIPseeker) 
dim(diff1_ChIPseeker_hyper) 
dim(diff1_ChIPseeker_hypo) 
print("#####################################################")


 




############### for BED
diff1_BED = cbind(diff1_all[,1:3],  diff1_name )   
colnames(diff1_BED) = c("chr", "start",   "end",  "name" )
diff1_BED[,3] = diff1_BED[,3] + 1

diff1_BED_hyper = diff1_BED[ diff1_ChIPseeker$meth.diff > 0, ]
diff1_BED_hypo  = diff1_BED[ diff1_ChIPseeker$meth.diff < 0, ]

print("#####################################################")
print("diff1_BED:")
dim(diff1_BED) 
dim(diff1_BED_hyper) 
dim(diff1_BED_hypo) 
print("#####################################################")







############### for BED (not plus 1 for end site)
diff1_BED_keptRaw = cbind(diff1_all[,1:3],  diff1_name )   
colnames(diff1_BED_keptRaw) = c("chr", "start",   "end",  "name" )

diff1_BED_keptRaw_hyper = diff1_BED_keptRaw[ diff1_ChIPseeker$meth.diff > 0, ]
diff1_BED_keptRaw_hypo  = diff1_BED_keptRaw[ diff1_ChIPseeker$meth.diff < 0, ]

print("#####################################################")
print("diff1_BED_keptRaw:")
dim(diff1_BED_keptRaw) 
dim(diff1_BED_keptRaw_hyper) 
dim(diff1_BED_keptRaw_hypo) 
print("#####################################################")







Part1_g = "ChIPseeker/1-rawFiles"
if( ! file.exists(Part1_g) )  { dir.create(Part1_g, recursive = TRUE) }
write.table(x=diff1_ChIPseeker,       file = paste(Part1_g, "diff1_all.txt", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff1_ChIPseeker_hyper, file = paste(Part1_g, "diff1_hyper.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff1_ChIPseeker_hypo,  file = paste(Part1_g, "diff1_hypo.txt", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )




Part1_g = "BED"
if( ! file.exists(Part1_g) )  { dir.create(Part1_g, recursive = TRUE) }
write.table(x=diff1_BED,       file = paste(Part1_g, "diff1_all.bed", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
write.table(x=diff1_BED_hyper, file = paste(Part1_g, "diff1_hyper.bed", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
write.table(x=diff1_BED_hypo,  file = paste(Part1_g, "diff1_hypo.bed", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )





Part1_g = "BED_keptRaw"
if( ! file.exists(Part1_g) )  { dir.create(Part1_g, recursive = TRUE) }
write.table(x=diff1_BED_keptRaw,       file = paste(Part1_g, "diff1_all.raw.bed", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
write.table(x=diff1_BED_keptRaw_hyper, file = paste(Part1_g, "diff1_hyper.raw.bed", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
write.table(x=diff1_BED_keptRaw_hypo,  file = paste(Part1_g, "diff1_hypo.raw.bed", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )









