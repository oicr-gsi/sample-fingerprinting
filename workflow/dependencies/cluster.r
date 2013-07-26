# R script for clustering genotypes by their similarity cumulative coverage results

cmd_args = commandArgs(trailingOnly = TRUE);

DATA<-read.table(cmd_args[1],header=T)

lim<-dim(DATA)[1]
HEAT<-heatmap(as.matrix(DATA[,1:lim]),Rowv=TRUE, Colv=TRUE, symm=TRUE, distfun=function(c) as.dist(1-cor(t(c), method="pearson")))
indexes = HEAT$rowInd
cat(colnames(DATA)[indexes],"\n")
