# R script for clustering genotypes by their similarity cumulative coverage results

cmd_args = commandArgs(trailingOnly = TRUE);

DATA<-read.table(cmd_args[1],header=T,check.names=FALSE)

# Sort the matrix, it saves a lot of headaches
ordcol<-order(colnames(DATA[,1:dim(DATA)[1]]))
ordrow<-order(rownames(DATA))
SORTED<-cbind(DATA[,ordcol],SNPs=DATA[,dim(DATA)[1] + 1])
DATA<-SORTED[ordrow,]

lim<-dim(DATA)[1]
HEAT<-heatmap(as.matrix(DATA[,1:lim]),Rowv=TRUE, Colv=TRUE, symm=TRUE, scale='none', distfun=function(c) as.dist(1-cor(t(c), method="pearson")))
indexes = HEAT$rowInd
cat(colnames(DATA)[indexes],"\n")
