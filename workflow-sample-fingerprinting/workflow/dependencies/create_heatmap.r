# R script for clustering genotypes by their similarity cumulative coverage results

cmd_args = commandArgs(trailingOnly = TRUE);

f       = cmd_args[1]             # matrix file
title   = cmd_args[2]             # title of the heatmap (list of samples)
maxSNPs = as.numeric(cmd_args[3]) # how many control 'hotspots' (ref SNPs)
pngname = cmd_args[4]             # Name of the pngfile for writing into

image_width  = as.numeric(cmd_args[5])
image_height = image_width

flagged = cmd_args[6]
dendro  = cmd_args[7]


plotfile<-function(f,title,maxSNPs,flagged=FALSE,dendro=FALSE) {
 maincol<-ifelse(flagged,"red","black")
 par(cex.main=0.7,col.main=maincol)
 DATA<-read.table(f,header=T)
 lim<-dim(DATA)[1]

 # Sort the matrix, it saves a lot of headaches
 ordcol<-order(colnames(DATA[,1:dim(DATA)[1]]))
 ordrow<-order(rownames(DATA))
 SORTED<-cbind(DATA[,ordcol],SNPs=DATA[,lim + 1],Color=DATA[,lim + 2])
 DATA<-SORTED[ordrow,]

 COLFUN <- colorRamp(c("white","lightblue","blue"),space="rgb")
 NUMFUN<-colorRamp(c("white","black"),space="rgb")
 normSNPs<-with(DATA,(SNPs) / max(SNPs))


 HEAT<-heatmap(as.matrix(DATA[,1:lim]),na.rm= TRUE,Rowv=TRUE,Colv=TRUE,col=rgb(COLFUN(seq(from=0,to=1,by=0.01)),maxColorValue=255),symm=TRUE, cexRow=0.9, cexCol=0.9, labRow=FALSE, labCol=FALSE, main=paste("Genotype proximity for ",title),RowSideColors=as.character(DATA$Color),ColSideColors=rgb(NUMFUN(normSNPs),maxColorValue=255),distfun=function(c) as.dist(1-cor(t(c), method="pearson")),margins=c(7,7),cex.main=0.6,bg="salmon",keep.dendro=TRUE)
 
 if (dendro) {
   utils::str(HEAT$Rowv)
 } else {
   indexes = HEAT$rowInd
   cat(rev(rownames(DATA)[indexes]),"\n")
 }

}

# Make png only if we have valid name (ending with png) - Useful when requesting only dendrogram
if (grepl(".png$", pngname)) { 
 png(filename=pngname,width=image_width,height=image_height,units="px",pointsize=15,bg="white")
}
 
plotfile(f,title,maxSNPs,flagged,dendro)

if (grepl(".png$", pngname)) {
 blah<-dev.off()
}