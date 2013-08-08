# R script for making fingerprint glyph
# needs *sample_interval_summary files from GATK

cmd_args = commandArgs(trailingOnly = TRUE);

# Need to get list of files and cluster colors assigned to these files

dir     = cmd_args[1] # Directory with data files
files   = cmd_args[2] # Comma-separated files with letter codes at 'hotspot' positions (preperaed with GATK DeapthOfCoverage results)
ccols   = cmd_args[3] # Cluster colors (From R)
snps    = as.numeric(cmd_args[4]) # Number of control 'hotspots'
pngname = cmd_args[5] # File name to wite to

print(paste("Got arguments: files=",files,"\ncolors=",ccols,"\nsnps=",snps,"\nfile name=",pngname))

names<-strsplit(files,",",fixed=T)
cclrs<-strsplit(ccols,",",fixed=T)
refcols<-rbind(c("#00960A","#0000FF","#D17105","#FF0000","#7C7C7C"))

if (length(names[[1]]) != length(cclrs[[1]])) {
 print("Error: we don't have as many files as cluster color values");
 q(status = 1)
}

# Find what is the longest name, in pixels (6px per char at cex 1.0)
maxchar = 0
for(n in names[[1]]) {
 maxchar<-ifelse(nchar(n)>maxchar,nchar(n),maxchar)
}

# Set the width of image an right margin
#off<-ifelse(maxchar*6.5>snps+10,maxchar*6.5,0)
off = 0

# Function for drawing colorstrips
colorstrip<-function(cls,ccol,base,height,margin=0) {
  rect(0, base+margin, 5, base+height-margin, border=NA, col=ccol, lty="solid")
 for(c in 1:length(cls)){
  rect(c+10, base+margin, c+11, base+height-margin, border=NA, col=cls[c], lty="solid")
 }
}

# Set the dimensions of the image
image_height = 20*length(names[[1]])

# Assign colnames for color array

colnames(refcols)<-c("A","C","G","T","M")
BLANC="#D7D7D7"

print(paste("We will have ",length(names[[1]]),"rows in our diagram"))

png(filename=pngname,width=ifelse(off==0,snps+10,off),height=image_height,units="px",pointsize=15,bg="white")
par(mar=c(0,0,0,0)+0.1,family="sans",mfrow=c(length(names[[1]]),1))

for(i in 1:length(names[[1]])){
 
 FP<-read.table(paste(dir,names[[1]][i],sep=""),header=T,sep="\t")
 DATA<-data.frame(FP$FLAG,rep(BLANC,length(FP$FLAG)))
 colnames(DATA)<-c("Flag","Color")
 C<-as.vector(DATA$Color)

 for( d in 1:length(DATA$Flag)) {
  for(c in 1:length(colnames(refcols))) {
   if(DATA$Flag[d]==colnames(refcols)[c]) {
       C[d]<-refcols[c]
       break
   }
  }
 }

 name = substr(names[[1]][i],0,nchar(names[[1]][i])-4)
 print(paste("Name for ",i," is ",name))
 
 plot(1, type="n", axes=F, xlab="", ylab="",xlim=c(0,ifelse(off==0,snps+10,off)),ylim=c(0,10),yaxs = 'i',xaxs = 'i')
 colorstrip(C,cclrs[[1]][i],0,10,0)
 #text(5,1,name,pos=4,cex=1,font=2)
 
}
dev.off()
