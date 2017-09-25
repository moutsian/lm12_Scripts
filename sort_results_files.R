trait="IBD"
for(chr in 1:22){
data=read.table(paste("results.",toupper(trait),".",chr,".txt",sep=""),head=T)
datao=data[order(as.numeric(as.character(data[,3]))),]
write.table(datao,paste("results.",toupper(trait),".",chr,".ordered.txt",sep=""),quote=F,row.names=F,col.names=T)
}
#END