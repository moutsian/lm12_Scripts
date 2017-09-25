dataset="newwave"
#CHR=20
for(CHR in 1:22){
full_file=paste(CHR,".",dataset,".tokeep.txt",sep="")
updated_id_file=paste(CHR,".",dataset,".tokeep.updatedIDs.txt",sep="")
outfile=paste(CHR,".",dataset,".tokeep.txt.forplink",sep="")
full=read.table(full_file, stringsAsFactors=FALSE)
updated=read.table(updated_id_file, stringsAsFactors=FALSE)
#for(i in 1:100){
for(i in 1:dim(updated)[1]){
idx=which(full[,2]==updated[i,2])
if(length(idx)>0){
full[idx,3]=updated[i,3]
}}
write.table(full[,3],outfile,quote=F,col.names=F,row.names=F)
}
#END