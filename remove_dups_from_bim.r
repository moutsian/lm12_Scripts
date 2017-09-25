# this takes in a bim file and removes duplicate entries, either at random (if the frequencies are very similar), or both (if the alleles or frequencies
# disagree while the two entries have the same rsIDs).
# Output is the new and old bim files, as well as the list of variants to remove  
bimfilename="ichip_b37.final.xMHC.bim"
bimfile=read.table(bimfilename,stringsAsFactors=F)
n_occur <- data.frame(table(bimfile[,2]))
dups=n_occur[n_occur$Freq>1,]
frq=read.table("plink.frq",head=T,stringsAsFactors=F)
list_to_remove=NULL
newbim=bimfile
for(i in 1:dim(dups)[1]){
dupindex=which(frq[,2]%in%dups[i,1])
#the following assumes there are no triplicates etc
if(frq[dupindex[1],3]==frq[dupindex[2],3] & frq[dupindex[1],4]==frq[dupindex[2],4] )
{
	#check if freqs are similar
	freqdiff=abs(frq[dupindex[1],5]-frq[dupindex[2],5])/max(frq[dupindex[1],5],frq[dupindex[2],5])
	if(is.nan(freqdiff) | freqdiff<0.15){
	#keep one at random
	newbim[dupindex[1],2]=paste(newbim[dupindex[2],2],"toremove",sep="")
	list_to_remove=c(list_to_remove,newbim[dupindex[1],2])
	}else{
	#remove both
	list_to_remove=c(list_to_remove,as.character(dups[i,1]))}
}else{
#just remove both
list_to_remove=c(list_to_remove,as.character(dups[i,1]))}
}
write.table(newbim,bimfilename,quote=F,col.names=F,row.names=F)
write.table(bimfile,paste(bimfilename,".with_dups",sep=""),quote=F,col.names=F,row.names=F) #just to save the old file too
write.table(list_to_remove,paste(bimfilename,".dups_to_remove",sep=""),sep="",quote=F,col.names=F,row.names=F)

#END