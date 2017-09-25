

###fix this!
all=read.table("step2.ALLTRAITS.dataset_B.all.round6.with_info_on_previous_hits_from_all_PAPERS.txt",head=T,sep="\t")
toremove=which(as.character(all[,19])=="")
allknown=all[-toremove,]
allknowno=allknown[order(as.numeric(as.character(allknown[,2])),as.numeric(as.character(allknown[,4]))),]
unique_locus=matrix(ncol=1,nrow=dim(allknowno)[1],NA)
unique_locus[1]=1
for(i in 2:dim(allknowno)[1]){
if( as.numeric(as.character(allknowno[i,2]))==as.numeric(as.character(allknowno[i-1,2])) & as.numeric(as.character(allknowno[i,8]))<=as.numeric(as.character(allknowno[i-1,9]))){
#print(paste("loci: ",allknowno[i,c(1:2,4:9)]," and ",allknowno[i-1,c(1:2,4:9)]," overlap",sep=""))
#unique_locus[i-1]=0
unique_locus[i]=0
allknowno[i,8]=min(as.numeric(as.character(allknowno[i-1,8])),as.numeric(as.character(allknowno[i,8])))
allknowno[i,9]=max(as.numeric(as.character(allknowno[i-1,9])),as.numeric(as.character(allknowno[i,9])))
}else{
unique_locus[i]=1
}
}

#unique_locus[i]=1

allknowno[which(unique_locus==1),c(1:2,4:9)]

write.table(allknowno[which(unique_locus==1),],"/nfs/users/nfs_l/lm12/IBD_conditional/all_replicated_loci_split_by_trait.datB",quote=F,col.names=T,row.names=F,sep="\t")

#END