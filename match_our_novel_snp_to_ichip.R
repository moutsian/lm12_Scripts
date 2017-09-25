novel=read.table("~lm12/IBD_conditional/step2.NOVEL.from_all_papers.txt",head=T,sep="\t")
ichip=read.table("~lm12/IBD_conditional/EUR.IBD.ichip_only.MMM.assoc.clean",head=T,sep="\t")
liftedover=read.table("~lm12/IBD_conditional/EUR.IBD.ichip_only.MMM.assoc.clean.liftedover",head=F)
ichip[,3]=liftedover[,2] #now ichip data is in the right build (same as our novel loci)
novel_plus=matrix(ncol=dim(novel)[2]+4,nrow=dim(novel)[1],NA)
novel_plus[,1:dim(novel)[2]]=as.matrix(novel)
for(i in 1:dim(novel)[1]){
	idx=which(as.numeric(as.character(ichip[,1]))==as.numeric(as.character(novel[i,2])) & as.numeric(as.character(ichip[,3]))==as.numeric(as.character(novel[i,4])) )
	if(length(idx)>0){
	novel_plus[i,20:23]=as.matrix(ichip[idx[1],6:9])
	}
	print(paste("i: ",i," idx: ",idx,sep=""))
}
colnames(novel_plus)[20:23]=c("FREQ1","OR","SE","P")
write.table(novel_plus,"~lm12/IBD_conditional/step2.NOVEL.from_all_papers.with_ichip_info.txt",quote=F,col.names=T,row.names=F,sep="\t")
#END