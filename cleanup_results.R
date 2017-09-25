
trait="UC"
dirin=paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",trait,"/",sep="")
#chr="1"
for(chr in 22:1){
if(file.exists(paste(dirin,"results.",trait,".",chr,".txt",sep=""))){
datain=read.table(paste(dirin,"results.",trait,".",chr,".txt",sep=""),head=T) #work on the .filtered data -> in new version we filter beforehand (23/03)
data_all=read.table(paste(dirin,chr,"-meta.filtered.LM.txt",sep=""),head=T)
to_include=matrix(ncol=1,nrow=dim(datain)[1],0)

	for(i in 1:dim(datain)[1]){
		if(as.character(datain[i,2])%in%as.character(data_all[,1])){
		to_include[i]=1
		}
	}
	
write.table(datain[which(to_include==1),],paste(dirin,"results.filtered.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)

}
}

#END
