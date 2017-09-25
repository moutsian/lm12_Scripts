
trait="IBD"
dirin=paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/",trait,"/",sep="")
#chr="1"
data_all=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/IBDseq.",trait,".5e-08.filtered_by_QCsummaries.final.may.info0.8.txt",sep=""),head=T)
for(chr in 22:1){
if(file.exists(paste(dirin,"results.",trait,".",chr,".txt",sep=""))){
datain=read.table(paste(dirin,"results.",trait,".",chr,".txt",sep=""),head=T) #work on the .filtered data -> in new version we filter beforehand (23/03)
#data=read.table("datasetA/IBDseq.CD.5e-08.filtered_by_QCsummaries.final.may.info0.8.txt",head=T)
data=data_all[as.numeric(as.character(data_all[,1]))==chr,]
to_include=matrix(ncol=1,nrow=dim(datain)[1],0)

	for(i in 1:dim(datain)[1]){
		if(as.character(datain[i,2])%in%as.character(data[,2])){
		to_include[i]=1
		}
	}
	
write.table(datain[which(to_include==1),],paste(dirin,"results.filtered.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)

}
}

#END
