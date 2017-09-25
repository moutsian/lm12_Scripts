#use this file to filter the results of datA: /lustre/scratch114/projects/crohns/RELEASE/QCsummaries
trait="IBD"
TRAIT=toupper(trait)
suggestive=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/IBDseq.",TRAIT,".5e-08.filtered.txt",sep=""),head=F)
suggestiveo=suggestive[order(as.numeric(as.character(suggestive[,1])),as.numeric(as.character(suggestive[,3]))),]

for(chr in 1:22){
print(paste("chr:",chr))
subdat=suggestiveo[suggestiveo[,1]==chr,]
qc_data=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/QC_svmfiltered.17M.",chr,".txt",sep=""),head=F)
tokeep=which(as.character(subdat[,2])%in%as.character(qc_data[,1]))
write.table(subdat[tokeep,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/IBDseq.",TRAIT,".5e-08.chr",chr,".filtered_by_QCsummaries.final.txt",sep=""),quote=F,row.names=F,col.names=T)
}
#END