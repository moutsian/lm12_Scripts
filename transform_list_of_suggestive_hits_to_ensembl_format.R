
trait="IBD"
utrait=toupper(trait)
data=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.txt",sep=""),head=T)
strand=rep("+",dim(data)[1]);
alleles=matrix(ncol=1,nrow=dim(data)[1],NA)
for(i in 1:dim(data)[1]){
alleles[i]=paste(as.character(data[i,4]),"/",as.character(data[i,5]),sep="")
}
ensembl=cbind(as.character(data[,1]),as.character(data[,3]),as.character(data[,3]),alleles,strand)
write.table(ensembl,paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.txt.ensembl",sep=""),quote=F,row.names=F,col.names=F)

#END