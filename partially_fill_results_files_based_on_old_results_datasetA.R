#updated to work with METAL output.

trait=tolower("uc")
TRAIT=toupper(trait)
new=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/IBDseq.",TRAIT,".5e-08.filtered_by_QCsummaries.final.txt",sep=""),head=F)
#final_meta_analysis_C.IBD.5e_08.filtered.txt
for(chr in 1:22){
datasub=new[as.numeric(as.character(new[,1]))==chr,]
Cres=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T)
old_B_res=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T)
LD=matrix(ncol=8,nrow=dim(datasub)[1],NA)
if(dim(datasub)[1]>0){
colnames(LD)=colnames(Cres)
LD[,1]=chr
LD[,c(2:6)]=as.matrix(datasub[,c(2,3,4,5,10)])
for(i in 1:dim(LD)[1]){
idx=which(as.numeric(as.character(Cres[,3]))==as.numeric(as.character(LD[i,3])))
if(length(idx)>0){
LD[i,c(7,8)]=as.matrix(Cres[idx[1],c(7,8)])
}else{
idx=which(as.character(old_B_res[,3])==as.character(LD[i,3]))
if(length(idx)>0){
LD[i,c(7,8)]=as.matrix(old_B_res[idx[1],c(7,8)])
}
#else{
#idx=which(as.character(oldC[,3])==as.character(LD[i,3]))
#if(length(idx)>0){
#LD[i,c(7,8)]=as.matrix(oldC[idx[1],c(7,8)])
#}
#}
}
}
write.table(LD,paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),col.names=T,row.names=F,quote=F)
}
}
#END