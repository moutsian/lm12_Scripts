#updated to work with METAL output.

trait=tolower("ibd")
TRAIT=toupper(trait)
new=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".5e_08.filtered.txt",sep=""),head=T)
#final_meta_analysis_C.IBD.5e_08.filtered.txt
for(chr in 1:22){
datasub=new[as.numeric(as.character(new[,1]))==chr,]
old=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results_beforePCcorrection/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T)
oldB=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results_beforePCcorrection/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T)
LD=matrix(ncol=8,nrow=dim(datasub)[1],NA)
if(dim(datasub)[1]>0){
colnames(LD)=colnames(old)
LD[,1]=chr
LD[,c(2:6)]=as.matrix(datasub[,c(3,2,4,5,8)])
for(i in 1:dim(LD)[1]){
idx=which(as.numeric(as.character(old[,3]))==as.numeric(as.character(LD[i,3])))
if(length(idx)>0){
LD[i,c(7,8)]=as.matrix(old[idx[1],c(7,8)])
}else{
idx=which(as.character(oldB[,3])==as.character(LD[i,3]))
if(length(idx)>0){
LD[i,c(7,8)]=as.matrix(oldB[idx[1],c(7,8)])
}
#else{
#idx=which(as.character(oldC[,3])==as.character(LD[i,3]))
#if(length(idx)>0){
#LD[i,c(7,8)]=as.matrix(oldC[idx[1],c(7,8)])
#}
#}
}
}
write.table(LD,paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),col.names=T,row.names=F,quote=F)
}
}
#END