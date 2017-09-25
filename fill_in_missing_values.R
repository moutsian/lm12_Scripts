
trait="IBD"
dirin=paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/",trait,"/",sep="")
dirout=paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",trait,"/",sep="")
dataout=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",trait,".5e-08.txt",sep=""),sep=" ",head=F)
#chr="1"
for(chr in 1:22){
datain=read.table(paste(dirin,"results.",trait,".",chr,".txt",sep=""),sep=" ",head=T)
datafilt=dataout[dataout[,1]==chr,1:8]
colnames(datafilt)=colnames(datain)
datafilt[,7]=NA
datafilt[,8]=NA

for(i in 1:dim(datafilt)[1]){
idx=which(as.numeric(datain[,3])==as.numeric(datafilt[i,3]))
if(length(idx)>0){
datafilt[i,7]=as.numeric(datain[idx[1],7])
datafilt[i,8]=as.numeric(datain[idx[1],8])
}
}

#outputfile
write.table(datafilt,paste(dirout,"results.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)
}

#END
