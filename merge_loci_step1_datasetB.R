
trait="IBD"
dirin=paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",trait,"/",sep="") 
#chr="1"
for(chr in 1:22){
if(file.exists(paste(dirin,"results.",trait,".",chr,".txt.filtered",sep=""))){
datain=read.table(paste(dirin,"results.",trait,".",chr,".txt.filtered",sep=""),sep=" ",head=T) #work on the .filtered data --ideally post-HWE after results are here !

if(dim(datain)[1]>0){
datao=datain[order(datain[,6]),]
dataout=matrix(ncol=dim(datao)[2]+2,nrow=1,1)

dataout[1,]=c(as.matrix(datao[1,]),0,as.character(datao[1,2]))
if(is.na(dataout[1,7])){dataout[1,7]=as.numeric(dataout[1,3])}
if(is.na(dataout[1,8])){dataout[1,8]=as.numeric(dataout[1,3])}

for(i in 2:dim(datao)[1]){
idx=which(as.numeric(datao[i,3])>=as.numeric(dataout[,7]) & as.numeric(datao[i,3])<=as.numeric(dataout[,8]))
if(length(idx)>=1){
#then add it to that locus. Else keep it separate
dataout[idx[1],9]=as.numeric(dataout[idx[1],9])+1; #number of gw signals in that locus
dataout[idx[1],10]=paste(dataout[idx[1],10],datao[i,2],sep=",")
}else{
dataout=rbind(dataout,c(as.matrix(datao[i,]),1,as.character(datao[i,2])))
if(is.na(dataout[dim(dataout)[1],7])){dataout[dim(dataout)[1],7]=as.numeric(dataout[dim(dataout)[1],3])}
if(is.na(dataout[dim(dataout)[1],8])){dataout[dim(dataout)[1],8]=as.numeric(dataout[dim(dataout)[1],3])}
}
}
colnames(dataout)=c(colnames(datain),"Total_GW_significant_variants_in_region","GW_variants_in_region")
#outputfile
write.table(dataout,paste(dirin,"step1.filtered.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)
}}
}

#END
