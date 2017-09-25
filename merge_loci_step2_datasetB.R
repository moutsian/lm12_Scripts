
trait="UC"
dirin=paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",trait,"/",sep="")
suggestive=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",trait,".1e_05.filtered.lowfreq.HWE.PASSED.txt",sep=""),head=F) 
# UK_only_analysis_B.CD.1e_05.filtered.HWE.PASSED.txt
for(chr in 1:22){
myfile=paste(dirin,"step1.filtered.",trait,".",chr,".txt",sep="")
if(file.exists(myfile)){
datain=read.table(myfile,sep="\t",head=T)
if(dim(datain)[1]>1){
datao=as.matrix(datain[order(as.numeric(datain[,7])),])
dataout=NULL
n_lead_var=1;
for(i in 2:dim(datao)[1]){

if(as.numeric(datao[i,7])<=as.numeric(datao[i-1,8]) || abs(as.numeric(datao[i,3])-as.numeric(datao[i-1,3]))<500000){
datao[i,2]=as.matrix(paste(as.character(datao[i,2]),as.character(datao[i-1,2]),sep=","))
datao[i,9]=as.numeric(datao[i,9])+as.numeric(datao[i-1,9])
datao[i,10]=paste(datao[i,10],datao[i-1,10],sep=",")
n_lead_var=n_lead_var+1;
if(as.numeric(datao[i,6])> as.numeric(datao[i-1,6]) ){
datao[i,3:6]=datao[i-1,3:6]
}
datao[i,7]=min(as.numeric(datao[i,7]),as.numeric(datao[i-1,7]))
datao[i,8]=max(as.numeric(datao[i,8]),as.numeric(datao[i-1,8]))
}else{
dataout=rbind(dataout,c(datao[i-1,],n_lead_var))
n_lead_var=1; #reset
}
}

if(n_lead_var>1){
dataout=rbind(dataout,c(datao[i,],n_lead_var))
n_lead_var=1; #reset
}
#outputfile
#add as final column the number of SNPs with suggestive evidence for association in each of the regions
sugg=suggestive[suggestive[,1]==chr,]
counts=matrix(ncol=1,nrow=dim(dataout)[1],0)
for(j in 1:dim(sugg)[1]){
idx=which(as.numeric(dataout[,8])>=as.numeric(as.character(sugg[j,3])) & as.numeric(dataout[,7])<=as.numeric(as.character(sugg[j,3])));
if(length(idx)>0){
counts[idx]=counts[idx]+1;
}
}
LDwin_size=as.numeric(dataout[,8])-as.numeric(dataout[,7])
dataout=cbind(dataout,LDwin_size,counts)
colnames(dataout)=c(colnames(datain),"number_of_loci_merged","LD_window_size","n_vars_with_p_lt_1e05_incl.lowfreq")
write.table(dataout,paste(dirin,"step2.filtered.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)
}else if(dim(datain)[1]==1){
#add as final column the number of SNPs with suggestive evidence for association in each of the regions
sugg=suggestive[suggestive[,1]==chr,]
counts=0;
for(j in 1:dim(sugg)[1]){
idx=which(as.numeric(datain[1,8])>=as.numeric(as.character(sugg[j,3])) & as.numeric(datain[1,7])<=as.numeric(as.character(sugg[j,3])));
if(length(idx)>0){
counts=counts+1;
}
}
LDwin_size=as.numeric(datain[,8])-as.numeric(datain[,7])
dataout=c(datain,1,LDwin_size[1],counts)
names(dataout)=c(colnames(datain),"number_of_loci_merged","LD_window_size","n_vars_with_p_lt_1e05")
write.table(dataout,paste(dirin,"step2.filtered.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)

}
}
}
#END
