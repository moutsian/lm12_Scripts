##updated on Aprilt 4th to account for the bug with the last entries.


trait="UC"
dirin=paste("/lustre/scratch115/teams/anderson/ibd_conditional//datasetA/",trait,"/",sep="")
suggestive=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/IBDseq.",trait,".1e-05.filtered.txt",sep=""),head=T) #change this to contain the filtered version once you have it!
for(chr in 1:22){
myfile=paste(dirin,"step1.",trait,".",chr,".txt",sep="")
if(file.exists(myfile)){
datain=read.table(myfile,sep="\t",head=T)
if(dim(datain)[1]>1){
datao=as.matrix(datain[order(as.numeric(datain[,7])),])
datao=rbind(datao,datao[dim(datao)[1],])
dataout=NULL
n_lead_var=1;
for(i in 2:dim(datao)[1]){

if(as.numeric(datao[i,7])<=as.numeric(datao[i-1,8]) || abs(as.numeric(datao[i,3])-as.numeric(datao[i-1,3]))<500000){
datao[i,2]=as.matrix(unique(paste(as.character(datao[i,2]),as.character(datao[i-1,2]),sep=",")))
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

dataout=rbind(dataout,c(datao[i,],n_lead_var-1))
for(i in 1:dim(dataout)[1]){
dataout[i,2]=as.matrix(paste(unique(unlist(strsplit(dataout[i,2],","))),sep=",",collapse=","))
dataout[i,10]=as.matrix(paste(unique(unlist(strsplit(dataout[i,10],","))),sep=",",collapse=","))
dataout[i,9]=as.matrix(length(unique(unlist(strsplit(dataout[i,10],",")))))
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

#pvals=matrix(ncol=4,nrow=dim(dataout)[1],NA);
#for(i in 1:dim(dataout)[1]){
#idx=which(as.numeric(as.character(sugg[,3]))==as.numeric(as.character(dataout[i,3])))
#if(length(idx)==1){
#pvals[i,]=as.matrix(sugg[idx,14:17])
#}else{
#print("lathos!")
#}
#}
#colnames(pvals)=colnames(sugg)[14:17]

#dataout=cbind(dataout,LDwin_size,counts,pvals)
#colnames(dataout)=c(colnames(datain),"number_of_loci_merged","LD_window_size","n_vars_with_p_lt_1e05",colnames(pvals))

dataout=cbind(dataout,LDwin_size,counts)
colnames(dataout)=c(colnames(datain),"number_of_loci_merged","LD_window_size","n_vars_with_p_lt_1e05")

write.table(dataout,paste(dirin,"step2.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)
}else if(dim(datain)[1]==1){
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
idx=which(as.numeric(as.character(sugg[,3]))==as.numeric(as.character(dataout[3])))
#pvals=as.matrix(sugg[idx,14:17])
#names(pvals)=colnames(sugg)[14:17]
#dataout=c(dataout,pvals)
#names(dataout)=c(colnames(datain),"number_of_loci_merged","LD_window_size","n_vars_with_p_lt_1e05",names(pvals))
names(dataout)=c(colnames(datain),"number_of_loci_merged","LD_window_size","n_vars_with_p_lt_1e05")
write.table(dataout,paste(dirin,"step2.",trait,".",chr,".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=T)

}
}
}
#END
