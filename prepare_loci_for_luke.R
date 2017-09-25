##

knownonly=1;


IBDf=read.table("~lm12/IBD_conditional/IBD.all.round6",head=T,sep="\t")
CDf=read.table("~lm12/IBD_conditional/CD.all.round6",head=T,sep="\t")
UCf=read.table("~lm12/IBD_conditional/UC.all.round6",head=T,sep="\t")

ALL=rbind(as.matrix(UCf),as.matrix(CDf),as.matrix(IBDf))
ALL=ALL[ALL[,2]!="CHR",]
fileout="known_and_novel"
if(knownonly==1){
fileout="known_only"
ALL=ALL[ALL[,20]!="",]

}
datao=as.matrix(ALL[order(as.numeric(as.character(ALL[,2])),as.numeric(as.character(ALL[,4]))),])
#datao=datao[datao[,1]!="CHR",]
dataout=NULL
n_lead_var=1;
#alltraits=matrix(ncol=2,nrow=dim(datao)[1],"")
#alltraits[,1]=as.character(datao[,14])
#alltraits[,2]=as.character(datao[,7])
#colnames(alltraits)=c("WhichTraits","Pvalues")
#datao=cbind(datao,alltraits)
for(i in 2:dim(datao)[1]){

if( as.numeric(datao[i,2])==as.numeric(datao[i-1,2]) & (as.numeric(datao[i,8])<=as.numeric(datao[i-1,9]) || abs(as.numeric(datao[i,4])-as.numeric(datao[i-1,4]))<500000) ){
datao[i,3]=as.matrix(paste(as.character(datao[i,3]),as.character(datao[i-1,3]),sep=","))
datao[i,20]=as.matrix(paste(as.character(datao[i,20]),as.character(datao[i-1,20]),sep=","))
datao[i,21]=as.matrix(paste(as.character(datao[i,21]),as.character(datao[i-1,21]),sep=","))

datao[i,10]=as.numeric(datao[i,10])+as.numeric(datao[i-1,10])
datao[i,11]=paste(datao[i,11],datao[i-1,11],sep=",")
datao[i,14]=as.numeric(datao[i,14])+as.numeric(datao[i-1,14])
n_lead_var=n_lead_var+1;
if(as.numeric(datao[i,7])> as.numeric(datao[i-1,7]) ){
datao[i,1:7]=datao[i-1,1:7]
}
datao[i,8]=min(as.numeric(datao[i,8]),as.numeric(datao[i-1,8]))
datao[i,9]=max(as.numeric(datao[i,9]),as.numeric(datao[i-1,9]))
}else{
dataout=rbind(dataout,c(datao[i-1,],n_lead_var))
n_lead_var=1; #reset
}
}

dataout=rbind(dataout,c(datao[i,],n_lead_var))

colnames(dataout)[22]=c("How_many_of_IBD_UC_CD_have_this_at_gw_significance")

write.table(dataout,paste("~lm12/IBD_conditional/list_of_loci_for_Luke.round6.",fileout,".txt",sep=""),col.names=T,quote=F,row.names=F,sep="\t")



#END
