
#this R file checks whether the hits we find in our study are previously reported or not. If they are not, it adds I2 and frequency info and outputs them separately
trait="UC"
previous=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/IBD.CD.UC.GWAS.Hits.WithPositions.hg19.25Jan2016.LM.txt",head=T,sep="\t")
#hitsC=as.matrix(read.table(paste("~lm12/IBD_conditional/step2.",trait,".all.round2",sep=""),sep="\t",head=T))
hitsC=as.matrix(read.table(paste("~lm12/IBD_conditional/step2.filtered.",trait,".all",sep=""),sep="\t",head=T))

#previousreported=previous[previous[,1]!="Novel",]
toinclude=which(as.character(previous[,1])=="26192919"|as.character(previous[,1])=="23128233"|as.character(previous[,1])=="21102463"|as.character(previous[,1])=="21297633"|as.character(previous[,1])=="FM")
previousreported=previous[toinclude,]

previousgw=previousreported[as.numeric(as.character(previousreported[,7]))<=5e-08,]
hitsCplus=matrix(ncol=dim(hitsC)[2]+4,nrow=dim(hitsC)[1],"")
hitsCplus[,1:dim(hitsC)[2]]=hitsC
colnames(hitsCplus)=c(colnames(hitsC),"pos_reported_hits","p_value_reported_hits","PMID","traits_reported")
for(i in 1:dim(previousgw)[1]){
idx=which(as.numeric(previousgw[i,3])==as.numeric(hitsCplus[,1]) & as.numeric(previousgw[i,4])>=as.numeric(hitsCplus[,7]) & as.numeric(previousgw[i,4])<=as.numeric(hitsCplus[,8]))
if(length(idx)>0){
hitsCplus[idx[1],13]=paste(hitsCplus[idx[1],13],as.character(previousgw[i,4]),sep=",")
hitsCplus[idx[1],14]=paste(hitsCplus[idx[1],14],as.character(previousgw[i,7]),sep=",")
hitsCplus[idx[1],15]=paste(hitsCplus[idx[1],15],as.character(previousgw[i,1]),sep=",")
hitsCplus[idx[1],16]=paste(hitsCplus[idx[1],16],as.character(previousgw[i,2]),sep=",")
}
}
write.table(hitsCplus,paste("~lm12/IBD_conditional/step2.",trait,".all.round2.with_info_on_previous_hits_from_5PAPERS.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#now output novel loci, with I2, freq_1KG and frequencies in our samples
novel=hitsCplus[which(hitsCplus[,1]!="CHR" & hitsCplus[,15]==""),1:12]
novel_plus=matrix(nrow=dim(novel)[1],ncol=dim(novel)[2]+5,"")
novel_plus[,1:dim(novel)[2]]=as.matrix(novel)
colnames(novel_plus)=c(colnames(novel),"MAF_IIBDGC","MAF_GWAS3","MAF_IBDseq","FREQ_in_1KG_CEU_GBR","I2")
for(i in 1:dim(novel)[1]){
shellcommand=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_freq_and_I2.sh ",toupper(trait),as.character(novel[i,1]),as.character(novel[i,3]),sep=" ")
shellout=unlist(strsplit(system(shellcommand,intern=T)," "))
shellfreqin=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_our_freq.sh ",toupper(trait),as.character(novel[i,1]),as.character(novel[i,3]),sep=" ")
shellfreqout=unlist(strsplit(system(shellfreqin,intern=T),","))
if(length(shellfreqout)==4){
novel_plus[i,dim(novel)[2]+1]=unlist(strsplit(shellfreqout[2],":"))[2]
novel_plus[i,dim(novel)[2]+2]=unlist(strsplit(shellfreqout[3],":"))[2]
novel_plus[i,dim(novel)[2]+3]=unlist(strsplit(shellfreqout[4],":"))[2]
}else{
print ("Something went wrong when getting frequencies from our data")
}

if(length(shellout==8)){
novel_plus[i,dim(novel)[2]+4]=shellout[7]
novel_plus[i,dim(novel)[2]+5]=shellout[8]
}else{
print("Something went wrong!")
}
}

Trait=rep(trait,dim(novel_plus)[1])
novel_plus=cbind(Trait,novel_plus)
write.table(novel_plus,paste("~lm12/IBD_conditional/step2.",trait,".NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
hitsCplus[10,c(1,3:9,11,12:16)]
#END