
#this R file merges the step2 files across traits, and checks whether the hits we find in our study are previously reported or not. If they are not, it adds I2 and frequency info and outputs them separately
trait="UC"
n_alldat=3
if(toupper(trait)=="IBD"){
n_alldat=4
}
previous=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/IBD.CD.UC.GWAS.Hits.WithPositions.hg19.25Jan2016.LM.txt",head=T,sep="\t")
shellcommand=paste(" cat /lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",toupper(trait),"/step2.filtered.",toupper(trait),"*.txt >  ~lm12/IBD_conditional/step2.filtered.",toupper(trait),".datasetB.all",sep="")
system(shellcommand)
hitsB=as.matrix(read.table(paste("~lm12/IBD_conditional/step2.filtered.",toupper(trait),".datasetB.all",sep=""),sep="\t",head=T))

previousreported=previous[previous[,1]!="Novel",]

toinclude=which(previousreported[,7]<5e-08)
#toinclude=which(as.character(previous[,1])=="26192919"|as.character(previous[,1])=="23128233"|as.character(previous[,1])=="21102463"|as.character(previous[,1])=="21297633"|as.character(previous[,1])=="FM")
previousreported=previousreported[toinclude,]

previousgw=previousreported[as.numeric(as.character(previousreported[,7]))<=5e-08,]
hitsBplus=matrix(ncol=dim(hitsB)[2]+4,nrow=dim(hitsB)[1],"")
hitsBplus[,1:dim(hitsB)[2]]=hitsB
colnames(hitsBplus)=c(colnames(hitsB),"pos_reported_hits","p_value_reported_hits","PMID","traits_reported")
for(i in 1:dim(previousgw)[1]){
idx=which(as.numeric(previousgw[i,3])==as.numeric(hitsBplus[,1]) & as.numeric(previousgw[i,4])>=as.numeric(hitsBplus[,7]) & as.numeric(previousgw[i,4])<=as.numeric(hitsBplus[,8]))
if(length(idx)>0){
hitsBplus[idx[1],dim(hitsB)[2]+1]=paste(hitsBplus[idx[1],dim(hitsB)[2]+1],as.character(previousgw[i,4]),sep=",")
hitsBplus[idx[1],dim(hitsB)[2]+2]=paste(hitsBplus[idx[1],dim(hitsB)[2]+2],as.character(previousgw[i,7]),sep=",")
hitsBplus[idx[1],dim(hitsB)[2]+3]=paste(hitsBplus[idx[1],dim(hitsB)[2]+3],as.character(previousgw[i,1]),sep=",")
hitsBplus[idx[1],dim(hitsB)[2]+4]=paste(hitsBplus[idx[1],dim(hitsB)[2]+4],as.character(previousgw[i,2]),sep=",")
}
}
hitsBplus_all=hitsBplus[which(hitsBplus[,1]!="CHR"),]
write.table(hitsBplus_all,paste("~lm12/IBD_conditional/step2.",trait,".all.with_info_on_previous_hits_from_all_PAPERS.datasetB.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#now output novel loci, with I2, freq_1KG and frequencies in our samples
novel=hitsBplus[which(hitsBplus[,1]!="CHR" & hitsBplus[,15]==""),1:dim(hitsB)[2]]
novel_plus=matrix(nrow=dim(novel)[1],ncol=dim(novel)[2]+5,"")
novel_plus[,1:dim(novel)[2]]=as.matrix(novel)
colnames(novel_plus)=c(colnames(novel),"MAF_GWAS1","MAF_GWAS2","MAF_GWAS3","MAF_IBDseq","FREQ_in_1KG_CEU_GBR","I2","datasets_used")
for(i in 1:dim(novel)[1]){
shellcommand=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_freq_and_I2_datasetB.sh ",toupper(trait),as.character(novel[i,1]),as.character(novel[i,3]),sep=" ")
shellout=unlist(strsplit(system(shellcommand,intern=T)," "))
shellfreqin=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_our_freq_datasetB.sh ",toupper(trait),as.character(novel[i,1]),as.character(novel[i,3]),sep=" ")
shellfreqout=unlist(strsplit(system(shellfreqin,intern=T),","))
if(length(shellfreqout)==4){
novel_plus[i,dim(novel)[2]+1]=unlist(strsplit(shellfreqout[2],":"))[2]
novel_plus[i,dim(novel)[2]+2]=unlist(strsplit(shellfreqout[3],":"))[2]
novel_plus[i,dim(novel)[2]+3]=unlist(strsplit(shellfreqout[4],":"))[2]
}else{
print ("Something went wrong when getting frequencies from our data")
}

if(length(shellout==9)){
data_plus[i,dim(data)[2]+1]=shellout[7]
data_plus[i,dim(data)[2]+2]=shellout[8]
data_plus[i,dim(data)[2]+3]=n_alldat - as.numeric(shellout[9])
}else{
print("Something went wrong!")
}
}

Trait=rep(trait,dim(novel_plus)[1])
novel_plus=cbind(Trait,novel_plus)
write.table(novel_plus,paste("~lm12/IBD_conditional/step2.",trait,".NOVEL.with_info_on_previous_hits_from_all_PAPERS.datasetB.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
hitsBplus[10,c(1,3:9,11,12:16)]
#END