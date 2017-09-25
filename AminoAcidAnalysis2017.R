# Analysis with AminoAcids for the DRB*15:01-DQA*01:01 interaction
# March 2015 - update June 2017


##installation from the tar.gz file for 3.3.0. For example:
#install.packages("/software/team152/Rpackages/S4Vectors_0.6.6.tar.gz",repos=NULL,type="source",lib="/software/team152/Rpackages/")


DATADIR="/lustre/scratch115/projects/ibdgwas/HLA/new_reference_panel/"
OUTDIR="/lustre/scratch115/projects/ibdgwas/HLA/AA/"
LOCUS="DRB1"
AADRB=read.table(paste("/lustre/scratch115/projects/ibdgwas/HLA/AA/AA.HLA_",LOCUS,".2017.txt",sep=""))
HLA=read.table(paste(DATADIR,"HLA_",LOCUS,"_hibag.imputations.txt",sep=""),head=T,stringsAsFactors=FALSE)
library(reshape2,lib="/software/team152/Rpackages/")

POSTERIOR_THRESHOLD=0.7
#but first put all AAs together:
AADRB_FULL=do.call(paste, c(AADRB[2:dim(AADRB)[2]], sep=""))#nice use of do.call
DRBSEQ=colsplit(AADRB_FULL,pattern="",names=seq(1:nchar(AADRB_FULL[1]))) #now we have a table that can be used
# for association analysis
AAcnts=matrix(ncol=1,nrow=dim(DRBSEQ)[2],0) #number of different alleles
for(i in 1:length(AAcnts)){
AAcnts[i]=length(table(DRBSEQ[,i]))
}

# Also prepare the HLA alleles in the first column of AADRB so that they match the cond_unique_alleles column used
# in the column names of the standard HLA files I use. (E.g. "0103" instead of "DRB1*01:03") 
AAalleles=matrix(ncol=1,nrow=dim(AADRB)[1],NA);
for(i in 1:length(AAalleles)){
tmp=strsplit(strsplit(as.character(AADRB[i,1]),split='*',fixed=T)[[1]][2],split=":")[[1]][1:2]
AAalleles[i]=paste(tmp[1],tmp[2],sep=":");
}


mainAA=NULL #which aminoacids we have at the specific position
list_of_variants=NULL
mainAminoacids=NULL
allotherAminoacids=NULL
AA_COUNTS=as.matrix(HLA[,1])
for(pos in 1:dim(DRBSEQ)[2]){
if(AAcnts[pos]==1){
#do nothing, monomorphic AA position
}else{
if(AAcnts[pos]==2){
if(names(table(DRBSEQ[pos]))[1]!="*"){
mainAA=names(table(DRBSEQ[pos]))[1]
allotherAminoacids=c(allotherAminoacids,"*")
}else{
mainAA=names(table(DRBSEQ[pos]))[2]
allotherAminoacids=c(allotherAminoacids,names(table(DRBSEQ[pos]))[1])

}}else{
mainAA=names(table(DRBSEQ[pos]))
}
#now find which HLA alleles have the main AA in that position
for(m in 1:length(mainAA)){
list_of_variants=c(list_of_variants,paste(pos,mainAA[m],sep=""))
mainAminoacids=c(mainAminoacids,mainAA[m])
if(length(mainAA)>1){
allotherAminoacids=c(allotherAminoacids,paste(mainAA[-m],collapse=""))
}
mainHLA=AAalleles[which(DRBSEQ[pos]==mainAA[m])] #here are all alleles which have this AA
counter = matrix(ncol=1,nrow=dim(HLA)[1],0)
hap1=which(HLA[,2]%in%mainHLA)
hap2=which(HLA[,3]%in%mainHLA)
counter[hap1,1]= counter[hap1,1]+1
counter[hap2,1]= counter[hap2,1]+1
AA_COUNTS=cbind(AA_COUNTS,counter)
}
}
}
colnames(AA_COUNTS)=c("SampleID",list_of_variants)
AA_COUNTS[which(HLA[,4]<POSTERIOR_THRESHOLD),2:dim(AA_COUNTS)[2]]="0 0 0" #so that we don't count them. (note that we don't take posterior probs into account, we just use them to extract individuals with 
#low imputation probabilities and then we treat the rest as equals)
AA_COUNTS_GEN=AA_COUNTS[,2:dim(AA_COUNTS)[2]]
AA_COUNTS_GEN[which(AA_COUNTS_GEN=="0")]="1 0 0"
AA_COUNTS_GEN[which(AA_COUNTS_GEN=="1")]="0 1 0"
AA_COUNTS_GEN[which(AA_COUNTS_GEN=="2")]="0 0 1"
options(scipen=999)

first_col= paste("HLA_",LOCUS,"_",1:dim(AA_COUNTS_GEN)[2],sep="")
second_col= paste("HLA_",LOCUS,"_",list_of_variants,sep="")
third_col=paste(100000*1:dim(AA_COUNTS_GEN)[2],sep="")

GENFIN=cbind(first_col,second_col,third_col,mainAminoacids,allotherAminoacids,t(AA_COUNTS_GEN))
write.table(GENFIN,paste(OUTDIR,"HLA_",LOCUS,"_hibag.AA.imputations_thres_",POSTERIOR_THRESHOLD,".gen",sep=""),quote=F,col.names=F,row.names=F)


#now turn this to .gen file

#END
