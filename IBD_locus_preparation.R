# This R script will prepare the IBD loci and help to declare novel loci as per Jeff's suggestion in Jan 2016.
# It takes Javi's compiled list as input (tab-delimited). Note that at some point it's using a unix command to get some LD info.
# input and output files in dir: /lustre/scratch115/teams/anderson/ibd_conditional/

mydir="/lustre/scratch115/teams/anderson/ibd_conditional/"
lddir="/lustre/scratch113/resources/1000g/pairwise_ld/"
datain= as.matrix(read.table(paste(mydir,"IBD.CD.UC.GWAS.Hits.WithPositions.hg19.25Jan2016.LM.txt",sep=""),head=T,sep="\t"))

dataout=matrix(nrow=dim(datain)[1],ncol=dim(datain)[2]+6,NA)
dataout[,1:dim(datain)[2]]=datain

FROM=1
TO=200
if(TO>dim(datain)[1]){
TO=dim(datain)[1]
}

for( i in FROM:TO){
shellcommand=paste("zgrep ",datain[i,4]," ",lddir,"chr",as.numeric(datain[i,3]),".EUR.ld.bgz | awk -v VAR=",as.numeric(datain[i,4])," '{if($6>=0.6 && ( $4==VAR || $2==VAR )){print $0;}}'",sep="")
print(paste(i,": ",shellcommand,sep=""));
#now convert to matrix in Rvar
#left=min(c(as.numeric(ldmatrix[,2]),as.numeric(ldmatrix[,4])))
#right=max(c(as.numeric(ldmatrix[,2]),as.numeric(ldmatrix[,4])))
#dataout[i,(dim(datain)[2]+1):(dim(datain)[2]+4)]=as.matrix(c(left,right,left-250000,right+250000))


shellout=system(shellcommand,intern=T)
if(identical(shellout,character(0))){
#just 250kb either side
dataout[i,(dim(datain)[2]+1):(dim(datain)[2]+4)]=as.matrix(c(NA,NA,max(0,as.numeric(datain[i,4])-250000),as.numeric(datain[i,4])+250000),NA,NA)
}else{
ldmatrix=matrix(unlist(strsplit(shellout,"\t")),ncol=7,byrow=T) #make sure to check what happens when Rvar is zero
left=min(c(as.numeric(ldmatrix[,2]),as.numeric(ldmatrix[,4])))
right=max(c(as.numeric(ldmatrix[,2]),as.numeric(ldmatrix[,4])))
idxleft=which(as.numeric(ldmatrix[,2])==left)
if(identical(idxleft,integer(0))){
idxleft=which(as.numeric(ldmatrix[,4])==left)
SNPLEFT=ldmatrix[idxleft[1],5]
}else{
SNPLEFT=ldmatrix[idxleft[1],3]
}
idxright=which(as.numeric(ldmatrix[,2])==right)
SNPRIGHT=ldmatrix[idxright[1],5]
if(identical(idxright,integer(0))){
idxright=which(as.numeric(ldmatrix[,4])==right)
SNPRIGHT=ldmatrix[idxright[1],3]
}
dataout[i,(dim(datain)[2]+1):(dim(datain)[2]+6)]=as.matrix(c(left,right,left-250000,right+250000,SNPLEFT,SNPRIGHT))
}

}

colnames(dataout)=c("PMID","Trait","Chr","Position","SNP_ID","Freq","Pval","LeftMostPos_in_LD","RightMostPos_in_LD","WindowLeft","WindowRight","LeftMostSNP_in_LD","RightMostSNP_in_LD")
write.table(dataout[FROM:TO,],paste(mydir,"loci_list_before_any_merging_r2_0.60.",FROM,"_",TO,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)

#END