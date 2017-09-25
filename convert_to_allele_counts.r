# convert_to_allele_counts
LOCUS="A"
POSTERIOR_THRESHOLD=0.5
DIR="/lustre/scratch115/projects/ibdgwas/HLA/preimpute_outputs/"
data=read.table(paste(DIR,"HLA_",LOCUS,"_hibag.imputations.txt",sep=""),head=T)
all_alleles=unique(c(as.character(data[,2]),as.character(data[,3])))
OUTPUT=matrix(nrow=dim(data)[1],ncol=length(all_alleles),0)
colnames(OUTPUT)=all_alleles
for(i in 1:dim(data)[1]){
if(data[i,4]>=POSTERIOR_THRESHOLD){
idx1=which(colnames(OUTPUT)==as.character(data[i,2]))
idx2=which(colnames(OUTPUT)==as.character(data[i,3]))
OUTPUT[i,idx1]=OUTPUT[i,idx1]+1
OUTPUT[i,idx2]=OUTPUT[i,idx2]+1
}else{
OUTPUT[i,]=NA
}
}

#now convert this to .gen format
GEN=t(OUTPUT)
GGEN=GEN
GGEN[which(GEN==1)]="0 1 0"
GGEN[which(GEN==0)]="1 0 0"
GGEN[which(GEN==2)]="0 0 1"
GGEN[which(is.na(GEN))]="0 0 0"
first_col= paste("HLA_",LOCUS,"_",1:dim(GGEN)[1],sep="")
second_col= paste("HLA_",LOCUS,"_",all_alleles,sep="")
options(scipen=999)
third_col=paste(100000*1:dim(GGEN)[1],sep="")
options(scipen=999)
GENFIN=cbind(first_col,second_col,third_col,"A","B",apply(GGEN,1,paste,collapse=" "))
write.table(GENFIN,paste(DIR,"HLA_",LOCUS,"_hibag.imputations_thres_",POSTERIOR_THRESHOLD,".gen",sep=""),quote=F,col.names=F,row.names=F)
#END