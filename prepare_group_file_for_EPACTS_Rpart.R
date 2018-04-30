#prepare_group_file_for_EPACTS_Rpart.R
chrom=21
for(chrom in 1:22){
print(paste("chrom:",chrom,sep=""))
varfile=paste("tili.i3.site_QCplus_binom.sample_QC.annot.vcf.",chrom,".variants",sep="")
annofile=paste("tili.i3.site_QC.sample_QC.vcf.annot.",chrom,".annot.fc.tmp",sep="")
outfile=paste("tili.i3.site_QC.sample_QC.",chrom,".annot.fc.genes.epacts",sep="")
var_table=read.table(varfile,stringsAsFactors=F)
anno_table=read.table(annofile,stringsAsFactors=F)
together=merge(anno_table,var_table,by.x="V1",by.y="V3")
# We can do this both for genes and for transcripts. I will do it for genes first.
# Note that I am not applying any MAF threshold here at present.
genes=unique(together[,4])
epacts_table=matrix(ncol=2,nrow=length(genes),"")
epacts_table[,1]=genes
for(i in 1:dim(together)[1]){
if(i%%1000==0){print(paste(i,"out of",dim(together)[1]))}
idx=which(genes==together[i,4])
tmp=epacts_table[idx,2]
epacts_table[idx,2]=paste(tmp," ",together[i,2],"_",together[i,17],"/",together[i,18],sep="")
}
#since we have entries by transcript there will be multiple entries per gene. Thus, unique and sort after you are done
test=lapply(epacts_table[,2],function(x) paste(sort(unique(unlist(strsplit((x),split=" ")))),collapse=' '))
epacts_table=cbind(genes,test)
write.table(epacts_table,outfile,col.names=F,row.names=F,quote=F,sep="\t")
}
#END