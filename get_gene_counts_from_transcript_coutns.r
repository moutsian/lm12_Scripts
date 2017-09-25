#collapse Salmon's transcript counts to gene counts.
dataset="21121"
samplenum=12
for(samplenum in 1:16){
filename=paste("/lustre/scratch115/projects/paxgene/SalmonOutput_after_filtering/",dataset,".",samplenum,".postfilter.newidx.ISR.raw/quant.sf",sep="")

index=read.table("/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/transcript2gene.txt",head=T)
transdatafile=read.table(filename,head=T) #file with transcript quantifications
index_for_transdatafile=read.table("/lustre/scratch115/projects/paxgene/GRCh38_15_plus_hs38d1.known.fa.tlst") #this contains the matching between the transcripts in transdatafile and their ENST IDs.
transquant=cbind(index_for_transdatafile[,2],transdatafile)
fullinfo_file=merge(index,transquant,by.x="target_id",by.y="index_for_transdatafile[, 2]")
unique_genes=unique(index[,c(2:4)])
NumGeneReads=matrix(ncol=2,nrow=dim(unique_genes)[1],0)
colnames(NumGeneReads)=c("n_transcripts","TotalNumReads")
#for(i in 1:100){
for(i in 1:dim(unique_genes)[1]){
transidx=which(fullinfo_file$ens_gene==unique_genes$ens_gene[i])
if(length(transidx>0)){
NumGeneReads[i,]=c(length(transidx),sum(fullinfo_file$NumReads[transidx]))
}
}
GeneReads=cbind(unique_genes,NumGeneReads)
outputfile=paste(filename,".genecounts",sep="")
write.table(GeneReads,outputfile,quote=F,col.names=T,row.names=F,sep="\t")
}

#memory seems fine for parallelization and moving to cluster.
#END