#collapse Salmon's transcript counts to gene counts.
dataset="21121"
library="ISR"
#samplenum=12
for(samplenum in 1:16){
print(samplenum)
filename=paste("/lustre/scratch115/projects/paxgene/SalmonOutput/",library,"_",dataset,"/",dataset,".",samplenum,".new_index.",library,".multi.raw/quant.sf",sep="")

index=read.table("/lustre/scratch113/teams/anderson/users/jga/UsefulFiles/transcript2gene.txt",head=T)
transdatafile=read.table(filename,head=T) #file with transcript quantifications
index_for_transdatafile=read.table("/lustre/scratch115/projects/paxgene/GRCh38_15_plus_hs38d1.known.fa.tlst") #this contains the matching between the transcripts in transdatafile and their ENST IDs.
transquant=cbind(index_for_transdatafile[,2],transdatafile)

transquant_out=transquant[,c(1,3:6)]
colnames(transquant_out)[1]="Name"
write.table(transquant_out,filename,quote=F,col.names=T,row.names=F,sep="\t")
}
#END