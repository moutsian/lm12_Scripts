
#Filters for dat A:
# INFO > 0.8
# MAF < 0.005
# HWE p > 1e-7 in IBDSeq ctrls

trait="IBD"

HWE_IBDseq=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.IBDseq.ibd.outofHWE.controls.txt",head=T)

for(chr in 22:1){
print(paste("chrom:",chr))
resIBDseq=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/",toupper(trait),"/IBDseq.",toupper(trait),".chr",chr,".1e-05.txt",sep=""),head=T)
to_keep=matrix(ncol=1,nrow=dim(resIBDseq)[1],0)

#filter for info
to_keep[which(resIBDseq$info>0.8)]=to_keep[which(resIBDseq$info>0.8)]+1
#filter out variants due to high frequency

low_freq=unique(c(as.character(resIBDseq$rsid[which(as.numeric(resIBDseq$all_maf)>0.001 & as.numeric(resIBDseq$all_maf)<0.999 )])) )
idx_lowfreq=which(as.character(resIBDseq[,2])%in%low_freq)
to_keep[idx_lowfreq]=to_keep[idx_lowfreq]+1

#filter out variants not in HWE
outIBDseq=which(as.character(resIBDseq$rsid)%in%as.character(HWE_IBDseq[,2]))
to_keep[outIBDseq]=to_keep[outIBDseq]-1

#now keep only dat A entries:
IBDseqA=resIBDseq[which(to_keep==2),c(1:2,4:6,9,29:31,42:45)]

write.table(IBDseqA,paste("/lustre/scratch115/teams/anderson/ibd_conditional/datasetA/",toupper(trait),"/IBDseq.",toupper(trait),".chr",chr,".1e-05.filtered.may.info0.8.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}

#END
