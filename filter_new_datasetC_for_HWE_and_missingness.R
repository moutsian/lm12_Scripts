
#results are already filtered for MAF 0.5%, INFO>0.4 and I2<80.
trait="IBD"
trait=tolower(trait)
TRAIT=toupper(trait)
data=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".5e_08.filtered.txt",sep=""),head=T)
#IBD=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.IBD.5e-08.filtered.txt",head=F)
#UC=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.UC.5e-08.filtered.txt",head=F)
GWAS3outofHWE=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS3.outofHWE.controls.txt",head=T)
IBDseqoutofHWE=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.IBDseq.ibd.outofHWE.controls.txt",head=T)

which(as.character(data[,1])==as.character(GWAS3outofHWE[,1]) && as.character(data[,3])%in%as.character(GWAS3outofHWE[,4]))
which(as.character(data[,2])%in%as.character(GWAS3outofHWE[,4]))
which(as.character(data[,2])%in%as.character(IBDseqoutofHWE[,4]))
datafilt=data[-which(as.character(data[,2])%in%as.character(IBDseqoutofHWE[,4])),]
write.table(datafilt,paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".5e_08.filtered.HWEpass.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

to_include=matrix(ncol=1,nrow=dim(data)[1],1)


#END