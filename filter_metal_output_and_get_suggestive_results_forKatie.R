
#final filters:
# present in all 3 cohorts
# INFO > 0.4 in all 3 cohorts (warning INFO < 0.8)   
# I2 < 90 (warning I2 > 80)
# p_meta < min(p_cohort)
# MAF > 0.001 (warning MAF < 0.005)
# HWE p > 1e-7 in GWAS3, IBDSeq ctrls

trait="IBD"
HWE_IBDseq=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.IBDseq.ibd.outofHWE.controls.txt",head=T)
HWE_GWAS3=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS3.outofHWE.controls.txt",head=T)

for(chr in 22:1){
metares=read.table(paste("/lustre/scratch113/projects/crohns/iibdgc_meta/results/",tolower(trait),"/",chr,"-meta.txt",sep=""),head=T)

#individual results:
resGWAS3=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/",tolower(trait),"/",chr,".assoc",sep=""),head=T)
resIIBDGC=read.table(paste("/lustre/scratch113/projects/crohns/iibdgc_meta/data/IIBDGC/",tolower(trait),"/",chr,".assoc",sep=""),head=T)
resIBDseq=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/",tolower(trait),"/",chr,".assoc",sep=""),head=T)

#Katie has filtered for info and maf so just double check
if(min(resGWAS3$info)<0.4 |min(resIIBDGC$info)<0.4 | min(resIBDseq$info)<0.4 ){
print("error!! There are entries where info is below the threshold!")
}else{
print("All entries are filtered for info score")}

if(min(resGWAS3$info)<0.4 |min(resIIBDGC$info)<0.4 | min(resIBDseq$info)<0.4 ){
print("error!! There are entries where info is below the threshold!")
}else{
print("All entries are filtered for info score")}

#filter out variants due to low frequency
toolow_freq=c(as.character(resIIBDGC$rsid[which( as.numeric(resIIBDGC$all_maf)<0.001)]),as.character(resGWAS3$rsid[which(as.numeric(resGWAS3$all_maf)<0.001)]), as.character(resIBDseq$rsid[which(as.numeric(resIBDseq$all_maf)<0.001)]))
toolow_freq=unique(c(toolow_freq,as.character(resIIBDGC$rsid[which(as.numeric(resIIBDGC$all_maf)>0.999)]),as.character(resGWAS3$rsid[which(as.numeric(resGWAS3$all_maf)>0.999)]), as.character(resIBDseq$rsid[which(as.numeric(resIBDseq$all_maf) >0.999)]) ))
idx_toolowfreq=which(as.character(metares[,1])%in%toolow_freq)
if(length(idx_toolowfreq)>0){
metares_f0=metares[-idx_toolowfreq,]
print("Warning.  Entries are not fully filtered for MAF")
}
else{
metares_f0=metares
}

#filter out variants not present in all three datasets or with I2>90:
metares_f1=metares_f0[!grepl("?",metares_f0[,7],fixed=T) & as.numeric(metares_f0[,8])<90,]
#filter out variants not in HWE
outIBDseq=which(as.character(metares_f1[,1])%in%as.character(HWE_IBDseq[,2]))
outGWAS3=which(as.character(metares_f1[,1])%in%as.character(HWE_GWAS3[,2]))
outall=unique(c(outGWAS3,outIBDseq))
metares_f2=metares_f1[-outall,]

#finally, filter out variants which don't have a lower p_meta than the p in each separate cohort.
idxGWAS3=which(as.character(resGWAS3[,2])%in%as.character(metares_f2[,1]))
GWAS3filt=resGWAS3[idxGWAS3,c(2,4,42)]
GWAS3filt=GWAS3filt[!duplicated(GWAS3filt[,1]),]
GWAS3filto=GWAS3filt[order(GWAS3filt[,1]),]

idxIIBDGC=which(as.character(resIIBDGC[,1])%in%as.character(metares_f2[,1]))
IIBDGCfilt=resIIBDGC[idxIIBDGC,c(1,3,6)]
IIBDGCfilt=IIBDGCfilt[!duplicated(IIBDGCfilt[,1]),]
IIBDGCfilto=IIBDGCfilt[order(IIBDGCfilt[,1]),]

idxIBDseq=which(as.character(resIBDseq[,2])%in%as.character(metares_f2[,1]))
IBDseqfilt=resIBDseq[idxIBDseq,c(2,4,42)]
IBDseqfilt=IBDseqfilt[!duplicated(IBDseqfilt[,1]),]
IBDseqfilto=IBDseqfilt[order(IBDseqfilt[,1]),]

metares_f2o=metares_f2[order(metares_f2[,1]),]

if(sum(as.character(IBDseqfilto[,1])==as.character(metares_f2o[,1]))==dim(metares_f2o)[1] & sum(as.character(IIBDGCfilto[,1])==as.character(metares_f2o[,1]))==dim(metares_f2o)[1] & sum(as.character(GWAS3filto[,1])==as.character(metares_f2o[,1]))==dim(metares_f2o)[1]){
metares_f2o=cbind(metares_f2o,IBDseqfilto[,3],IIBDGCfilto[,3],GWAS3filto[,3],apply(cbind(IBDseqfilto[,3],IIBDGCfilto[,3],GWAS3filto[,3]),1,min))
}else{
print("something went wrong")
}

#lowpval=which(metares_f2o[,6]>metares_f2o[,15])
#metares_fin=metares_f2o[-lowpval,]
colnames(metares_f2o)[12:15]=c("Pval_IBDseq","Pval_IIBDGC","Pval_GWAS3","Min_single_cohort_pval")
write.table(metares_f2o,paste("/lustre/scratch113/projects/crohns/iibdgc_meta/results/",tolower(trait),"/",chr,"-meta.filtered.LM.without_pmeta_vs_pcohort.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}

#END