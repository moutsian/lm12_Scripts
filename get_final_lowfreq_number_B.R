
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
HWE_GWAS2=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS2.outofHWE.controls.txt",head=T)
HWE_GWAS1=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS1.outofHWE.controls.txt",head=T)

for(chr in 22:1){
print(paste("chrom:",chr))
metares=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/UKIBDGC_meta_B/results/",tolower(trait),"/",chr,"-meta.txt",sep=""),head=T)

#individual results:
resGWAS3=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/",tolower(trait),"/",chr,".assoc",sep=""),head=T)
resIBDseq=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq/",tolower(trait),"/",chr,".assoc",sep=""),head=T)
if(tolower(trait)!="cd"){
resGWAS2=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS2/uc/",chr,".assoc",sep=""),head=T)
}else{
resGWAS2=resGWAS3}
if(tolower(trait)!="uc"){
resGWAS1=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS1/cd/",chr,".assoc",sep=""),head=T)
}else{
resGWAS1=resGWAS3
}

#Katie has filtered for info and maf so just double check
if( min(resGWAS1$info)<0.4 | min(resGWAS2$info)<0.4 | min(resGWAS3$info)<0.4 | min(resIBDseq$info)<0.4 ){
print("error!! There are entries where info is below the threshold!")
}else{
print("All entries are filtered for info score")}


#filter out variants due to low frequency
to_keepfreq=c(as.character(resGWAS1$rsid[which( as.numeric(resGWAS1$all_maf)>=0.001 && as.numeric(resGWAS1$all_maf)<0.05)]),as.character(resGWAS2$rsid[which( as.numeric(resGWAS2$all_maf)>=0.001 && as.numeric(resGWAS2$all_maf)<0.05)]),as.character(resGWAS3$rsid[which(as.numeric(resGWAS3$all_maf)>=0.001 && as.numeric(resGWAS3$all_maf)<0.05)]), as.character(resIBDseq$rsid[which(as.numeric(resIBDseq$all_maf)>=0.001 && as.numeric(resIBDseq$all_maf)<0.05)]))
idx_tokeepfreq=which(as.character(metares[,1])%in%to_keepfreq)
#toolow_freq=c(as.character(resGWAS1$rsid[which( as.numeric(resGWAS1$all_maf)<0.001)]),as.character(resGWAS2$rsid[which( as.numeric(resGWAS2$all_maf)<0.001)]),as.character(resGWAS3$rsid[which(as.numeric(resGWAS3$all_maf)<0.001)]), as.character(resIBDseq$rsid[which(as.numeric(resIBDseq$all_maf)<0.001)]))
#toolow_freq=c(toolow_freq,as.character(resGWAS1$rsid[which( as.numeric(resGWAS1$all_maf)>0.05)]),as.character(resGWAS2$rsid[which( as.numeric(resGWAS2$all_maf)<0.05)]),as.character(resGWAS3$rsid[which(as.numeric(resGWAS3$all_maf)<0.001)]), as.character(resIBDseq$rsid[which(as.numeric(resIBDseq$all_maf)<0.001)]))
#toolow_freq=unique(c(toolow_freq,as.character(resGWAS1$rsid[which(as.numeric(resGWAS1$all_maf)>0.999)]),as.character(resGWAS2$rsid[which(as.numeric(resGWAS2$all_maf)>0.999)]),as.character(resGWAS3$rsid[which(as.numeric(resGWAS3$all_maf)>0.999)]), as.character(resIBDseq$rsid[which(as.numeric(resIBDseq$all_maf) >0.999)]) ))
#idx_toolowfreq=which(as.character(metares[,1])%in%toolow_freq)
if(length(idx_tokeepfreq)>0){
#metares_f0=metares[-idx_toolowfreq,]
metares_f0=metares[idx_tokeepfreq,]
print("Warning.  Entries are not fully filtered for MAF")
}
else{
metares_f0=metares
}

#filter out variants not present in all three datasets 
metares_f1=metares_f0[!grepl("?",metares_f0[,7],fixed=T) ,]
#filter out variants not in HWE
#outIBDseq=which(as.character(metares_f1[,1])%in%as.character(HWE_IBDseq[,2]))
#outGWAS3=which(as.character(metares_f1[,1])%in%as.character(HWE_GWAS3[,2]))
#outGWAS1=which(as.character(metares_f1[,1])%in%as.character(HWE_GWAS1[,2]))
#outGWAS2=which(as.character(metares_f1[,1])%in%as.character(HWE_GWAS2[,2]))
#outall=unique(c(outGWAS1,outGWAS2,outGWAS3,outIBDseq))
#metares_f2=metares_f1[-outall,]

#metares_f2o=metares_f2[order(metares_f2[,1]),]

if(sum(as.character(IBDseqfilto[,1])==as.character(metares_f2o[,1]))==dim(metares_f2o)[1] & sum(as.character(GWAS1filto[,1])==as.character(metares_f2o[,1]))==dim(metares_f2o)[1] & sum(as.character(GWAS2filto[,1])==as.character(metares_f2o[,1]))==dim(metares_f2o)[1] & sum(as.character(GWAS3filto[,1])==as.character(metares_f2o[,1]))==dim(metares_f2o)[1]){
	if(tolower(trait)=="uc"){
	metares_f2o=cbind(metares_f2o,rep(NA,dim(GWAS1filto)[1]),GWAS2filto[,3],GWAS3filto[,3],IBDseqfilto[,3],apply(cbind(GWAS2filto[,3],GWAS3filto[,3],IBDseqfilto[,3]),1,min))
	}else if(tolower(trait)=="cd"){
	metares_f2o=cbind(metares_f2o,GWAS1filto[,3],rep(NA,dim(GWAS2filto)[1]),GWAS3filto[,3],IBDseqfilto[,3],apply(cbind(GWAS1filto[,3],GWAS3filto[,3],IBDseqfilto[,3]),1,min))
	}else if(tolower(trait)=="ibd"){
	metares_f2o=cbind(metares_f2o,GWAS1filto[,3],GWAS2filto[,3],GWAS3filto[,3],IBDseqfilto[,3],apply(cbind(GWAS1filto[,3],GWAS2filto[,3],GWAS3filto[,3],IBDseqfilto[,3]),1,min))
	}
}else{
print("something went wrong")
}

lowpval=which(metares_f2o[,6]>metares_f2o[,16])
metares_fin=metares_f2o[-lowpval,]
colnames(metares_fin)[12:16]=c("Pval_GWAS1","Pval_GWAS2","Pval_GWAS3","Pval_IBDseq","Min_single_cohort_pval")
write.table(metares_fin,paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",toupper(trait),"/",chr,"-meta.filtered.lowfreq.LM.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}

#END