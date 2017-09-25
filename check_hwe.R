## updated so it uses control files.

LOWFREQ=0;
trait="uc"
utrait=toupper(trait);
ltrait=tolower(trait);
#a HWE threshold of 1e-07 has been used. Note that what is reported in the HWE column is -log10P for HWE.
GWAS1out=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS1.outofHWE.controls.txt",head=T)
GWAS2out=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS2.outofHWE.controls.txt",head=T)
GWAS3out=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS3.outofHWE.controls.txt",head=T)
if(LOWFREQ==1){
ourdata=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.txt",sep=""),head=T)
}else{
ourdata=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.txt",sep=""),head=T)
}
toremove=NULL
REMOVED=NULL
ourdata_idx=which(as.character(ourdata[,3])%in%as.character(GWAS1out[,4]))
ourdata_idx=c(ourdata_idx,which(as.character(ourdata[,3])%in%as.character(GWAS2out[,4])))
ourdata_idx=c(ourdata_idx,GWAS3idx=which(as.character(ourdata[,3])%in%as.character(GWAS3out[,4])))
ourdata_idx=unique(ourdata_idx)
for(i in 1:length(ourdata_idx)){
gwas1_idx=which( as.character(GWAS1out[,1])==as.character(ourdata[ourdata_idx[i],1]) & as.character(GWAS1out[,4])==as.character(ourdata[ourdata_idx[i],3]) )
gwas2_idx=which( as.character(GWAS2out[,1])==as.character(ourdata[ourdata_idx[i],1]) & as.character(GWAS2out[,4])==as.character(ourdata[ourdata_idx[i],3]) )
gwas3_idx=which( as.character(GWAS3out[,1])==as.character(ourdata[ourdata_idx[i],1]) & as.character(GWAS3out[,4])==as.character(ourdata[ourdata_idx[i],3]) )
if(length(gwas1_idx)>0){
toremove=c(toremove,ourdata_idx[i])
REMOVED=rbind(REMOVED,c("GWAS1",as.matrix(ourdata[ourdata_idx[i],]),as.matrix(GWAS1out[gwas1_idx,])))
}else if(length(gwas2_idx)>0){
toremove=c(toremove,ourdata_idx[i])
REMOVED=rbind(REMOVED,c("GWAS2",as.matrix(ourdata[ourdata_idx[i],]),as.matrix(GWAS2out[gwas2_idx,])))
}else if(length(gwas3_idx)>0){
toremove=c(toremove,ourdata_idx[i])
REMOVED=rbind(REMOVED,c("GWAS3",as.matrix(ourdata[ourdata_idx[i],]),as.matrix(GWAS3out[gwas3_idx,])))
}
}

if(utrait=="IBD"){
tmpnames=c("chr","rsid","pos","allele_A","allele_B","P_value","coded_af","beta","se","Q","P_heterogeneity","I2","P_cohort_1","P_cohort_2","P_cohort_3","P_cohort_4")
}else{
tmpnames=c("chr","rsid","pos","allele_A","allele_B","P_value","coded_af","beta","se","Q","P_heterogeneity","I2","P_cohort_1","P_cohort_2","P_cohort_3")
}
colnames(REMOVED)=c("cohort_out_of_HWE",tmpnames,colnames(GWAS1out))
#filtered output
ourdata_filt=ourdata[-toremove,]
if(LOWFREQ==1){
write.table(ourdata_filt,paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.HWE.PASSED.txt",sep=""),quote=F,row.names=F,col.names=F)
write.table(REMOVED,paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.HWE.REMOVED.txt",sep=""),quote=F,row.names=F,col.names=T)
}else{
write.table(ourdata_filt,paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.HWE.PASSED.txt",sep=""),quote=F,row.names=F,col.names=F)
write.table(REMOVED,paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.HWE.REMOVED.txt",sep=""),quote=F,row.names=F,col.names=T)
}


#END
