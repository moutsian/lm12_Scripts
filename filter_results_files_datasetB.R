# filter dataset B results files according to Katie's QC files. - updated to also filter out SNPs with deviations from HWE
#chr=3
trait=tolower("ibd")
TRAIT=toupper(trait)
GWAS1outofHWE=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.filtered.lowfreq.HWE.REMOVED.txt",sep=""),head=T)
GWAS2outofHWE=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.filtered.lowfreq.HWE.REMOVED.txt",sep=""),head=T)
GWAS3outofHWE=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.filtered.lowfreq.HWE.REMOVED.txt",sep=""),head=T)
GWASoutofHWE=NULL
if(TRAIT=="UC"){
GWASoutofHWE=rbind(GWAS2outofHWE,GWAS3outofHWE)
}else if(TRAIT=="IBD"){
GWASoutofHWE=rbind(GWAS1outofHWE,GWAS2outofHWE,GWAS3outofHWE)
}else if(TRAIT=="CD"){
GWASoutofHWE=rbind(GWAS1outofHWE,GWAS3outofHWE)
}else{
print("something went wrong when selecting trait")
}

if(TRAIT=="IBD"){
for(chr in 1:22){
if(file.exists(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""))){
results=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T, sep="\t")
if(dim(results)[1]>1){
to_include=matrix(ncol=1,nrow=dim(results)[1],1)
qc_fail=0;
hwe_fail=0;
for(i in 1:length(to_include)){
print(paste("i: ",i,sep=""))
shellcmd=paste("zgrep ",paste(chr,":",as.character(results[i,3]),"_",sep="")," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/QC/",chr,".txt.gz",sep="")

#the header for the QC file is:
# VARIANTID GWAS1_PASS GWAS2_PASS GWAS3_CD_PASS GWAS3_UC_PASS GWAS3_IBD_PASS IBDSEQ_CD_PASS IBDSEQ_UC_PASS IBDSEQ_IBD_PASS

shellout=unlist(strsplit(system(shellcmd,intern=T)," "))
if(sum(as.numeric(shellout[2:length(shellout)]),na.rm=T)<8){
to_include[i]=0;
qc_fail=qc_fail+1;
}
#last addition: remove if out of HWE
idx=which(as.character(results[i,3])==as.character(GWASoutofHWE[,4]) & as.character(results[i,1])==as.character(GWASoutofHWE[,2]) )
if(length(idx)>0){
to_include[i]=0;
hwe_fail=hwe_fail+1;
print(paste("hwe_fail:",hwe_fail,sep=""))}
}
write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt.filtered",sep=""),quote=F,row.names=F,col.names=T)
}}
}
}else if(TRAIT=="CD"){
for(chr in 1:22){
if(file.exists(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""))){
results=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T, sep="\t")
if(dim(results)[1]>1){
to_include=matrix(ncol=1,nrow=dim(results)[1],1)
qc_fail=0;
hwe_fail=0;
for(i in 1:length(to_include)){
print(paste("i: ",i,sep=""))
shellcmd=paste("zgrep ",paste(chr,":",as.character(results[i,3]),"_",sep="")," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/QC/",chr,".txt.gz",sep="")

#the header for the QC file is:
# VARIANTID GWAS1_PASS GWAS2_PASS GWAS3_CD_PASS GWAS3_UC_PASS GWAS3_IBD_PASS IBDSEQ_CD_PASS IBDSEQ_UC_PASS IBDSEQ_IBD_PASS

shellout=unlist(strsplit(system(shellcmd,intern=T)," "))
if(sum(as.numeric(shellout[c(2,4,7)]),na.rm=T)<3){
to_include[i]=0;
qc_fail=qc_fail+1;
}
#last addition: remove if out of HWE
idx=which(as.character(results[i,3])==as.character(GWASoutofHWE[,4]) & as.character(results[i,1])==as.character(GWASoutofHWE[,2]) )
if(length(idx)>0){
to_include[i]=0;
hwe_fail=hwe_fail+1;
print(paste("hwe_fail:",hwe_fail,sep=""))}
}
write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt.filtered",sep=""),quote=F,row.names=F,col.names=T)
}}
}
}else if(TRAIT=="UC"){
for(chr in 1:22){
if(file.exists(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""))){
results=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T, sep="\t")
if(dim(results)[1]>1){
to_include=matrix(ncol=1,nrow=dim(results)[1],1)
qc_fail=0;
hwe_fail=0;
for(i in 1:length(to_include)){
print(paste("i: ",i,sep=""))
shellcmd=paste("zgrep ",paste(chr,":",as.character(results[i,3]),"_",sep="")," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/QC/",chr,".txt.gz",sep="")

#the header for the QC file is:
# VARIANTID GWAS1_PASS GWAS2_PASS GWAS3_CD_PASS GWAS3_UC_PASS GWAS3_IBD_PASS IBDSEQ_CD_PASS IBDSEQ_UC_PASS IBDSEQ_IBD_PASS

shellout=unlist(strsplit(system(shellcmd,intern=T)," "))
if(sum(as.numeric(shellout[c(3,5,8)]),na.rm=T)<3){
to_include[i]=0;
qc_fail=qc_fail+1;
}
#last addition: remove if out of HWE
idx=which(as.character(results[i,3])==as.character(GWASoutofHWE[,4]) & as.character(results[i,1])==as.character(GWASoutofHWE[,2]) )
if(length(idx)>0){
to_include[i]=0;
hwe_fail=hwe_fail+1;
print(paste("hwe_fail:",hwe_fail,sep=""))}
}
write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/",TRAIT,"/results.",TRAIT,".",chr,".txt.filtered",sep=""),quote=F,row.names=F,col.names=T)
}}
}
}else{
print("kati paixtike...")
}
#END

