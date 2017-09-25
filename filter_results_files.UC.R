# filter CD results files according to Katie's QC files.
#chr=3
trait=tolower("uc")
TRAIT=toupper(trait)
for(chr in 1:22){
if(file.exists(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""))){
results=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/",TRAIT,"/results.",TRAIT,".",chr,".txt",sep=""),head=T, sep=" ")
if(dim(results)[1]>1){
to_include=matrix(ncol=1,nrow=dim(results)[1],1)
info_fail=0;
maf_fail=0;
for(i in 1:length(to_include)){
print(paste("i: ",i,sep=""))
shellcmdGWAS3=paste("zgrep ",unlist(strsplit(as.character(results[i,2]),"_"))[1]," /lustre/scratch113/projects/crohns/iibdgc_meta/GWAS3/",trait,"/",chr,".maf",sep="")
shellcmdIIBDGC=paste("zgrep ",unlist(strsplit(as.character(results[i,2]),"_"))[1]," /lustre/scratch113/projects/crohns/iibdgc_meta/IIBDGC/",trait,"/",chr,".maf",sep="")
shellcmdIBDseq=paste("zgrep ",unlist(strsplit(as.character(results[i,2]),"_"))[1]," /lustre/scratch113/projects/crohns/iibdgc_meta/IBD_seq/",trait,"/",chr,".maf",sep="")

#the headers for the files are:
#GWAS3 (all three files):
#alternate_ids info cases_maf controls_maf
#IBDseq (all three files):
#rsid info cases_maf controls_maf
#IIBDGC (all three files):
#IBD - CHR:BP_A1_A2    FRQ_A_12882    FRQ_U_21770    INFO
#CD - CHR:BP_A1_A2    FRQ_A_5956    FRQ_U_14927    INFO
#UC - CHR:BP_A1_A2    FRQ_A_6968    FRQ_U_20464    INFO

shelloutGWAS3=unlist(strsplit(system(shellcmdGWAS3,intern=T)," "))
shelloutIIBDGC=unlist(strsplit(system(shellcmdIIBDGC,intern=T),"\t"))
shelloutIBDseq=unlist(strsplit(system(shellcmdIBDseq,intern=T)," "))

if(as.numeric(shelloutIIBDGC[4])<0.4 | as.numeric(shelloutGWAS3[2])<0.4 | as.numeric(shelloutIBDseq[2])<0.4 ){
to_include[i]=0;
info_fail=info_fail+1;
}

if(as.numeric(shelloutIIBDGC[3])<0.005 | as.numeric(shelloutGWAS3[4])<0.005 | as.numeric(shelloutIBDseq[4])<0.005 | as.numeric(shelloutIIBDGC[3])>0.995 | as.numeric(shelloutGWAS3[4])>0.995 | as.numeric(shelloutIBDseq[4])>0.995 ){
to_include[i]=0;
maf_fail=maf_fail+1;
}
}
write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/",TRAIT,"/results.",TRAIT,".",chr,".txt.filtered",sep=""),quote=F,row.names=F,col.names=T)
}}
}
#END

