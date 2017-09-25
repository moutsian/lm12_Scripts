# filter CD results files according to Katie's QC files.
#chr=3
trait=tolower("ibd")
TRAIT=toupper(trait)
#final_meta_analysis_C.UC.1e_05.txt
if(file.exists(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".1e_05.txt",sep=""))){
results=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".1e_05.txt",sep=""),head=F, sep=" ")
resultso=results[order(as.numeric(as.character(results[,1])),decreasing=F),]
to_include=matrix(ncol=1,nrow=dim(results)[1],1)
info_fail=0;
maf_fail=0;
for(chr in 1:22){
resultso_subset=resultso[as.numeric(as.character(resultso[,1]))==chr,]
GWAS3subset=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/",trait,"/",chr,".assoc",sep=""),head=T)
IIBDGCsubset=read.table(paste("/lustre/scratch113/projects/crohns/iibdgc_meta/data/IIBDGC/",trait,"/",chr,".assoc",sep=""),head=T)
IBDseqsubset=read.table(paste("/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/",trait,"/",chr,".assoc",sep=""),head=T)

GWAS3subset=GWAS3subset[c(4,9,31)]
IBDseqsubset=IBDseqsubset[c(4,9,31)]
IIBDGCsubset=IIBDGCsubset[c(3,7,10)]

IIBDGCsubset_filt=IIBDGCsubset[as.numeric(as.character(IIBDGCsubset[,3]))>0.005 & as.numeric(as.character(IIBDGCsubset[,3]))<0.995 ,]
IBDseqsubset_filt=IBDseqsubset[as.numeric(as.character(IBDseqsubset[,3]))>0.005 & as.numeric(as.character(IBDseqsubset[,3]))<0.995 ,]
GWAS3subset_filt=GWAS3subset[as.numeric(as.character(GWAS3subset[,3]))>0.005 & as.numeric(as.character(GWAS3subset[,3]))<0.995 ,]

********* WORK IN PROGRESS **** THE FOR LOOP BELOW WILL BE REMOVED

for(i in 1: dim(resultso_subset)[1]){

if(is.null(shelloutIIBDGC) & is.null(shelloutGWAS3) & is.null(shellcmdIBDseq)){
print("something is wrong - this variant is not in any of the datasets")
}
if(is.null(shelloutGWAS3)){
	if(is.null(shelloutIIBDGC)){
	shelloutGWAS3=shelloutIBDseq
	}else{
	shelloutGWAS3=shelloutIIBDGC
	}
}
if(is.null(shelloutIIBDGC)){
	if(is.null(shelloutGWAS3)){
	shelloutIIBDGC=shelloutIBDseq
	}else{
	shelloutIIBDGC=shelloutGWAS3
	}
}
if(is.null(shelloutIBDseq)){
	if(is.null(shelloutGWAS3)){
	shelloutIBDseq=shelloutIIBDGC
	}else{
	shelloutIBDseq=shelloutGWAS3
	}
}

if(as.numeric(shelloutIIBDGC[1])<0.4 | as.numeric(shelloutGWAS3[1])<0.4 | as.numeric(shelloutIBDseq[1])<0.4 ){
to_include[i]=0;
info_fail=info_fail+1;
}

if(as.numeric(shelloutIIBDGC[2])<0.005 | as.numeric(shelloutGWAS3[2])<0.005 | as.numeric(shelloutIBDseq[2])<0.005 | as.numeric(shelloutIIBDGC[2])>0.995 | as.numeric(shelloutGWAS3[2])>0.995 | as.numeric(shelloutIBDseq[2])>0.995 ){
to_include[i]=0;
maf_fail=maf_fail+1;
}
}
}
write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".1e_05.filtered.txt",sep=""),quote=F,row.names=F,col.names=T)
}

#END

