# filter CD results files according to Katie's QC files.
#chr=3
trait=tolower("cd")
TRAIT=toupper(trait)
#final_meta_analysis_C.UC.1e_05.txt
if(file.exists(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".1e_05.txt",sep=""))){
results=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".1e_05.txt",sep=""),head=F, sep="\t")
to_include=matrix(ncol=1,nrow=dim(results)[1],1)
info_fail=0;
maf_fail=0;
for(i in 1: length(to_include)){
print(paste("i: ",i,sep=""))

#if(i==1){
#headercmdGWAS3=paste("head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/",trait,"/",as.character(results[i,1]),".assoc",sep="")
#headerGWAS3=unlist(strsplit(system(headercmdGWAS3,intern=T)," "))
#headercmdIIBDGC=paste("head -n1 /lustre/scratch113/projects/crohns/iibdgc_meta/data/IIBDGC/",trait,"/",as.character(results[i,1]),".assoc",sep="")
#headerIIBDGC=unlist(strsplit(system(headercmdIIBDGC,intern=T),"\t"))
#headercmdIBDseq=paste("head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/",trait,"/",as.character(results[i,1]),".assoc",sep="")
#headerIBDseq=unlist(strsplit(system(headercmdIBDseq,intern=T)," "))
#}

POS=results[i,2]
CHROM=as.character(results[i,1])
shellcmdGWAS3=paste("grep ",POS," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/",trait,"/",CHROM,".assoc",sep="")
shellcmdIIBDGC=paste("grep ",POS," /lustre/scratch113/projects/crohns/iibdgc_meta/data/IIBDGC/",trait,"/",CHROM,".assoc",sep="")
shellcmdIBDseq=paste("grep ",POS," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/",trait,"/",CHROM,".assoc",sep="")

#the headers for the files are:
#GWAS3 (all three files):
#alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total cases_AA cases_AB cases_BB cases_NULL cases_total controls_AA controls_AB controls_BB controls_NULL controls_total all_maf cases_maf controls_maf missing_data_proportion het_OR het_OR_lower het_OR_upper hom_OR hom_OR_lower hom_OR_upper all_OR all_OR_lower all_OR_upper frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1 comment
#IBDseq (all three files):
#alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total cases_AA cases_AB cases_BB cases_NULL cases_total controls_AA controls_AB controls_BB controls_NULL controls_total all_maf cases_maf controls_maf missing_data_proportion het_OR het_OR_lower het_OR_upper hom_OR hom_OR_lower hom_OR_upper all_OR all_OR_lower all_OR_upper frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1 comment
#IIBDGC (all three files):
#/lustre/scratch113/projects/crohns/iibdgc_meta/data/IIBDGC/cd/
#rsid    chr     pos     allele_A        allele_B        P_value info    all_maf cases_maf       controls_maf    beta    se


shelloutGWAS3=unlist(strsplit(system(shellcmdGWAS3,intern=T)," "))
shelloutIIBDGC=unlist(strsplit(system(shellcmdIIBDGC,intern=T),"\t"))
shelloutIBDseq=unlist(strsplit(system(shellcmdIBDseq,intern=T)," "))

shelloutGWAS3=shelloutGWAS3[c(9,31)]
shelloutIBDseq=shelloutIBDseq[c(9,31)]
shelloutIIBDGC=shelloutIIBDGC[c(7,10)]


#if(is.null(shelloutIIBDGC) & is.null(shelloutGWAS3) & is.null(shellcmdIBDseq)){
#print("something is wrong - this variant is not in any of the datasets")
#}
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
write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.",TRAIT,".1e_05.filtered.txt",sep=""),quote=F,row.names=F,col.names=T)
}

#END

