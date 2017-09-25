## adds info scores to the vepres loci from Dataset C
 

headercmdGWAS3="head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/ibd/21.assoc"
headerGWAS3=unlist(strsplit(system(headercmdGWAS3,intern=T)," "))
headercmdIIBDGC="head -n1 /lustre/scratch113/projects/crohns/iibdgc_meta/data/IIBDGC/ibd/21.assoc"
headerIIBDGC=unlist(strsplit(system(headercmdIIBDGC,intern=T),"\t"))
headercmdIBDseq="head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/ibd/21.assoc"
headerIBDseq=unlist(strsplit(system(headercmdIBDseq,intern=T)," "))
vepres= read.table("~lm12/IBD_conditional/Vep_results_datasetB.round6.ALL.high_and_moderate_impact_only.txt",head=T,sep="\t")
info=matrix(nrow=dim(vepres)[1],ncol=3,NA)
for(i in 1: dim(vepres)[1]){
print(paste("i: ",i,sep=""))

if(tolower(as.character(vepres[i,1]))=="all"){
vepres[i,1]="ibd"}

shellcmdGWAS3=paste("zgrep ",unlist(strsplit(as.character(vepres[i,5]),"_"))[1]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/",tolower(as.character(vepres[i,2])),"/",as.character(vepres[i,3]),".assoc",sep="")
shellcmdIIBDGC=paste("zgrep ",unlist(strsplit(as.character(vepres[i,5]),"_"))[1]," /lustre/scratch113/projects/crohns/iibdgc_meta/data/IIBDGC/",tolower(as.character(vepres[i,2])),"/",as.character(vepres[i,3]),".assoc",sep="")
shellcmdIBDseq=paste("zgrep ",unlist(strsplit(as.character(vepres[i,5]),"_"))[1]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/",tolower(as.character(vepres[i,2])),"/",as.character(vepres[i,3]),".assoc",sep="")

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

infoGWAS3=shelloutGWAS3[9]
infoIBDseq=shelloutIBDseq[9]
infoIIBDGC=shelloutIIBDGC[7]
if(is.null(infoIIBDGC)){
infoIIBDGC="-"
}
if(is.null(infoIBDseq)){
infoIBDseq="-"
}
if(is.null(infoGWAS3)){
infoGWAS3="-"
}
info[i,1:3]=as.matrix(c(infoIIBDGC,infoGWAS3,infoIBDseq))
}
colnames(info)=c("infoIIBDGC","infoGWAS3","infoIBDseq")
vepresplus=cbind(vepres,info)
write.table(vepresplus,"~lm12/IBD_conditional/Vep_results_datasetB.round6.ALL.high_and_moderate_impact_only.txt.with_info.txt",sep="\t",quote=F,col.names=T,row.names=F)
#END