## adds info scores to the novel loci from Dataset B
 

headercmdGWAS1="head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS1/cd/21.assoc"
headerGWAS1=unlist(strsplit(system(headercmdGWAS1,intern=T)," "))
headercmdGWAS2="head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS2/uc/21.assoc"
headerGWAS2=unlist(strsplit(system(headercmdGWAS2,intern=T)," "))
headercmdGWAS3="head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/ibd/21.assoc"
headerGWAS3=unlist(strsplit(system(headercmdGWAS3,intern=T)," "))
headercmdIBDseq="head -n1 /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/ibd/21.assoc"
headerIBDseq=unlist(strsplit(system(headercmdIBDseq,intern=T)," "))

novel=read.table("~lm12/IBD_conditional/step2.ALL.dataset_B.NOVEL.round6.with_info_on_previous_hits_from_all_PAPERS.auto.txt",head=T,sep="\t")
#novel=read.table("~lm12/IBD_conditional/step2.ALL.NOVEL.round6.with_info_on_previous_hits_from_all_PAPERS.txt",head=T,sep="\t")
info=matrix(nrow=dim(novel)[1],ncol=4,NA)
for(i in 1: dim(novel)[1]){
print(paste("i: ",i,sep=""))

if(tolower(as.character(novel[i,1]))=="all"){
novel[i,1]="ibd"}

shellcmdGWAS1=paste("zgrep ",unlist(strsplit(as.character(novel[i,4]),"_"))[1]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS1/cd/",as.character(novel[i,2]),".assoc",sep="")
shellcmdGWAS2=paste("zgrep ",unlist(strsplit(as.character(novel[i,4]),"_"))[1]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS2/uc/",as.character(novel[i,2]),".assoc",sep="")
shellcmdGWAS3=paste("zgrep ",unlist(strsplit(as.character(novel[i,4]),"_"))[1]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/",tolower(as.character(novel[i,1])),"/",as.character(novel[i,2]),".assoc",sep="")
shellcmdIBDseq=paste("zgrep ",unlist(strsplit(as.character(novel[i,4]),"_"))[1]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/IBDSeq_no_IIBDGC/",tolower(as.character(novel[i,1])),"/",as.character(novel[i,2]),".assoc",sep="")

shelloutGWAS1=unlist(strsplit(system(shellcmdGWAS1,intern=T)," "))
shelloutGWAS2=unlist(strsplit(system(shellcmdGWAS2,intern=T)," "))
shelloutGWAS3=unlist(strsplit(system(shellcmdGWAS3,intern=T)," "))
shelloutIBDseq=unlist(strsplit(system(shellcmdIBDseq,intern=T)," "))

infoGWAS1=shelloutGWAS1[9]
infoGWAS2=shelloutGWAS2[9]
infoGWAS3=shelloutGWAS3[9]
infoIBDseq=shelloutIBDseq[9]
if(is.null(infoGWAS1)){
infoGWAS1="-"
}
if(is.null(infoGWAS2)){
infoGWAS2="-"
}
if(is.null(infoIBDseq)){
infoIBDseq="-"
}
if(is.null(infoGWAS3)){
infoGWAS3="-"
}
info[i,1:4]=as.matrix(c(infoGWAS1,infoGWAS2,infoGWAS3,infoIBDseq))
}
colnames(info)=c("infoGWAS1","infoGWAS2","infoGWAS3","infoIBDseq")
novelplus=cbind(novel,info)
write.table(novelplus,"~lm12/IBD_conditional/step2.ALLTRAITS.dataset_B.NOVEL.round6.with_info_on_previous_hits_from_all_PAPERS.auto.with_info.txt",sep="\t",quote=F,col.names=T,row.names=F)
#END