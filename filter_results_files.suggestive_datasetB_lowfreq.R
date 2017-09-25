# filter CD results files according to Katie's QC files.
#chr=3
trait=tolower("ibd")
TRAIT=toupper(trait)
#UK_only_analysis_B.UC.1e_05.txt
if(file.exists(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.txt",sep=""))){
results=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.txt",sep=""),head=F, sep=" ")
to_include=matrix(ncol=1,nrow=dim(results)[1],1)
qc_fail=0;
if(TRAIT=="IBD"){
	for(i in 1: length(to_include)){
	print(paste("i: ",i,sep=""))
	shellcmd=paste("zgrep ",as.character(results[i,3])," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/QC_01/",as.character(results[i,1]),".txt.gz|grep ",as.character(results[i,4])," | grep ",as.character(results[i,5]),sep="")
	#the header for the files is:
	# VARIANTID GWAS1_PASS GWAS2_PASS GWAS3_CD_PASS GWAS3_UC_PASS GWAS3_IBD_PASS IBDSEQ_CD_PASS IBDSEQ_UC_PASS IBDSEQ_IBD_PASS
	shellout=unlist(strsplit(system(shellcmd,intern=T)," "))

	if(sum(as.numeric(shellout[2:length(shellout)]),na.rm=T)<8){
		to_include[i]=0;
		qc_fail=qc_fail+1;
		}
	}
	write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.filtered.lowfreq.txt",sep=""),quote=F,row.names=F,col.names=T)
	
}else if(TRAIT=="UC"){

for(i in 1: length(to_include)){
	print(paste("i: ",i,sep=""))
	shellcmd=paste("zgrep ",as.character(results[i,3])," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/QC_01/",as.character(results[i,1]),".txt.gz|grep ",as.character(results[i,4])," | grep ",as.character(results[i,5]),sep="")
	shellout=unlist(strsplit(system(shellcmd,intern=T)," "))

	if(sum(as.numeric(shellout[c(3,5,8)]),na.rm=T)<3){
		to_include[i]=0;
		qc_fail=qc_fail+1;
	}
	}
	write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.filtered.lowfreq.txt",sep=""),quote=F,row.names=F,col.names=T)

}else if(TRAIT=="CD"){

for(i in 1: length(to_include)){
	print(paste("i: ",i,sep=""))
	shellcmd=paste("zgrep ",as.character(results[i,3])," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/QC_01/",as.character(results[i,1]),".txt.gz|grep ",as.character(results[i,4])," | grep ",as.character(results[i,5]),sep="")
	shellout=unlist(strsplit(system(shellcmd,intern=T)," "))

	if(sum(as.numeric(shellout[c(2,4,7)]),na.rm=T)<3){
		to_include[i]=0;
		qc_fail=qc_fail+1;
	}
	}
	write.table(results[to_include==1,],paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",TRAIT,".1e_05.filtered.lowfreq.txt",sep=""),quote=F,row.names=F,col.names=T)

	}else{
print("Something wrong!");
}
}

#END

