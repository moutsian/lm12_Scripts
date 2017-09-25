michel=read.table("~lm12/IBD_conditional/Michel_variants_2016_03_01.txt",sep="\t",head=T)
DATAOUT=matrix(nrow=dim(michel)[1],ncol=15,NA)
DATAOUT[,1:3]=as.matrix(michel[,c(1,2,22)])
for(i in 1:dim(michel)[1]){

	commandA=paste("grep ",michel[i,2]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/IBDseq/CD/",michel[i,1],".assoc",sep="")
	commandB=paste("grep ",michel[i,2]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/CD/",michel[i,1],"-meta.txt",sep="")
	commandC=paste("grep ",michel[i,2]," /lustre/scratch113/projects/crohns/iibdgc_meta/results/cd-meta-filtered.out",sep="")
	headerA=paste("grep chromosome /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/IBDseq/CD/",michel[i,1],".assoc",sep="")
	headerB=paste("grep chr /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/CD/",michel[i,1],"-meta.txt",sep="")
	headerC=paste("grep chr /lustre/scratch113/projects/crohns/iibdgc_meta/results/cd-meta-filtered.out",sep="")


	headerA=unlist(strsplit(system(headerA,intern=T)," "))
	shelloutA=unlist(strsplit(system(commandA,intern=T)," "))
	if(!is.null(shelloutA)){
	DATAOUT[i,4]=shelloutA[5] #alA
	DATAOUT[i,5]=shelloutA[6] #alB
	DATAOUT[i,6]=shelloutA[44] #beta
	DATAOUT[i,7]=shelloutA[45] #se
	DATAOUT[i,8]=shelloutA[42] #pval
	DATAOUT[i,9]=shelloutA[9] #info
	}
	headerB=unlist(strsplit(system(headerB,intern=T)," "))
	shelloutB=unlist(strsplit(system(commandB,intern=T)," "))
	if(!is.null(shelloutB)){
	DATAOUT[i,10]=shelloutB[4]
	DATAOUT[i,11]=shelloutB[5]	
	DATAOUT[i,12]=shelloutB[8]	
	DATAOUT[i,13]=shelloutB[9]	
	DATAOUT[i,14]=shelloutB[6]	
	DATAOUT[i,15]=shelloutB[12]	
	}
	###For now we are only reporting datasets A and B
	#headerC=unlist(strsplit(system(headerC,intern=T)," "))
	#shelloutC=unlist(strsplit(system(commandC,intern=T)," "))
	#if(!is.null(shelloutC)){
	#DATAOUT[i,16]=shelloutC[4]
	#DATAOUT[i,17]=shelloutC[5]	
	#DATAOUT[i,18]=shelloutC[7]	
	#DATAOUT[i,19]=shelloutC[8]	
	#DATAOUT[i,20]=shelloutC[6]	
	#DATAOUT[i,21]=shelloutC[11]	
	#}
}

colnames(DATAOUT)=c("Chr","Pos","Michel_pvalue","DatasetA_alleleA","DatasetA_alleleB","DatasetA_beta","DatasetA_se","DatasetA_pval","DatasetA_info","DatasetB_alleleA","DatasetB_alleleB","DatasetB_beta","DatasetB_se","DatasetB_pvalue","DatasetB_I2")
write.table(DATAOUT,"~lm12/IBD_conditional/variant_lookup_for_michel_datasetsABC.updated.txt",col.names=T,row.names=F,quote=F,sep="\t")
#END