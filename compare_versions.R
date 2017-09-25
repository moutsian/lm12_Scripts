#this is for dataset B

new=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.CD.1e_05.txt")
oldfilt=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results_beforePCcorrection/UK_only_analysis_B.CD.1e_05.lowfreq.txt",head=T)
#chr=20
for(chr in 1:22){
oldsub=oldfilt[which(as.numeric(as.character(oldfilt[,1]))==chr & as.numeric(oldfilt[,12])<80 & as.numeric(oldfilt[,12])!=-1),]
newsub=new[which(as.numeric(as.character(new[,1]))==chr & as.numeric(new[,12])<80 & as.numeric(new[,12])!=-1),]
novel=which(!(as.character(newsub[,3])%in%as.character(oldsub[,3])))
disappeared=which(!(as.character(oldsub[,3])%in%as.character(newsub[,3])))
tocheck_new=which(as.numeric(as.character(newsub[novel,6]))<5e-08)
tocheck_old=which(as.numeric(as.character(oldsub[disappeared,6]))<5e-08)

results=NULL
if(length(tocheck_new)>0){
for(i in 1:length(tocheck_new)){
shellcmdOLD=paste("zgrep ",newsub[novel[tocheck_new[i]],3]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL/CD/",newsub[novel[tocheck_new[i]],1],"-meta.txt | head -n1",sep="")
shelloutOLD=unlist(strsplit(system(shellcmdOLD,intern=T)," "))
results=rbind(results,c(newsub[novel[tocheck_new[i]],6],shelloutOLD[2],shelloutOLD[3],shelloutOLD[6],shelloutOLD[12]))
}
}

#now the other way around
results_disappeared=NULL
if(length(tocheck_old)>0){
for(i in 1:length(tocheck_old)){
shellcmdNEW=paste("zgrep ",oldsub[disappeared[tocheck_old[i]],3]," /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/UKIBDGC_meta_B/cd/",oldsub[disappeared[tocheck_old[i]],1],"-meta.txt | head -n1",sep="")
shelloutNEW=unlist(strsplit(system(shellcmdNEW,intern=T)," "))
results_disappeared=rbind(results_disappeared,c(oldsub[disappeared[tocheck_old[i]],6],shelloutNEW[2],shelloutNEW[3],shelloutNEW[6],shelloutNEW[12]))
}
}

write.table(results_disappeared,paste("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/check_hits_between_versions_before_and_after_PC_correction/CD.chr",chr,".disappeared",sep=""),quote=F,col.names=F)
write.table(results,paste("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/check_hits_between_versions_before_and_after_PC_correction/CD.chr",chr,".novel",sep=""),quote=F,col.names=F)
}
#END