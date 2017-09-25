
novel=read.table("/nfs/users/nfs_l/lm12/IBD_conditional/novel_hits_to_get_extended_ld_info_for.txt",sep=" ")
previous=read.table("/nfs/users/nfs_l/lm12/IBD_conditional/step2.ALL.round6.with_info_on_previous_hits_from_all_PAPERS.final.txt",head=T,sep="\t")
previousreported=previous[previous[,20]!="",]

ldcheck=NULL
for(chr in 1:22){
print(paste("chr: ",chr,sep=""))
positions=novel[which(novel[,2]==chr),4]
if(length(positions)>0){
for(i in 1:length(positions)){
oldidx=which(as.numeric(previousreported[,2])==chr & abs(as.numeric(as.character(previousreported[,4]))-as.numeric(positions[i]))<5000000)
	if(length(oldidx)>0){
	unique_old=unique(as.character(previousreported[oldidx,4]))
	minpos=min(as.numeric(positions[i]),as.numeric(unique_old))
	maxpos=max(as.numeric(positions[i]),as.numeric(unique_old))	
	shellcommand1=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/check_ld_between_two_loci.sh",chr,minpos,maxpos,sep=" ")
	shellout=system(shellcommand1,intern=F)
	for(j in 1:length(unique_old)){
	shellcommand2=paste(" grep ",positions[i]," /lustre/scratch115/teams/anderson/ibd_conditional/plink_ld.out.",chr,".",minpos,".",maxpos,".ld |grep ",unique_old[j],sep="")
	shellout=system(shellcommand2,intern=T)
		if(length(shellout)==0){
		ldcheck=rbind(ldcheck,c(positions[i],unique_old[j],"<0.1"))
		}else{
		ldcheck=rbind(ldcheck,c(positions[i],unique_old[j],unlist(strsplit(shellout,"\\s+"))[8]))
		}
	}
	}else{
		ldcheck=rbind(ldcheck,c(positions[i],"Nothing_in_a_5Mb_window","0"))
	}
}

}else{
print(paste("no putatively novel hits in chromosome ",chr,sep=""))
}
}
#END