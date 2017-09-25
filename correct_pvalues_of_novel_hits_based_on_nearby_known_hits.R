
#All hits:
IBD=read.table("~lm12/IBD_conditional/step2.IBD.all.round2.with_info_on_previous_hits_from_5PAPERS.txt",head=T,sep="\t")
CD=read.table("~lm12/IBD_conditional/step2.CD.all.round2.with_info_on_previous_hits_from_5PAPERS.txt",head=T,sep="\t")
UC=read.table("~lm12/IBD_conditional/step2.UC.all.round2.with_info_on_previous_hits_from_5PAPERS.txt",head=T,sep="\t")

known_hits=rbind(IBD,CD,UC)
known_hits=known_hits[known_hits[,16]!="",]

#Novel hits:

Trait="UC"
NOVEL=read.table(paste("~lm12/IBD_conditional/step2.",Trait,".NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.txt",sep=""),head=T,sep="\t")

CLOSER_LOCI=matrix(ncol=4,nrow=dim(NOVEL)[1],NA)
for(i in 1:dim(NOVEL)[1]){
idx=which(as.numeric(NOVEL[i,2])==as.numeric(known_hits[,1]))
known_hits_samechr=known_hits[idx,]
distance_right=as.numeric(known_hits_samechr[,3])-as.numeric(NOVEL[i,4])
distance_left=as.numeric(NOVEL[i,4])-as.numeric(known_hits_samechr[,3])
distance_left[distance_left<0]=NA
distance_right[distance_right<0]=NA

pos_left=NA
if(sum(is.na(distance_left))!=length(distance_left)){
CLOSER_LOCI[i,3]=as.character(known_hits_samechr[which.min(distance_left),3]) #variant_to_include_on the left
}
pos_right=NA
if(sum(is.na(distance_right))!=length(distance_right)){
CLOSER_LOCI[i,4]=as.character(known_hits_samechr[which.min(distance_right),3]) #variant_to_include_on the right
}
CLOSER_LOCI[i,1]=as.numeric(NOVEL[i,2]) #chr
CLOSER_LOCI[i,2]=as.numeric(NOVEL[i,4]) #novel locus main variant
#print(paste("i: ",i,", for ",NOVEL[i,2],":",NOVEL[i,4],", left= ",as.character(known_hits_samechr[pos_left,3]),", right: ",as.character(known_hits_samechr[pos_right,3]),sep=""))
}

colnames(CLOSER_LOCI)=c("CHR","POS_NOVEL_LOCUS_VARIANT","POS_KNOWN_VAR_LEFT","POS_KNOWN_VAR_RIGHT")
#CLOSER_LOCI=CLOSER_LOCI[order(as.numeric(CLOSER_LOCI[,1]),as.numeric(CLOSER_LOCI[,2])),]

R2LEFT=matrix(ncol=1,nrow=dim(CLOSER_LOCI)[1],NA)
R2RIGHT=matrix(ncol=1,nrow=dim(CLOSER_LOCI)[1],NA)
for(i in 1:dim(CLOSER_LOCI)[1]){

	#now add r2 values based on 1KG CEU & GBR 
	if(!is.na(CLOSER_LOCI[i,3])){
	shellcommandLEFT=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/calc_pairwise_ld.sh ",as.character(CLOSER_LOCI[i,1]),as.character(CLOSER_LOCI[i,2]),as.character(CLOSER_LOCI[i,3]),sep=" ")
	system(shellcommandLEFT)
	ldfile=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/plink_ld.out.",as.character(CLOSER_LOCI[i,1]),".",as.character(min(as.numeric(CLOSER_LOCI[i,2]),as.numeric(CLOSER_LOCI[i,3]))),".",as.character(max(as.numeric(CLOSER_LOCI[i,2]),as.numeric(CLOSER_LOCI[i,3]))),".ld",sep=""),head=T)
	if(dim(ldfile)[1]>=1){
	R2LEFT[i]=max(as.numeric(ldfile[,7]))
	}
	}

	if(!is.na(CLOSER_LOCI[i,4])){
	shellcommandRIGHT=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/calc_pairwise_ld.sh ",as.character(CLOSER_LOCI[i,1]),as.character(CLOSER_LOCI[i,2]),as.character(CLOSER_LOCI[i,4]),sep=" ")
	system(shellcommandRIGHT)
	ldfile=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/plink_ld.out.",as.character(CLOSER_LOCI[i,1]),".",as.character(min(as.numeric(CLOSER_LOCI[i,2]),as.numeric(CLOSER_LOCI[i,4]))),".",as.character(max(as.numeric(CLOSER_LOCI[i,2]),as.numeric(CLOSER_LOCI[i,4]))),".ld",sep=""),head=T)
	if(dim(ldfile)[1]>=1){
	R2RIGHT[i]=max(as.numeric(ldfile[,7]))
	}
	}
}


RRIGHT=sqrt(R2RIGHT)
RLEFT=sqrt(R2LEFT) #make sure I haven't forgotten one downstream

ZSCORENOVEL=qnorm(1-(as.numeric(NOVEL[,7]))) #this is for our novel hits.
ZSCORELEFT=matrix(ncol=1,nrow=dim(NOVEL)[1],NA)
ZSCORERIGHT=matrix(ncol=1,nrow=dim(NOVEL)[1],NA)
for(i in 1:length(ZSCORELEFT)){
idx=which(known_hits[,1]==CLOSER_LOCI[i,1] & known_hits[,3]==CLOSER_LOCI[i,3] )
if(!identical(idx,integer(0))){
ZSCORELEFT[i]=-qnorm((min(as.numeric(known_hits[idx,6]))))
}

idxr=which(known_hits[,1]==CLOSER_LOCI[i,1] & known_hits[,3]==CLOSER_LOCI[i,4] )
if(!identical(idxr,integer(0))){
ZSCORERIGHT[i]=-qnorm((min(as.numeric(known_hits[idxr,6]))))
}
}

ZSCORECORRECTED=matrix(ncol=1,nrow=dim(NOVEL)[1],NA)
for(i in 1:dim(NOVEL)[1]){
if(!is.na(ZSCORELEFT[i]) & !is.na(ZSCORERIGHT[i])){
ZSCORECORRECTED[i]=ZSCORENOVEL[i]-RLEFT[i]*ZSCORELEFT[i]-RRIGHT[i]*ZSCORERIGHT[i]
}else if(!is.na(ZSCORELEFT[i])){
ZSCORECORRECTED[i]=ZSCORENOVEL[i]-RLEFT[i]*ZSCORELEFT[i]
}else if(!is.na(ZSCORERIGHT[i])){
ZSCORECORRECTED[i]=ZSCORENOVEL[i]-RRIGHT[i]*ZSCORERIGHT[i]
}
}

PVALCORRECTED=2*pnorm(-abs(as.numeric(ZSCORECORRECTED)))

FINAL=cbind(NOVEL,CLOSER_LOCI,ZSCORENOVEL,ZSCORELEFT,ZSCORERIGHT,RRIGHT,RLEFT,ZSCORECORRECTED,PVALCORRECTED)
write.table(FINAL,paste("~lm12/IBD_conditional/NOVEL_variants_",Trait,".corrected.pvalues",sep=""),quote=F,col.names=T,row.names=F,sep="\t")
#END