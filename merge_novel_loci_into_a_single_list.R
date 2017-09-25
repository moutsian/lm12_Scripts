#final merging for eye checking

#shellcommand=paste(" cat ~lm12/IBD_conditional/step2.*.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.txt >  ~lm12/IBD_conditional/step2.ALL.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.txt",sep="")

##********* CAREFUL not to include entries from dataset B . Some renaming is probably needed!! *********

shellcommand=paste(" cat ~lm12/IBD_conditional/step2.*NOVEL.round6.with_info_on_previous_hits_from_all_PAPERS.txt >  ~lm12/IBD_conditional/step2.ALL.NOVEL.round6.with_info_on_previous_hits_from_all_PAPERS.txt",sep="")

system(shellcommand)
#results=as.matrix(read.table("~lm12/IBD_conditional/step2.ALL.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.txt",sep="\t",head=T))
#results=as.matrix(read.table("~lm12/IBD_conditional/step2.ALL.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.with_corr_pvals.txt",sep="\t",head=T))
results=as.matrix(read.table("~lm12/IBD_conditional/step2.ALL.NOVEL.round6.with_info_on_previous_hits_from_all_PAPERS.txt",sep="\t",head=T))

results=results[results[,1]!="Trait",]
resultso=results[order(as.numeric(results[,2]),as.numeric(results[,4])),]
resultso_original=resultso

dataout=NULL
n_lead_var=1;
more_info=matrix(nrow=dim(resultso)[1],ncol=3,"") # trait, position, p-val
END=dim(resultso)[1]
for(i in 2:END){
if (i> END){
break}
if(as.numeric(resultso[i,2])==as.numeric(resultso[i-1,2]) & as.numeric(resultso[i,8])<=as.numeric(resultso[i-1,9])){
#merge and keep strongest p-value.

if(as.numeric(resultso[i,7])<as.numeric(resultso[i-1,7])){
#keep results from i
more_info[i,1]=paste(more_info[i,1],resultso[i-1,1],sep=",")
more_info[i,2]=paste(more_info[i,2],resultso[i-1,4],sep=",")
more_info[i,3]=paste(more_info[i,3],resultso[i-1,7],sep=",")
resultso[i,8]=min(as.numeric(resultso[i,8]),as.numeric(resultso[i-1,8]))
resultso[i,9]=max(as.numeric(resultso[i,9]),as.numeric(resultso[i-1,9]))

resultso=resultso[-(i-1),]
more_info=more_info[-(i-1),]
END=END-1
}else{
#keep results from i-1 
more_info[i-1,1]=paste(more_info[i-1,1],resultso[i,1],sep=",")
more_info[i-1,2]=paste(more_info[i-1,2],resultso[i,4],sep=",")
more_info[i-1,3]=paste(more_info[i-1,3],resultso[i,7],sep=",")
resultso[i-1,8]=min(as.numeric(resultso[i,8]),as.numeric(resultso[i-1,8]))
resultso[i-1,9]=max(as.numeric(resultso[i,9]),as.numeric(resultso[i-1,9]))

resultso=resultso[-i,]
more_info=more_info[-i,]
END=END-1
}

}
}
resultso_plus=cbind(resultso,more_info)
colnames(resultso_plus)=c(colnames(resultso),"OtherTrait","position","p-value")
write.table(resultso_plus,"~lm12/IBD_conditional/step2.ALL.NOVEL.round6.with_info_on_previous_hits_from_all_PAPERS.auto.txt",quote=F,col.names=T,row.names=F,sep="\t")

#END