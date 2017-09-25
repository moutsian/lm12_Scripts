# add_freq_and_I2 to the list -- note that the update of the function check_if_previously_reported.R does it automatically.
#novel=read.table( "~lm12/IBD_conditional/IIBDGC-meta-analysis-novel-loci-20160211.txt",sep="\t",head=T)
novel=read.table( "~lm12/IBD_conditional/step2.ALL.round2.novel_based_on_info_on_previous_hits_from_5PAPERS.txt",sep="\t",head=T)
novel_plus=matrix(nrow=dim(novel)[1],ncol=dim(novel)[2]+2,"")
novel_plus[,1:dim(novel)[2]]=as.matrix(novel)
colnames(novel_plus)=c(colnames(novel),"FREQ_in_1KG_CEU_GBR","I2")
for(i in 1:dim(novel)[1]){
shellcommand=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_freq_and_I2.sh ",as.character(novel[i,1]),as.character(novel[i,2]),as.character(novel[i,4]),sep=" ")
shellout=unlist(strsplit(system(shellcommand,intern=T)," "))
if(length(shellout==8)){
novel_plus[i,dim(novel)[2]+1]=shellout[7]
novel_plus[i,dim(novel)[2]+2]=shellout[8]
}else{
print("Something went wrong!")
}
}
write.table(novel_plus,"~lm12/IBD_conditional/Step2.ALL.round2.novel_based_on_info_on_previous_hits_from_5PAPERS.txt.plus",quote=F,sep="\t",row.names=F,col.names=T)

#END