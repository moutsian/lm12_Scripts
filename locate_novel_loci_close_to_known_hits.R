
#All hits:
IBD=read.table("~lm12/IBD_conditional/step2.IBD.all.round2.with_info_on_previous_hits_from_5PAPERS.txt",head=T,sep="\t")
CD=read.table("~lm12/IBD_conditional/step2.CD.all.round2.with_info_on_previous_hits_from_5PAPERS.txt",head=T,sep="\t")
UC=read.table("~lm12/IBD_conditional/step2.UC.all.round2.with_info_on_previous_hits_from_5PAPERS.txt",head=T,sep="\t")

known_hits=rbind(IBD,CD,UC)
known_hits=known_hits[known_hits[,16]!="",]

#Novel hits:

Trait="IBD"
NOVEL=read.table(paste("~lm12/IBD_conditional/step2.",Trait,".NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.txt",sep=""),head=T,sep="\t")

MARGIN=500000 #kb
for_gcta=NULL
for(i in 1:dim(NOVEL)[1]){
idx= which(as.numeric(NOVEL[i,2])==as.numeric(known_hits[,1]) & ( abs(as.numeric(known_hits[,8])-as.numeric(NOVEL[i,8]))<500000 | abs(as.numeric(known_hits[,7])-as.numeric(NOVEL[i,9]))<500000 )  )
if(length(idx)>0){
all_variants=""
for(j in 1:length(idx)){
all_variants=paste(all_variants,as.character(known_hits[idx[j],2]),sep=",")
}
print(paste("i= ",i,". For ",NOVEL[i,2],":",NOVEL[i,4],", check ",known_hits[idx,2],sep=""))
tmp=unique(unlist(strsplit(all_variants,",")))
tmp=tmp[2:length(tmp)]
for_gcta=rbind(for_gcta,c(Trait,as.numeric(NOVEL[i,2]),as.numeric(NOVEL[i,4]),as.numeric(NOVEL[i,8]),as.numeric(NOVEL[i,9]),as.numeric(known_hits[idx[1],3]),as.numeric(known_hits[idx[1],7]),as.numeric(known_hits[idx[1],8]),paste(tmp, sep="", collapse=",")))
}
}
colnames(for_gcta)=c("TRAIT","CHR","Novel_Locus","Novel_Locus_WindowLeft","Novel_Locus_WindowRight","Main_nearby_known_variant","Known_Locus_Window_Left","Known_Locus_Window_Right","Other_gw_nearby_variants")
write.table(for_gcta,paste("~lm12/IBD_conditional/Loci_to_consider_for_conditional_analysis.",Trait,".txt",sep=""),sep="\t",quote=F,col.names=T,row.names=F)
#END