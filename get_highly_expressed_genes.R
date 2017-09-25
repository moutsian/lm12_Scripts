

#following is gene quantification data from kallisto - estimated counts
data=read.table("/lustre/scratch115/projects/pneumo-inf/TranscriptCountFromfastqFiles.kallisto/kallisto.100Bootstrap.GeneQuantification/AllSamplesGeneQuant",head=T)
info=read.table("/nfs/users/nfs_l/lm12/RNAextractions.PneumoRNAseqPilot.txt",head=T,sep="\t")
gene=as.matrix(as.character(data[,1]))
data=data[,-1]
datablood=data[,which(info$Tissue=="Blood")]
datanasal=data[,which(info$Tissue=="Nasal")]
databloodnorm=datablood
for(i in 1:dim(datablood)[2]){
databloodnorm[,i]=datablood[,i]-mean(datablood[,i])
}
meanbefore=apply(datablood,2,mean)
meanafter=apply(databloodnorm,2,mean)
normasums=apply(databloodnorm,1,sum)
idx100=which(normasums>1000000)
idx10=which(normasums>4500000)
gene[idx10,1]

#repeat for only blood -5 samples
databloodmin5=data[,which(info$Tissue=="Blood" & info$Condition=="Minus5")]
databloodmin5norm=databloodmin5
for(i in 1:dim(databloodmin5)[2]){
databloodmin5norm[,i]=databloodmin5[,i]-mean(databloodmin5[,i])
}
meanbefore=apply(databloodmin5,2,mean)
meanafter=apply(databloodmin5norm,2,mean)
normasums_min5=apply(databloodmin5norm,1,sum)
idxm5_100=which(normasums_min5>500000)
idxm5_10=which(normasums_min5>2300000)
gene[idxm5_10,1]


#repeat for nasal samples
datanasalnorm=datanasal
for(i in 1:dim(datanasal)[2]){
datanasalnorm[,i]=datanasal[,i]-mean(datanasal[,i])
}
meanbefore=apply(datanasal,2,mean)
meanafter=apply(datanasalnorm,2,mean)
normasumsn=apply(datanasalnorm,1,sum)
idx100n=which(normasumsn>450000)
idx10n=which(normasumsn>2300000)
gene[idx10n,1]

#END