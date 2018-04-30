# Prioritisation for follow up in our group 
# IBD signals. Checking for variants which are strong candidates for driving the signal we have observed in IBD, then checking whether they are in enhancer regions.

fm=read.table("finemapped_loci_nature2017.csv",sep=",",head=T)
enhancers=read.table("human_permissive_enhancers_phase_1_and_2.bed",stringsAsFactors=F)

vie=matrix(ncol=2,nrow=dim(fm),0)
for(i in 1:dim(fm)[1]){
varpos=as.numeric(as.character(fm$position.lead[i]))
chrom=paste("chr",fm$chr[i],sep="")
idx=which(enhancers[,2]<=varpos & enhancers[,3]>=varpos & enhancers[,1]==chrom)
if(length(idx)>0)
vie[i,]=c(i,idx[1])
}

oven=cbind(fm[vie[,1]>0,],enhancers[vie[,2],])

### now get the TPM for these enhancers into a file (the whole one is too big and unneccessary to load), then load it and have a look at TPM in T-cells

data=read.table("overlapping_enhancers_subset.txt",head=T)
Tcells = c("CNhs11858","CNhs11867","CNhs13195","CNhs13223","CNhs13202","CNhs13215","CNhs13203","CNhs13204","CNhs13205","CNhs13206","CNhs13512","CNhs13513","CNhs13538","CNhs13539","CNhs13813","CNhs13920","CNhs13918","CNhs13811","CNhs13814","CNhs13921","CNhs13919","CNhs13812","CNhs13914","CNhs13915","CNhs13238","CNhs13239","CNhs13235","CNhs13237")
Tcell_idx=which(colnames(data)%in%Tcells)
data[,Tcell_idx]
write.table(data[,c(1,Tcell_idx)],"overlapping_enhancers_Tcell_TPM.txt",quote=F,sep="\t",row.names=F,col.names=T)
#END


###make some histograms of the expression of that enhancer across tissues
enh=read.table("C:/Academic/SANGER/IBD 2017 ONWARDS/TARGET PRIORITISATION/overlapping_enhancers_subset.txt",head=T)


ENH=3
for (ENH in 1:dim(enh)[1]){
Filename=paste(enh[ENH,1],"_histo.png",sep="")
myname=gsub(":", "_", Filename)
png(filename=myname)
hist(as.numeric(as.character(enh[ENH,c(2:dim(enh)[2])])),breaks=40,ylim=c(0,100),main=paste("expression of enhancer ",enh[ENH,1]," across tissues",sep=""),xlab="TPM")
dev.off()
}