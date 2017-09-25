##this is the same as pca_TILI.R but where a subset of the samples have been removed (through smartpca, iteratively)
#pca=read.table("TILI_merged.qc3.pruned.corrected.pca.outrem.30.4.4.evec",head=F,stringsAsFactors=F)
pca=read.table("TILI_merged.qc3.pruned.corrected.pca.30.10.6.evec",head=F,stringsAsFactors=F)
fam=read.table("TILI_merged.qc3.pruned.corrected.fam",stringsAsFactors=F)

 pc1=3
 pc2=4

postqc_samples=unlist(strsplit(pca[,1],":"))[seq(1,2*dim(pca)[1],2)]

pca_cases=fam[which(fam[,6]==2),1]
pca_ctrls=fam[which(fam[,6]==1),1]

pca_cases_idx=which(postqc_samples%in%pca_cases)
pca_ctrls_idx=which(postqc_samples%in%pca_ctrls)

plot(pca[,pc1+1],pca[,pc2+1],pch=19,xlab=paste("PC",pc1,sep=""),ylab=paste("PC",pc2,sep=""),main=list(paste("PC",pc1," vs PC",pc2,sep=""),cex=1.1))
points(pca[pca_cases_idx,pc1+1],pca[pca_cases_idx,pc2+1],pch=19,col="darkred")
points(pca[pca_ctrls_idx,pc1+1],pca[pca_ctrls_idx,pc2+1],pch=19,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)
abline(h=-0.15,col="darkgrey",lty=2)
abline(v=-0.15,col="darkgrey",lty=2)
abline(h=0.075,col="darkgrey",lty=2)
abline(v=0.075,col="darkgrey",lty=2)
abline(h=-0.075,col="darkgrey",lty=2)
abline(v=-0.075,col="darkgrey",lty=2)




#END