pca=read.table("TILI_merged.qc3.pruned.corrected.pca.evec",head=F,stringsAsFactors=F)
fam=read.table("TILI_merged.qc3.pruned.corrected.fam",stringsAsFactors=F)
pcafam=cbind(pca,fam[,c(1,2,6)])
pca_cases_idx=which(fam[,6]==2)
pca_ctrls_idx=which(fam[,6]==1)

plot(pca[,1+1],pca[,2+1],pch=19,xlab="PC1",ylab="PC2",main=list("PCA 1 vs 2\nTILI_merged.qc3.pruned.filtered",cex=1.1))
points(pca[pca_cases_idx,1+1],pca[pca_cases_idx,2+1],pch=19,col="darkred")
points(pca[pca_ctrls_idx,1+1],pca[pca_ctrls_idx,2+1],pch=19,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)

 toremove1=which(abs(pca[,1+1])>0.15)
 toremove2=which(abs(pca[,2+1])>0.15)
 toremove=unique(c(toremove1,toremove2))
 pcafilt=pcafam[-toremove,]
 pcafilt_cases=pcafilt[which(pcafilt[,15]==2),]
 pcafilt_ctrls=pcafilt[which(pcafilt[,15]==1),]
 
plot(pcafilt[,1+1],pcafilt[,2+1],pch=19,xlab="PC1",ylab="PC2",main=list(paste("PC1 vs PC2, after removal of ",length(toremove)," samples",sep=""),cex=1.1))
points(pcafilt_cases[,1+1],pcafilt_cases[,2+1],pch=19,col="darkred")
points(pcafilt_ctrls[,1+1],pcafilt_ctrls[,2+1],pch=19,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)
abline(h=-0.15,col="darkgrey",lty=2)
abline(v=-0.15,col="darkgrey",lty=2)
abline(h=0.075,col="darkgrey",lty=2)
abline(v=0.075,col="darkgrey",lty=2)
abline(h=-0.075,col="darkgrey",lty=2)
abline(v=-0.075,col="darkgrey",lty=2)

 torem_limit=0.95
 pc1=7
 pc2=8
 toremove1=which(abs(pca[,2])>torem_limit)
 toremove2=which(abs(pca[,3])>torem_limit)
 toremove=unique(c(toremove1,toremove2))
 if(length(toremove)>0){
 pcafilt=pcafam[-toremove,]
 }else{
 pcafilt=pcafam
 }
 pcafilt_cases=pcafilt[which(pcafilt[,15]==2),]
 pcafilt_ctrls=pcafilt[which(pcafilt[,15]==1),]
 
 
plot(pcafilt[,pc1+1],pcafilt[,pc2+1],pch=19,xlab=paste("PC",pc1,sep=""),ylab=paste("PC",pc2,sep=""),main=list(paste("PC",pc1," vs PC",pc2," after removal of ",length(toremove)," samples",sep=""),cex=1.1))
points(pcafilt_cases[,pc1+1],pcafilt_cases[,pc2+1],pch=19,col="darkred")
points(pcafilt_ctrls[,pc1+1],pcafilt_ctrls[,pc2+1],pch=19,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)
abline(h=-0.15,col="darkgrey",lty=2)
abline(v=-0.15,col="darkgrey",lty=2)
abline(h=0.075,col="darkgrey",lty=2)
abline(v=0.075,col="darkgrey",lty=2)
abline(h=-0.075,col="darkgrey",lty=2)
abline(v=-0.075,col="darkgrey",lty=2)
#following will just colour the sanger samples in a different colour, to check pc5-pc6 plot.
sanger=grepl("_",pcafilt[,14])
points(pcafilt[sanger,pc1+1],pcafilt[sanger,pc2+1],pch=19,col="orange")




#END