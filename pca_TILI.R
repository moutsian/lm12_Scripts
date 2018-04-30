pca=read.table("tili.i3.site_QCplus_binom.sample_QC.annot.pruned.30.evec",head=F,stringsAsFactors=F)
fam=read.table("tili.i3.site_QCplus_binom.sample_QC.annot.pruned.fam",stringsAsFactors=F)
pcafam=cbind(pca,fam[,c(1,2,6)])
pca_cases_idx=which(fam[,6]==1)
pca_ctrls_idx=which(fam[,6]==2)

PC1=7
PC2=8
TITLE=paste("PC ",PC1," vs ",PC2,"\nTILI exomeseq with all ctrls, post QC",sep="")
plot(pca[,PC1+1],pca[,PC2+1],pch=19,xlab=paste("PC",PC1,sep=""),ylab=paste("PC",PC2,sep=""),main=list(TITLE,cex=1.1))
points(pca[pca_cases_idx,PC1+1],pca[pca_cases_idx,PC2+1],pch=19,col="darkred")
points(pca[pca_ctrls_idx,PC1+1],pca[pca_ctrls_idx,PC2+1],pch=19,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)

 toremove1=which(abs(pca[,1+1])>0.05)
 toremove2=which(abs(pca[,2+1])>0.05)
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
 pc1=3
 pc2=4
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