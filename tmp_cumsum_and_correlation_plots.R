plot(tmp[1:500,1],pch=19,ylim=c(5e06,7e07),xlab=list("no of genes (ordered by counts)",cex=1.1),main=list("Dataset 21364, library IU",cex=1.1),ylab=list("Total counts (4 samples per concentration)",cex=1.1))
points(tmp[1:500,2],pch=19,col="darkblue")
points(tmp[1:500,3],pch=19,col="orange")
points(tmp[1:500,4],pch=19,col="yellow")
legend(x=400,y=3e07,colnames(tmp),col=c("black","darkblue","orange","yellow"),pch=19)

plot(tmp[1:50,1],pch=19,ylim=c(5e06,7e07),xlab=list("no of genes (ordered by counts)",cex=1.1),main=list("Dataset 21364, library IU",cex=1.1),ylab=list("Total counts (4 samples per concentration)",cex=1.1))
points(tmp[1:50,2],pch=19,col="darkblue")
points(tmp[1:50,3],pch=19,col="orange")
points(tmp[1:50,4],pch=19,col="yellow")
legend(x=40,y=6.7e07,colnames(tmp),col=c("black","darkblue","orange","yellow"),pch=19)




plot(tmp21121[1:500,1],pch=19,ylim=c(5e06,7e07),xlab=list("no of genes (ordered by counts)",cex=1.1),main=list("Dataset 21121, library IU",cex=1.1),ylab=list("Total counts (4 samples per concentration)",cex=1.1))
points(tmp21121[1:500,2],pch=19,col="darkblue")
points(tmp21121[1:500,3],pch=19,col="orange")
points(tmp21121[1:500,4],pch=19,col="yellow")
legend(x=400,y=3e07,colnames(tmp21121),col=c("black","darkblue","orange","yellow"),pch=19)

plot(tmp21121[1:50,1],pch=19,ylim=c(5e06,7e07),xlab=list("no of genes (ordered by counts)",cex=1.1),main=list("Dataset 21121, library IU",cex=1.1),ylab=list("Total counts (4 samples per concentration)",cex=1.1))
points(tmp21121[1:50,2],pch=19,col="darkblue")
points(tmp21121[1:50,3],pch=19,col="orange")
points(tmp21121[1:50,4],pch=19,col="yellow")
legend(x=40,y=2.5e07,colnames(tmp21121),col=c("black","darkblue","orange","yellow"),pch=19)

cortable21364=matrix(ncol=16,nrow=16,0)
cortable21121=matrix(ncol=16,nrow=16,0)
for(i in 1:16){
for(j in 1:16){
cortable21121[i,j]=cor(iu21121[,i],iu21121[,j])
cortable21364[i,j]=cor(iu21364[,i],iu21364[,j])
}}


###only for genes in refseq

cortable21364refseq=matrix(ncol=16,nrow=16,0)
cortable21121refseq=matrix(ncol=16,nrow=16,0)
for(i in 1:16){
for(j in 1:16){
cortable21121refseq[i,j]=cor(iu21121[index21121,i],iu21121[index21121,j])
cortable21364refseq[i,j]=cor(iu21364[index21364,i],iu21364[index21364,j])
}}

###only for genes with at least two counts per sample
##this gives me almost no difference
new21121=iu21121
new21121[new21121<=2]=NA
new21364=iu21364
new21364[new21364<=2]=NA

cortable21364nz=matrix(ncol=16,nrow=16,0)
cortable21121nz=matrix(ncol=16,nrow=16,0)
for(i in 1:16){
for(j in 1:16){
cortable21121nz[i,j]=cor(new21121[,i],new21121[,j],use="complete.obs")
cortable21364nz[i,j]=cor(new21364[,i],new21364[,j],use="complete.obs")
}}


###only for genes with at least 5 counts per sample
##this gives me almost no difference
new21121=iu21121
new21121[new21121<=5]=NA
new21364=iu21364
new21364[new21364<=5]=NA

cortable21364nz5=matrix(ncol=16,nrow=16,0)
cortable21121nz5=matrix(ncol=16,nrow=16,0)
for(i in 1:16){
for(j in 1:16){
cortable21121nz5[i,j]=cor(new21121[,i],new21121[,j],use="complete.obs")
cortable21364nz5[i,j]=cor(new21364[,i],new21364[,j],use="complete.obs")
}}
#end
#