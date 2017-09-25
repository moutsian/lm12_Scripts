ISR21121=read.table("ISR_21121_sorted.transcript_to_gene.counts",head=T)
IU21364=read.table("IU_21364_sorted.transcript_to_gene.counts",head=T)
ISR22219=read.table("ISR_22219_sorted.transcript_to_gene.counts",head=T)

CEIL=5000
IU21364[IU21364>CEIL]=CEIL
ISR21121[ISR21121>CEIL]=CEIL
ISR22219[ISR22219>CEIL]=CEIL

plot(ISR21121[,1],ISR21121[,2],pch=19,xlab="2.5ML ORIG",ylab="2.5ML MOD",main=list("first donor",cex=1.1))
plot(ISR21121[,1],ISR21121[,3],pch=19,xlab="2.5ML ORIG",ylab="1ML MOD",main=list("first donor",cex=1.1))
plot(ISR21121[,1],ISR21121[,4],pch=19,xlab="2.5ML ORIG",ylab="0.5ML MOD",main=list("first donor",cex=1.1))
plot(ISR21121[,2],ISR21121[,3],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("first donor",cex=1.1))
plot(ISR21121[,2],ISR21121[,4],pch=19,xlab="2.5ML MOD",ylab="0.5ML MOD",main=list("first donor",cex=1.1))
plot(ISR21121[,3],ISR21121[,4],pch=19,xlab="1ML MOD",ylab="0.5ML MOD",main=list("first donor",cex=1.1))

plot(ISR21121[,30],ISR21121[,24],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("third donor",cex=1.1))
plot(ISR21121[,27],ISR21121[,21],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("fourth donor",cex=1.1))


plot(IU21364[,1],IU21364[,2],pch=19,xlab="2.5ML ORIG",ylab="2.5ML MOD",main=list("first donor",cex=1.1))
plot(IU21364[,1],IU21364[,3],pch=19,xlab="2.5ML ORIG",ylab="1ML MOD",main=list("first donor",cex=1.1))
plot(IU21364[,1],IU21364[,4],pch=19,xlab="2.5ML ORIG",ylab="0.5ML MOD",main=list("first donor",cex=1.1))
plot(IU21364[,2],IU21364[,3],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("first donor",cex=1.1))
plot(IU21364[,2],IU21364[,4],pch=19,xlab="2.5ML MOD",ylab="0.5ML MOD",main=list("first donor",cex=1.1))
plot(IU21364[,3],IU21364[,4],pch=19,xlab="1ML MOD",ylab="0.5ML MOD",main=list("first donor",cex=1.1))


plot(ISR22219[,3],ISR22219[,1],pch=19,xlab="2.5ML MOD",ylab="0.5ML MOD",main=list("Sample 893, kapa",cex=1.1))
plot(ISR22219[,3],ISR22219[,2],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("Sample 893, kapa",cex=1.1))
plot(ISR22219[,2],ISR22219[,1],pch=19,xlab="1ML MOD",ylab="0.5ML MOD",main=list("Sample 893, kapa",cex=1.1))


plot(ISR22219[,6],ISR22219[,4],pch=19,xlab="2.5ML MOD",ylab="0.5ML MOD",main=list("Sample 908, kapa",cex=1.1))
plot(ISR22219[,6],ISR22219[,5],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("Sample 908, kapa",cex=1.1))
plot(ISR22219[,5],ISR22219[,4],pch=19,xlab="1ML MOD",ylab="0.5ML MOD",main=list("Sample 908, kapa",cex=1.1))

plot(ISR22219[,9],ISR22219[,7],pch=19,xlab="2.5ML MOD",ylab="0.5ML MOD",main=list("Sample 903, kapa",cex=1.1))
plot(ISR22219[,9],ISR22219[,8],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("Sample 903, kapa",cex=1.1))
plot(ISR22219[,8],ISR22219[,7],pch=19,xlab="1ML MOD",ylab="0.5ML MOD",main=list("Sample 903, kapa",cex=1.1))

plot(ISR22219[,12],ISR22219[,10],pch=19,xlab="2.5ML MOD",ylab="0.5ML MOD",main=list("Sample 899, kapa",cex=1.1))
plot(ISR22219[,12],ISR22219[,11],pch=19,xlab="2.5ML MOD",ylab="1ML MOD",main=list("Sample 899, kapa",cex=1.1))
plot(ISR22219[,11],ISR22219[,10],pch=19,xlab="1ML MOD",ylab="0.5ML MOD",main=list("Sample 899, kapa",cex=1.1))


plot(ISR22219[,2],ISR22219[,4],pch=19,xlab="2.5ML MOD",ylab="0.5ML MOD",main=list("first donor",cex=1.1))
plot(ISR22219[,3],ISR22219[,4],pch=19,xlab="1ML MOD",ylab="0.5ML MOD",main=list("first donor",cex=1.1))



CEIL=1000
IU21364[IU21364>CEIL]=CEIL
ISR21121[ISR21121>CEIL]=CEIL
highISR=which(ISR21121[,1]>100 & ISR21121[,5]>100) #2.5Original

plot(ISR21121[highISR,1],ISR21121[highISR,2],pch=19,xlab="2.5ML ORIG",ylab="2.5ML MOD",main=list("first donor",cex=1.1))
abline(v=100,lty=2,col="darkred")
abline(h=100,lty=2,col="darkred")

plot(ISR21121[highISR,1],ISR21121[highISR,3],pch=19,xlab="2.5ML ORIG",ylab="1ML MOD",main=list("first donor",cex=1.1))
abline(v=100,lty=2,col="darkred")
abline(h=100,lty=2,col="darkred")
abline(h=50,lty=2,col="purple")
abline(h=25,lty=2,col="gold")

plot(ISR21121[highISR,1],ISR21121[highISR,4],pch=19,xlab="2.5ML ORIG",ylab="0.5ML MOD",main=list("first donor",cex=1.1))
abline(v=100,lty=2,col="darkred")
abline(h=100,lty=2,col="darkred")
abline(h=50,lty=2,col="purple")
abline(h=25,lty=2,col="gold")

hist(ISR21121[highISR,4],breaks=40,main="Histogram for counts of 0.5ML sample\n for highly expressed (>100counts) genes in 2.5MLORIG")
abline(v=100,lty=2,col="darkred",lwd=2)
abline(v=50,lty=2,col="purple")
abline(v=25,lty=2,col="gold")

hist(ISR21121[highISR,3],breaks=40,main="Histogram for counts of 1ML sample\n for highly expressed (>100counts) genes in 2.5MLORIG")
abline(v=100,lty=2,col="darkred",lwd=2)
abline(v=50,lty=2,col="purple")
abline(v=25,lty=2,col="gold")


highIU=which(IU21364[,1]>100 & IU21364[,5]>100) #2.5Original

plot(IU21364[highIU,1],IU21364[highIU,3],pch=19,xlab="2.5ML ORIG",ylab="1ML MOD",main=list("first donor",cex=1.1))
abline(v=100,lty=2,col="darkred")
abline(h=100,lty=2,col="darkred")
abline(h=50,lty=2,col="purple")
abline(h=25,lty=2,col="gold")

plot(IU21364[highIU,1],IU21364[highIU,4],pch=19,xlab="2.5ML ORIG",ylab="0.5ML MOD",main=list("first donor",cex=1.1))
abline(v=100,lty=2,col="darkred")
abline(h=100,lty=2,col="darkred")
abline(h=50,lty=2,col="purple")
abline(h=25,lty=2,col="gold")

hist(IU21364[highIU,4],breaks=40,main="Histogram for counts of 0.5ML sample\n for highly expressed (>100counts) genes in 2.5MLORIG")
abline(v=100,lty=2,col="darkred",lwd=2)
abline(v=50,lty=2,col="purple")
abline(v=25,lty=2,col="gold")

hist(IU21364[highIU,3],breaks=40,main="Histogram for counts of 1ML sample\n for highly expressed (>100counts) genes in 2.5MLORIG")
abline(v=100,lty=2,col="darkred",lwd=2)
abline(v=50,lty=2,col="purple")
abline(v=25,lty=2,col="gold")

COR22219=matrix(ncol=12,nrow=12,0)
for(i in 1:12){for (j in 1:12){
COR22219[i,j]=cor(ISR22219[,i],ISR22219[,j])
}
}

corISR22219=matrix(ncol=40,nrow=1,0)
for (i in 1:40){ corISR22219[1,i]=cor(ISR21121[,4],ISR22219[,i])}
colnames(corISR22219)=colnames(ISR22219)
corISR22219
corISR22219[1,which.max(corISR22219)]
#END






