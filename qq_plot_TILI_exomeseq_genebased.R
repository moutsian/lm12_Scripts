# QQ plot and Manhattan plot

#epactsres=read.table("C:/Academic/SANGER/TILI/exomeseq/test.gene.allchr.fc.skato.epacts",head=F,stringsAsFactors=F) 
#epactsres=read.table("C:/Academic/SANGER/TILI/exomeseq/test.gene.allchr.PD.lenient.0.01.withcov.skato.epacts",head=F,stringsAsFactors=F) 
#epactsres=read.table("C:/Academic/SANGER/TILI/exomeseq/tili.i2.qc3binom.hepmix_vs_ctrls.allchr.PD.lenient.withcov.0.01.skato.epacts",head=F,stringsAsFactors=F)

annofilter="PD.strict"
test="emmaxVT" # "emmaxVT" or "skato.withcov" or "skato.nocov"
epactsres=NULL#
#tili.i3.allchr.PD.strict.0.01.emmaxVT.epacts
epactsres=read.table(paste("C:/Academic/SANGER/TILI/exomeseq/tili.i3.allchr.",annofilter,".0.01.",test,".epacts",sep=""),head=F,stringsAsFactors=F)
epactsres_filt=NULL
if(test=="emmaxVT"){
epactsres_filt=epactsres[which(as.numeric(as.character(epactsres[,6]))>1),] #this is for the skat-o output, for the emmax you need to change the index as the respective column is 6

}else{
epactsres_filt=epactsres[which(as.numeric(as.character(epactsres[,8]))>1),] #this is for the skat-o output, for the emmax you need to change the index as the respective column is 6
}
PVAL=-log10(as.numeric(as.character(epactsres_filt[,10]))) 
N <- length(PVAL) ## number of p-values

## create the null distribution 
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(PVAL,null),na.rm=T)
#MAX=15

## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)

## the jth order statistic from a 
## uniform(0,1) sample 
## has a beta(j,n-j+1) distribution 
## (Casella & Berger, 2002, 
## 2nd edition, pg 230, Duxbury)

for(i in 1:N){
c95[i] <- qbeta(0.95,i,N-i+1)
c05[i] <- qbeta(0.05,i,N-i+1)
}

## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",lty=2,
axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", lty=2,
axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red",lwd=2)
par(new=T)

## add the qqplot
qqplot(null,PVAL, ylim=c(0,MAX),xlim=c(0,MAX), pch=19,main=list(paste("TILI exome seq, epacts ",test," results, with all ctrls (i3)\n(",annofilter,", 1% MAF cutoff)",sep=""),cex=1.3),
xlab=expression(Theoretical~~-log[10](italic(p))), ylab= expression(Observed~~-log[10](italic(p))))

##plot with only REFSEQ sequences
epactsres_refseq=epactsres_filt[grepl("ENSG",epactsres_filt[,4]),]
#need to filter on deviations from hwe in cases too, as this filters out some weird variants
PVAL=-log10(as.numeric(as.character(epactsres_refseq[,10]))) 
N <- length(PVAL) ## number of p-values

## create the null distribution 
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(PVAL,null),na.rm=T)
#MAX=15

## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)

## the jth order statistic from a 
## uniform(0,1) sample 
## has a beta(j,n-j+1) distribution 
## (Casella & Berger, 2002, 
## 2nd edition, pg 230, Duxbury)

for(i in 1:N){
c95[i] <- qbeta(0.95,i,N-i+1)
c05[i] <- qbeta(0.05,i,N-i+1)
}

## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",lty=2,
axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", lty=2,
axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red",lwd=2)
par(new=T)

## add the qqplot
qqplot(null,PVAL, ylim=c(0,MAX),xlim=c(0,MAX), pch=19,main=list(paste("TILI exome seq, epacts ",test," results, with all ctrls (i3) \n(",annofilter,", 1% REFSEQ only, MAF cutoff)",sep=""),cex=1.3),
xlab=expression(Theoretical~~-log[10](italic(p))), ylab= expression(Observed~~-log[10](italic(p))))


#END
