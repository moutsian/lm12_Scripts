# QQ plot


PVAL=-log10(RES[,TOPLOT])
N <- length(PVAL) ## number of p-values

## create the null distribution 
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(PVAL,null),na.rm=T)


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
qqplot(null,PVAL, ylim=c(0,MAX),xlim=c(0,MAX), pch=19,main=list("chrom20, K11",cex=1.3),
xlab=expression(Theoretical~~-log[10](italic(p))), ylab= expression(Observed~~-log[10](italic(p))))


#END
