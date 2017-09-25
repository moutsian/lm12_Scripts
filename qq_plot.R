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



################################################################################
## COMPARISON OF THE EFFECT ON SNP 								##
################################################################################
MAX=max(c(PRESENCE,ABSENCE),na.rm=T)
MIN=min(c(PRESENCE,ABSENCE),na.rm=T)

PRESENCE=exp(RES[,IDX_MAIN]+RES[,IDX_MAIN+1])
ABSENCE=exp(RES[,IDX_MAIN])
plot(PRESENCE,ABSENCE,pch=19,ylim=c(MIN,MAX),xlim=c(MIN,MAX), ylab=list(paste("SNP Effect in the absence of ",allele,sep=""),cex=1.2),
xlab=list(paste("SNP Effect in the presence of ",allele,sep=""),cex=1.2),
main=list(paste("Effect of the presence/absence of ",allele,"on SNP Effects, ", cohort,sep=""),cex=1.3))
abline(0,1,col="red",lwd=2)
binom.test(sum(PRESENCE>ABSENCE,na.rm=T),length(PRESENCE),alternative="two.sided")





#compare logit and probit (note that the two sets of results for the respective models will need to be saved in 
# RES_LOGIT and RES_PROBIT) :
#plot(-log10(RES_PROBIT[,TOPLOT]),-log10(RES_LOGIT[,TOPLOT]),pch=19,xlab= expression(-log[10](italic(p))~~ "link=probit"),
#ylab= expression(-log[10](italic(p))~~ "link=logit"), main=list(paste("Comparison of the p-values between logit and probit models,\n",cohort,", interaction term with ",allele,sep=""),cex=1.3) )
#abline(0,1,col="red",lwd=2)

#END
