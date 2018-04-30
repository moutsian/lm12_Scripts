# QQ plot and Manhattan plot

#res=read.table("TILI_merged.qc4.gemma_output.assoc.txt",head=T,stringsAsFactors=F)
#res=read.table("TILI_merged.qc4_v2.gemma_output.assoc.txt",head=T,stringsAsFactors=F)
#snptestres=read.table("TILI_merged.qc4.snptest_output.txt",head=T,stringsAsFactors=F)
#snptestres=read.table("TILI_merged.qc4_v2.snptest_output.txt",head=T,stringsAsFactors=F) #this is the default one
#snptestres=read.table("TILI_merged.qc4_v2.nospecial.snptest_output.txt",head=T,stringsAsFactors=F) #this is the one without the individuals which were tolerant to 6MP when rechallenged
#snptestres=read.table("TILI_merged.qc4_v2.definite.snptest_output.txt",head=T,stringsAsFactors=F) #this is the one with only definite TILI cases.
#snptestres=read.table("TILI_merged.qc4_v2.hepatocellular.snptest_output.txt",head=T,stringsAsFactors=F) #this is the one with only hepatocellular TILI cases.
snptestres=read.table("TILI_merged.qc1e.snptest_output.oct.txt",head=T,stringsAsFactors=F) #this is with the rerun of the QC, with added differential missingness and control control association test checks
#snptestres=read.table("TILI_merged.qc1e.hepmix.snptest_output.txt",head=T,stringsAsFactors=F) #this is with the rerun of the QC, with added differential missingness and control control association test checks
#PVAL=-log10(as.numeric(as.character(res[,9]))) #from GEMMA
PVAL=-log10(as.numeric(as.character(snptestres[,42]))) #from SNPTEST2
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
qqplot(null,PVAL, ylim=c(0,MAX),xlim=c(0,MAX), pch=19,main=list("TILI, SNPTEST2 results, hepmix vs ctrls\n (4 PCs as covariates)",cex=1.3),
xlab=expression(Theoretical~~-log[10](italic(p))), ylab= expression(Observed~~-log[10](italic(p))))

##genomic control
chisq1 <- qchisq(1-as.numeric(as.character(res[,9])),1)
lambda1 = median(chisq1)/qchisq(0.5,1)

chisq2 <- qchisq(1-as.numeric(as.character(snptestres[,42])),1)
lambda2 = median(chisq2,na.rm=T)/qchisq(0.5,1)

#Not that GenABEL gives higher lambda estimates
library(GenABEL)
lambda_snptest  = estlambda(as.numeric(as.character(snptestres[,42])))
lambda_gemma= estlambda(as.numeric(as.character(res[,9])))

# ##check some manhattan plots per chromosome
 # plotchrom = function(chrom,YLIM){
 # chrom=chrom
 # par(mfrow=c(2,1))
 # plot(snptestres[which(snptestres[,1]==chrom),4],-log10(snptestres[which(snptestres[,1]==chrom),42]),ylab="-log10P",xlab=paste("chrom ",chrom,sep=""),main="SNPTEST2 with 5 PCs",pch=20,ylim=c(0,YLIM))
 # abline(h=7.3,lty=2,col="red")
 # plot(res[which(res[,1]==chrom),3],-log10(res[which(res[,1]==chrom),9]),ylab="-log10P",xlab=paste("chrom ",chrom,sep=""),main="GEMMA",pch=20,ylim=c(0,YLIM))
 # abline(h=7.3,lty=2,col="red")
# }
# plotchrom(17,10)


library(qqman)
snptestres_filtered=snptestres[-which(is.na(snptestres$frequentist_add_pvalue)),]
manhattan(snptestres_filtered,chr="alternate_ids",bp="position",p="frequentist_add_pvalue",snp="rsid")
snptestres_filtered[which(snptestres_filtered$frequentist_add_pvalue<=1e-15),42]=1e-15
par(mfrow=c(1,1))
manhattan(snptestres_filtered,chr="alternate_ids",bp="position",p="frequentist_add_pvalue",snp="rsid")

#have a look at some specific hits
snptestres_filtered[which(snptestres_filtered[,42]<1e-05),c(1:2,4,30:31,42:44)]
write.table(snptestres_filtered[which(snptestres_filtered[,42]<1e-05),c(1:2,4,30:31,42:44)],"putative_hits.txt",sep="\t",quote=F,col.names=T,row.names=F)


#END
