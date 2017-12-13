# QQ plot and Manhattan plot

res=read.table("C:/Academic/SANGER/TILI/exomeseq/tili.i2.qc3.cholmix_vs_ctrls.gemma_output.assoc.txt",head=T,stringsAsFactors=F)
#res=read.table("C:/Academic/SANGER/TILI/exomeseq/tili.i2.qc.sample_qc3.rsID.annot.gemma_output.assoc.txt",head=T,stringsAsFactors=F)
#res=read.table("C:/Academic/SANGER/TILI/exomeseq/tili.i2.qc.sample_qc3.output.txt",head=T,stringsAsFactors=F) 
snptestres=read.table("C:/Academic/SANGER/TILI/exomeseq/tili.i2.qc3.cholmix_vs_ctrls.output.txt",head=T,stringsAsFactors=F) #this is with the rerun of the QC, with added differential missingness and control control association test checks
res[which(res[,13]==0),13]=1e-100
snptestres_filt=snptestres[which(snptestres[,31]>=0.005),]
hwe=read.table("C:/Academic/SANGER/TILI/exomeseq/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.hwe",head=T,stringsAsFactors=F)
snptestmerged=merge(snptestres_filt,hwe,by.x="rsid",by.y="SNP")
snptestres_filt=snptestmerged[which(snptestmerged[,54]>1e-15),]
#need to filter on deviations from hwe in cases too, as this filters out some weird variants
resmerged=merge(res,hwe,by.x="rs",by.y="SNP")
res_filt=resmerged[which(resmerged[,22]>1e-7),] #
res_filt2=merge(res_filt,snptestres,by.x="rs",by.y="rsid") #to filter n MAF>0.005 in both cases and controls
res_filt=res_filt2[which(res_filt2[,52]>=0.005),]
PVAL=-log10(as.numeric(as.character(res_filt[,13]))) #from GEMMA

#PVAL=-log10(as.numeric(as.character(snptestres_filt[,42]))) #from SNPTEST2
N <- length(PVAL) ## number of p-values

## create the null distribution 
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(PVAL,null),na.rm=T)
#MAX=20

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
qqplot(null,PVAL, ylim=c(0,MAX),xlim=c(0,MAX), pch=19,main=list("TILI exome seq, hep + mix vs ctrls,\n GEMMA   ",cex=1.3),
xlab=expression(Theoretical~~-log[10](italic(p))), ylab= expression(Observed~~-log[10](italic(p))))

# ##genomic control
# chisq1 <- qchisq(1-as.numeric(as.character(res[,9])),1)
# lambda1 = median(chisq1)/qchisq(0.5,1)

chisq2 <- qchisq(1-as.numeric(as.character(snptestres_filt[,42])),1)
lambda2 = median(chisq2,na.rm=T)/qchisq(0.5,1)

#Not that GenABEL gives higher lambda estimates
library(GenABEL)
lambda_snptest  = estlambda(as.numeric(as.character(snptestres_filt[,42])))
lambda_gemma= estlambda(as.numeric(as.character(res_filt[,13])))

##check some manhattan plots per chromosome
 plotchrom = function(chrom,YLIM){
 chrom=chrom
 par(mfrow=c(2,1))
 plot(snptestres_filt[which(snptestres_filt[,1]==chrom),4],-log10(snptestres_filt[which(snptestres_filt[,1]==chrom),42]),ylab="-log10P",xlab=paste("chrom ",chrom,sep=""),main="SNPTEST2 with 5 PCs",pch=20,ylim=c(0,YLIM))
 abline(h=7.3,lty=2,col="red")
 # plot(res[which(res[,1]==chrom),3],-log10(res[which(res[,1]==chrom),9]),ylab="-log10P",xlab=paste("chrom ",chrom,sep=""),main="GEMMA",pch=20,ylim=c(0,YLIM))
 # abline(h=7.3,lty=2,col="red")
}
plotchrom(17,10)


library(qqman)
snptestres_filt_toplot=snptestres_filt[!is.na(snptestres_filt$frequentist_add_pvalue),]
manhattan(snptestres_filt_toplot,chr="alternate_ids",bp="position",p="frequentist_add_pvalue",snp="rsid")
#snptestres_filt_toplot[which(snptestres_filt_toplot$frequentist_add_pvalue<=1e-15),42]=1e-15
par(mfrow=c(1,1))
snptestres_filt_auto=snptestres_filt_toplot[which(snptestres_filt_toplot[,2]>=1 & snptestres_filt_toplot[,2]<=22),]
manhattan(snptestres_filt_auto,chr="alternate_ids",bp="position",p="frequentist_add_pvalue",snp="rsid")

check a bit the results
snptestres_filt[which(snptestres_filt[,2]==chr & snptestres_filt[,42]<1e-06),c(1:2,4,30,31,42:44)]


par(mfrow=c(1,1))
TEST=14
res_auto=res_filt[which(res_filt[,2]>=1 & res_filt[,2]<=22),]
res_auto[which(res_auto[,TEST]<1e-15),TEST]=1e-15
manhattan(res_auto,chr="chr",bp="ps",p="p_score",snp="rs")

par(mfrow=c(1,1))
TEST=14
res_auto=res_filt[which(res_filt[,2]>=1 & res_filt[,2]<=22),]
res_auto[which(res_auto[,TEST]<1e-15),TEST]=1e-15
manhattan(res_auto,chr="chr",bp="ps",p="p_score",snp="rs")

chr=15
res_auto[which(res_auto[,2]==chr & res_auto[,13]<1e-06),]

##check concordance between snptest2 and gemma_outputs


#END
