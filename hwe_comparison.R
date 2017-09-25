setwd("U:/R_workspace/")

GAII_CASES_SD = read.table("../data/high_load.GAII.cases.hwe",head=T)
GAII_CTRLS_SD = read.table("../data/high_load.GAII.ctrls.hwe",head=T)
HiSEQ_CASES_SD = read.table("../data/high_load.HiSEQ.cases.hwe",head=T)
HiSEQ_CTRLS_SD = read.table("../data/high_load.HiSEQ.ctrls.hwe",head=T)

NOHL_GAII_CASES_SD = read.table("../data/no_high_load.GAII.cases.hwe",head=F,skip=1)
NOHL_GAII_CTRLS_SD = read.table("../data/no_high_load.GAII.ctrls.hwe",head=F,skip=1)
NOHL_HiSEQ_CASES_SD = read.table("../data/no_high_load.HiSEQ.cases.hwe",head=F,skip=1)
NOHL_HiSEQ_CTRLS_SD = read.table("../data/no_high_load.HiSEQ.ctrls.hwe",head=F,skip=1)


par(mfrow=c(2,4))
hist(breaks=40,HiSEQ_CASES_SD[,6],pch=19,col="darkred",main=list("p-value of dev from HWE, HiSEQ cases (high load SNPs)",cex=1.2))
hist(breaks=40,HiSEQ_CTRLS_SD[,6],pch=19,col="darkblue",main=list("p-value of dev from HWE, HiSEQ ctrls (high load SNPs)",cex=1.2))
hist(breaks=40,GAII_CASES_SD[,6],pch=19,col="darkred",main=list("p-value of dev from HWE, GAII cases (high load SNPs)",cex=1.2))
hist(breaks=40,GAII_CTRLS_SD[,6],pch=19,col="darkblue",main=list("p-value of dev from HWE, GAII ctrls (high load SNPs)",cex=1.2))

hist(breaks=40,NOHL_HiSEQ_CASES_SD[,6],pch=19,col="darkred",main=list("p-value of dev from HWE, HiSEQ cases (rest)",cex=1.2))
hist(breaks=40,NOHL_HiSEQ_CTRLS_SD[,6],pch=19,col="darkblue",main=list("p-value of dev from HWE, HiSEQ ctrls (rest)",cex=1.2))
hist(breaks=40,NOHL_GAII_CASES_SD[,6],pch=19,col="darkred",main=list("p-value of dev from HWE, GAII cases (rest)",cex=1.2))
hist(breaks=40,NOHL_GAII_CTRLS_SD[,6],pch=19,col="darkblue",main=list("p-value of dev from HWE, GAII ctrls (rest)",cex=1.2))



#check the same for the rest.


NOHL_GAII_CASES_SD = read.table("../data/no_high_load.GAII.cases.ldepth",head=T)
NOHL_GAII_CTRLS_SD = read.table("../data/no_high_load.GAII.ctrls.ldepth",head=T)
NOHL_HiSEQ_SD = read.table("../data/no_high_load.HiSEQ.ldepth",head=T)



MIN_PVAL = apply(cbind(HiSEQ_CTRLS_SD[,6],GAII_CTRLS_SD[,6]),1,min)




#