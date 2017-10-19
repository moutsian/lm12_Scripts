#Split HLA files into cases and ctrls


#GWAS1
caseinfo=read.table("/lustre/scratch115/projects/ibdgwas/pre_imputation/qc/GWAS1/CD.inds",stringsAsFactors=F)
ctrlinfo=read.table("/lustre/scratch115/projects/ibdgwas/pre_imputation/qc/GWAS1/CTRL.inds",stringsAsFactors=F)
HLADIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/GWAS1/"
PREFIX="wtccc1_hg19.final.xMHC.HLA_"
SUFFIX=".4D.bestguess.txt"
LOCI=c("A","C","B","DRB1","DQA1","DQB1","DPB1")
for(loc in 1:length(LOCI)){
LOCUS=LOCI[loc]
data=read.table(paste(HLADIR,PREFIX,LOCUS,SUFFIX,sep=""),head=T,stringsAsFactors=F)
casetmp=paste(caseinfo[,1],"_",caseinfo[,2],sep="")
ctrltmp=paste(ctrlinfo[,1],"_",ctrlinfo[,2],sep="")
caseidx=which(data[,1]%in%casetmp)
ctrlidx=which(data[,1]%in%ctrltmp)
CASEOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".cases",sep="")
CTRLOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".ctrls",sep="")
write.table(data[caseidx,],CASEOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
write.table(data[ctrlidx,],CTRLOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
}


#GWAS2
fam=read.table("/lustre/scratch115/projects/ibdgwas/pre_imputation/qc/GWAS2/wtccc2_hg19.fam")
ctrlidx=which(fam[,6]==1)
caseidx=which(fam[,6]==2)
HLADIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/GWAS2/"
PREFIX="wtccc2_hg19.final.xMHC.HLA_"
SUFFIX=".4D.bestguess.txt"
LOCI=c("A","C","B","DRB1","DQA1","DQB1","DPB1")
for(loc in 1:length(LOCI)){
LOCUS=LOCI[loc]
data=read.table(paste(HLADIR,PREFIX,LOCUS,SUFFIX,sep=""),head=T,stringsAsFactors=F)
#casetmp=paste(caseinfo[,1],"_",caseinfo[,2],sep="")
#ctrltmp=paste(ctrlinfo[,1],"_",ctrlinfo[,2],sep="")
CASEOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".cases",sep="")
CTRLOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".ctrls",sep="")
write.table(data[caseidx,],CASEOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
write.table(data[ctrlidx,],CTRLOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
}


#GWAS3
fam=read.table("/lustre/scratch115/projects/ibdgwas/pre_imputation/qc/GWAS3/coreex_gaibdc_usgwas_qc.fam")
ctrlidx=which(fam[,6]==1)
caseidx=which(fam[,6]==2)
HLADIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/GWAS3/"
PREFIX="gwas3_final.unique.xMHC.HLA_"
SUFFIX=".4D.bestguess.txt"
LOCI=c("A","C","B","DRB1","DQA1","DQB1","DPB1")
for(loc in 1:length(LOCI)){
LOCUS=LOCI[loc]
data=read.table(paste(HLADIR,PREFIX,LOCUS,SUFFIX,sep=""),head=T,stringsAsFactors=F)
#casetmp=paste(caseinfo[,1],"_",caseinfo[,2],sep="")
#ctrltmp=paste(ctrlinfo[,1],"_",ctrlinfo[,2],sep="")
CASEOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".cases",sep="")
CTRLOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".ctrls",sep="")
write.table(data[caseidx,],CASEOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
write.table(data[ctrlidx,],CTRLOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
}

#new wave
fam=read.table("/lustre/scratch115/projects/ibdgwas/new_wave/new_wave_qc3b.final.unique.xMHC.fam")
ctrlidx=which(fam[,6]==1)
caseidx=which(fam[,6]==2)
HLADIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/new_wave/"
PREFIX="new_wave_qc3b.final.unique.xMHC.HLA_"
SUFFIX=".4D.bestguess.txt"
LOCI=c("A","C","B","DRB1","DQA1","DQB1","DPB1")
for(loc in 1:length(LOCI)){
LOCUS=LOCI[loc]
data=read.table(paste(HLADIR,PREFIX,LOCUS,SUFFIX,sep=""),head=T,stringsAsFactors=F)
#casetmp=paste(caseinfo[,1],"_",caseinfo[,2],sep="")
#ctrltmp=paste(ctrlinfo[,1],"_",ctrlinfo[,2],sep="")
CASEOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".cases",sep="")
CTRLOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".ctrls",sep="")
write.table(data[caseidx,],CASEOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
write.table(data[ctrlidx,],CTRLOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
}

#ichip

#new wave
fam=read.table("/lustre/scratch115/projects/ibdgwas/ichip/b37loukas/ichip_b37.final.xMHC.fam")
ctrlidx=which(fam[,6]==1)
caseidx=which(fam[,6]==2)
HLADIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/ichip/"
PREFIX="ichip_b37.final.xMHC.HLA_"
SUFFIX=".4D.bestguess.txt"
LOCI=c("A","C","B","DRB1","DQA1","DQB1","DPB1")
for(loc in 1:length(LOCI)){
LOCUS=LOCI[loc]
data=read.table(paste(HLADIR,PREFIX,LOCUS,SUFFIX,sep=""),head=T,stringsAsFactors=F)
#casetmp=paste(caseinfo[,1],"_",caseinfo[,2],sep="")
#ctrltmp=paste(ctrlinfo[,1],"_",ctrlinfo[,2],sep="")
CASEOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".cases",sep="")
CTRLOUTFILE=paste(HLADIR,PREFIX,LOCUS,SUFFIX,".ctrls",sep="")
write.table(data[caseidx,],CASEOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
write.table(data[ctrlidx,],CTRLOUTFILE,quote=F,col.names=T,row.names=F,sep="\t")
}


#check DRB1*01:03 / DRB*01:01
ALLELE="01:03"
LOCUS="DRB1"
DATASET="ichip"
DIR=NULL
THRES=0.7
if(DATASET=="GWAS1"){DIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/GWAS1/wtccc1_hg19.final.xMHC.HLA_"}
if(DATASET=="GWAS2"){DIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/GWAS2/wtccc2_hg19.final.xMHC.HLA_"}
if(DATASET=="GWAS3"){DIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/GWAS3/gwas3_final.unique.xMHC.HLA_"}
if(DATASET=="new_wave"){DIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/new_wave/new_wave_qc3b.final.unique.xMHC.HLA_"}
if(DATASET=="ichip"){DIR="/lustre/scratch115/projects/ibdgwas/HLA/HLA_imputations/4D/best_guess/ichip/ichip_b37.final.xMHC.HLA_"}

ctrls=read.table(paste(DIR,LOCUS,".4D.bestguess.txt.ctrls",sep=""),head=T)
cases=read.table(paste(DIR,LOCUS,".4D.bestguess.txt.cases",sep=""),head=T)

#checking posteriors for DRB1*01:03 carriers
summary(cases[c(which(cases[,2]==ALLELE),which(cases[,3]==ALLELE)),4])
summary(ctrls[c(which(ctrls[,2]==ALLELE),which(ctrls[,3]==ALLELE)),4])

#checking posteriors for all carriers
summary(cases[,4])
summary(ctrls[,4])


ctrls=ctrls[ctrls[,4]>=THRES,] #applying the threshold
ctrls0103=matrix(ncol=1,nrow=dim(ctrls)[1],0)
ctrls0103[which(ctrls[,2]==ALLELE),]=1
ctrls0103[which(ctrls[,3]==ALLELE),]=ctrls0103[which(ctrls[,3]==ALLELE),]+1

cases=cases[cases[,4]>=THRES,] #applying the threshold
cases0103=matrix(ncol=1,nrow=dim(cases)[1],0)
cases0103[which(cases[,2]==ALLELE),]=1
cases0103[which(cases[,3]==ALLELE),]=cases0103[which(cases[,3]==ALLELE),]+1

TABCASE=table(cases0103)
TABCTRL=table(ctrls0103)
if(length(TABCASE)==2){TABCASE=c(TABCASE,0)}
if(length(TABCTRL)==2){TABCTRL=c(TABCTRL,0)}
(TABCTRL[2]+2*TABCTRL[3])/(2*sum(TABCTRL))
(TABCASE[2]+2*TABCASE[3])/(2*sum(TABCASE))

#END