#Oct 2013
#setwd("C:/Academic/MS2012/Rfiles")
args<-commandArgs(TRUE)
print(args)
from=as.numeric(args[1])
to=as.numeric(args[2])
COHORT=args[3]
INT_TERM=args[4]; #1 for 15:01, 2 for A02:01, ...,9 for DQB03:02 (as in the rank below with COND_LOCUS,COND_LOCUS2,...,COND_LOCUS9)
print(paste("Printing for SNPs: ",from," - ",to,".",sep=""))
DATADIR="/nfs/team152/loukas/from_oxford/MS2012/"
RESDIR="/nfs/team152/loukas/from_oxford/MS2012/RES/SEPT16/"
SCRIPTDIR="/nfs/users/nfs_l/lm12/Scripts/"
source(paste(SCRIPTDIR,"MS_full_model_sept16_SNP.R",sep=""))

#LOCUS="DPB"
#COHORT="UK"

COND_LOCUS="DRB"
COND_TYPE="1501"
COND_LOCUS2="A"
COND_TYPE2="0201"
COND_LOCUS3="DRB"
COND_TYPE3="1303"
COND_LOCUS4="DRB"
COND_TYPE4="0301"
COND_LOCUS5="DRB"
COND_TYPE5="0801"
COND_LOCUS6="B"
COND_TYPE6="3801"
COND_LOCUS7="B"
COND_TYPE7="4402"
COND_LOCUS8="B"
COND_TYPE8="5501"
COND_LOCUS9="DQB"
COND_TYPE9="0302"
COND_LOCUS10="DQA"
COND_TYPE10="0101"
COND_LOCUS11="DQB"
COND_TYPE11="0301"
COND_SNP="rs9277565"
COND_SNP2="rs2229092"
INTERACTOR="DRB1501"
if(INT_TERM==2){
INTERACTOR="A0201"
}else if(INT_TERM==3){
INTERACTOR="DRB1303"
}else if(INT_TERM==4){
INTERACTOR="DRB0301"
}else if(INT_TERM==5){
INTERACTOR="DRB0801"
}else if(INT_TERM==6){
INTERACTOR="B3801"
}else if(INT_TERM==7){
INTERACTOR="B4402"
}else if(INT_TERM==8){
INTERACTOR="B5501"
}else if(INT_TERM==9){
INTERACTOR="DQB0302"
}


casefile=paste("snps.cases.ms2012.",COHORT,sep="")
ctrlfile=paste("snps.ctrls.ms2012.",COHORT,sep="")
cond_snp_casefile=paste("snps.cases.ms2012.",COHORT,sep="")
cond_snp_ctrlfile=paste("snps.ctrls.ms2012.",COHORT,sep="")
cond_casefile=paste(COND_LOCUS,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile=paste(COND_LOCUS,".ctrls.ms2012.",COHORT,sep="")
cond_casefile2=paste(COND_LOCUS2,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile2=paste(COND_LOCUS2,".ctrls.ms2012.",COHORT,sep="")
cond_casefile3=paste(COND_LOCUS3,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile3=paste(COND_LOCUS3,".ctrls.ms2012.",COHORT,sep="")
cond_casefile4=paste(COND_LOCUS4,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile4=paste(COND_LOCUS4,".ctrls.ms2012.",COHORT,sep="")
cond_casefile5=paste(COND_LOCUS5,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile5=paste(COND_LOCUS5,".ctrls.ms2012.",COHORT,sep="")
cond_casefile6=paste(COND_LOCUS6,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile6=paste(COND_LOCUS6,".ctrls.ms2012.",COHORT,sep="")
cond_casefile7=paste(COND_LOCUS7,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile7=paste(COND_LOCUS7,".ctrls.ms2012.",COHORT,sep="")
cond_casefile8=paste(COND_LOCUS8,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile8=paste(COND_LOCUS8,".ctrls.ms2012.",COHORT,sep="")
cond_casefile9=paste(COND_LOCUS9,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile9=paste(COND_LOCUS9,".ctrls.ms2012.",COHORT,sep="")
cond_casefile10=paste(COND_LOCUS10,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile10=paste(COND_LOCUS10,".ctrls.ms2012.",COHORT,sep="")
cond_casefile11=paste(COND_LOCUS11,".cases.ms2012.",COHORT,sep="")
cond_ctrlfile11=paste(COND_LOCUS11,".ctrls.ms2012.",COHORT,sep="")

pc_cases=paste("PCA.cases.ms2012.",COHORT,sep="")
pc_ctrls=paste("PCA.ctrls.ms2012.",COHORT,sep="")
print(paste("full model analysis - to check SNPs in the DDR1 region for associations. w/o gender, w/o uncertainty, with last term being an interaction with ",INTERACTOR,", with 5 PCs, Printing from ",from," to ",to,", for SNPs in ",COHORT," cohort.",sep=""))

OUTPUTFILE=paste("full_model_sept16_SNP.",from,"_",to,".int_with.",INTERACTOR,".",COHORT,sep="")

logreg_14p_nounc_nogender_sept16_SNP_INT(casefile,ctrlfile,
cond_casefile,cond_ctrlfile,COND_LOCUS,COND_TYPE,"4D",
cond_casefile2,cond_ctrlfile2,COND_LOCUS2,COND_TYPE2,"4D",
cond_casefile3,cond_ctrlfile3,COND_LOCUS3,COND_TYPE3,"4D",
cond_casefile4,cond_ctrlfile4,COND_LOCUS4,COND_TYPE4,"4D",
cond_casefile5,cond_ctrlfile5,COND_LOCUS5,COND_TYPE5,"4D",
cond_casefile6,cond_ctrlfile6,COND_LOCUS6,COND_TYPE6,"4D",
cond_casefile7,cond_ctrlfile7,COND_LOCUS7,COND_TYPE7,"4D",
cond_casefile8,cond_ctrlfile8,COND_LOCUS8,COND_TYPE8,"4D",
cond_casefile9,cond_ctrlfile9,COND_LOCUS9,COND_TYPE9,"4D",
cond_casefile10,cond_ctrlfile10,COND_LOCUS10,COND_TYPE10,"4D",
cond_casefile11,cond_ctrlfile11,COND_LOCUS11,COND_TYPE11,"4D",
cond_snp_casefile,cond_snp_ctrlfile,COND_SNP,COND_SNP2,INT_TERM,
OUTPUTFILE,pc_cases,pc_ctrls,from,to,DATADIR,RESDIR,SCRIPTDIR);

#END
