# FOR R 3.3.0 #

#basefile="/lustre/scratch115/projects/ibdgwas/new_wave/new_wave_qc3b.final.unique.xMHC"  #new wave
#basefile="/lustre/scratch115/projects/ibdgwas/pre_imputation/qc/GWAS3/gwas3_final.unique.xMHC" #GWAS3
#basefile="/lustre/scratch115/projects/ibdgwas/pre_imputation/qc/GWAS1/wtccc1_hg19.final.xMHC" #GWAS1
basefile="/lustre/scratch115/projects/ibdgwas/pre_imputation/qc/GWAS2/wtccc2_hg19.final.xMHC" #GWAS2

##installation from the tar.gz file for 3.3.0. For example:
#install.packages("/software/team152/Rpackages/S4Vectors_0.6.6.tar.gz",repos=NULL,type="source",lib="/software/team152/Rpackages/")

##After installation, for the latest R version (3.3.0) --note that I installed the packages by downloading the tar.gz files and then installing from them.
library(HIBAG,lib="/software/team152/Rpackages/")
#model.list=get(load("/lustre/scratch115/projects/ibdgwas/HLA/imputation_panels/HumanCoreExome-European-HLA4-hg19.RData")) #I have downloaded the prefit parameter set from the HIBAG website.
#model.list=get(load("/lustre/scratch115/projects/ibdgwas/HLA/imputation_panels/Affy500K-European-HLA4-hg19.RData")) #this should be the one for GWAS1.
model.list=get(load("/lustre/scratch115/projects/ibdgwas/HLA/imputation_panels/AffySNPv6-European-HLA4-hg19.RData")) #this should be the one for GWAS2.

#note that I use assembly hg19 (check whether hg19 should be used instead?)
geno19=hlaBED2Geno(bed.fn=paste(basefile,".bed",sep=""), fam.fn=paste(basefile,".fam",sep=""), bim.fn=paste(basefile,".bim",sep=""),assembly="hg19")

# HLA imputation at HLA-A
hla.id <- "A"
model <- hlaModelFromObj(model.list[[hla.id]])
summary(model)
# HLA allele frequencies
cbind(frequency = model$hla.freq)


# SNPs in the model
head(model$snp.id)
head(model$snp.position)


# best-guess genotypes and all posterior probabilities
#pred.guess18 <- predict(model, geno18, type="response+prob")
pred.guess19_pos <- predict(model, geno19, type="response+prob",match.type="Position") #this is safe to use as in our list of SNPs(from the MS study) we still have variants that were not assigned rsIDs back then.
summary(pred.guess19_pos)
pred.guess19_pos$value
pred.guess19_pos$postprob

# HLA imputation at HLA-C
hla.id <- "C"
model <- hlaModelFromObj(model.list[[hla.id]])
pred.HLAC.guess19_pos <- predict(model, geno19, type="response+prob",match.type="Position") #this is safe to use as in our list of SNPs(from the MS study) we still have variants that were not assigned rsIDs back then.

# HLA imputation at HLA-B
hla.id <- "B"
model <- hlaModelFromObj(model.list[[hla.id]])
pred.HLAB.guess19_pos <- predict(model, geno19, type="response+prob",match.type="Position") #this is safe to use as in our list of SNPs(from the MS study) we still have variants that were not assigned rsIDs back then.


# HLA imputation at HLA-DRB1
hla.id <- "DRB1"
model <- hlaModelFromObj(model.list[[hla.id]])
pred.HLADRB1.guess19_pos <- predict(model, geno19, type="response+prob",match.type="Position") #this is safe to use as in our list of SNPs(from the MS study) we still have variants that were not assigned rsIDs back then.


# HLA imputation at HLA-DQA1
hla.id <- "DQA1"
model <- hlaModelFromObj(model.list[[hla.id]])
pred.HLADQA1.guess19_pos <- predict(model, geno19, type="response+prob",match.type="Position") #this is safe to use as in our list of SNPs(from the MS study) we still have variants that were not assigned rsIDs back then.


# HLA imputation at HLA-DQB1
hla.id <- "DQB1"
model <- hlaModelFromObj(model.list[[hla.id]])
pred.HLADQB1.guess19_pos <- predict(model, geno19, type="response+prob",match.type="Position") #this is safe to use as in our list of SNPs(from the MS study) we still have variants that were not assigned rsIDs back then.

# HLA imputation at HLA-DPB1
hla.id <- "DPB1"
model <- hlaModelFromObj(model.list[[hla.id]])
pred.HLADPB1.guess19_pos <- predict(model, geno19, type="response+prob",match.type="Position") #this is safe to use as in our list of SNPs(from the MS study) we still have variants that were not assigned rsIDs back then.


#now save - buest guess
write.table(pred.guess19_pos$value,paste(basefile,".HLA_A.4D.bestguess.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(pred.HLAC.guess19_pos$value,paste(basefile,".HLA_C.4D.bestguess.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(pred.HLAB.guess19_pos$value,paste(basefile,".HLA_B.4D.bestguess.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(pred.HLADRB1.guess19_pos$value,paste(basefile,".HLA_DRB1.4D.bestguess.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(pred.HLADQA1.guess19_pos$value,paste(basefile,".HLA_DQA1.4D.bestguess.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(pred.HLADQB1.guess19_pos$value,paste(basefile,".HLA_DQB1.4D.bestguess.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(pred.HLADPB1.guess19_pos$value,paste(basefile,".HLA_DPB1.4D.bestguess.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)

#now save - all posterior probabilities
write.table(pred.guess19_pos$postprob,paste(basefile,".HLA_A.4D.postprob.txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
write.table(pred.HLAC.guess19_pos$postprob,paste(basefile,".HLA_C.4D.postprob.txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
write.table(pred.HLAB.guess19_pos$postprob,paste(basefile,".HLA_B.4D.postprob.txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
write.table(pred.HLADRB1.guess19_pos$postprob,paste(basefile,".HLA_DRB1.4D.postprob.txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
write.table(pred.HLADQA1.guess19_pos$postprob,paste(basefile,".HLA_DQA1.4D.postprob.txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
write.table(pred.HLADQB1.guess19_pos$postprob,paste(basefile,".HLA_DQB1.4D.postprob.txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
write.table(pred.HLADPB1.guess19_pos$postprob,paste(basefile,".HLA_DPB1.4D.postprob.txt",sep=""),sep="\t",col.names=T,row.names=T,quote=F)


#END

#END