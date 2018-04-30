# Run it with  R-3.3.0

setwd("/lustre/scratch115/realdata/mdt3/projects/psc2020/inhouse_data/rnaseq_data/from_dafni/")

library(BiocGenerics,lib="/software/team152/Rpackages/")
library(Biobase,lib="/software/team152/Rpackages/")
library(S4Vectors,lib="/software/team152/Rpackages/")
library(IRanges,lib="/software/team152/Rpackages/")
library(GenomeInfoDb,lib="/software/team152/Rpackages/")
library(BiocGenerics,lib="/software/team152/Rpackages/")
library(GenomicRanges,lib="/software/team152/Rpackages/")
library(SummarizedExperiment,lib="/software/team152/Rpackages/")

dafni=readRDS("treg_salmon_gene_abundances.rds")

#Dafni suggested I focus on day0 samples. This will require filtering on dafni$process
day0=which(dafni$process=="day0")
rest=which(dafni$process=="rest")
PMA=which(dafni$process=="PMA")
#this is helpful wrt SummarizedExperiment object: https://kasperdanielhansen.github.io/genbioconductor/html/SummarizedExperiment.html
colData(dafni)
assays(dafni)
counts_day0=assay(dafni, "counts")[,day0] 
counts_rest=assay(dafni, "counts")[,rest] 
counts_PMA=assay(dafni, "counts")[,PMA] 

liz=read.table("/lustre/scratch115/projects/psc2020/inhouse_data/rnaseq_data/run23806/SalmonOutput/23806_ISR.transcript_to_gene.counts")

#now ensure you are working on the shared genes
shared=which(rownames(counts_rest)%in%rownames(liz))

mrg1=merge(liz,counts_rest,by="row.names")
mrg2=merge(mrg1,counts_day0,by.x="Row.names",by.y="row.names")
mrg=merge(mrg2,counts_PMA,by.x="Row.names",by.y="row.names")
write.table(mrg,"merged_data_from_dafni_and_23806.txt",sep="\t",quote=F,col.names=T,row.names=F)

# pca time ##############################################
#(I do the following in Windows)

setwd("C:/Academic/SANGER/PSC2020/RNASEQ/RNASEQ_RUN23806/")
rnaseq=read.table("merged_data_from_dafni_and_23806.txt",head=T)

rnaseq_pca=rnaseq[,-1] #no gene names
# standardising the rnaseq data first:
scaled_rnaseq=scale(rnaseq_pca)
pca_sc <- prcomp(t(scaled_rnaseq))
pca_sc_prop=summary(pca_sc)
pca_sc.proportionvar <- ((pca_sc$sdev^2) / (sum(pca_sc$sdev^2)))*100

barplot(pca_sc.proportionvar, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca_sc$sdev)), ylab="Proportion of variation (%)", main="Scree plot (Gene counts, run 23806)", ylim=c(0,100))

#in order to label the PCA points, I load the sample list and then re-order it so it matches the order of the samples in the PC data.

PC_A=1
PC_B=2
VAR_A=round(pca_sc.proportionvar[PC_A], 2)
VAR_B=round(pca_sc.proportionvar[PC_B], 2)
samples_list=read.table("samples_with_corresponding_Cram_files.txt",sep="\t",head=T,stringsAsFactors=F)
tmp1=matrix(unlist(strsplit(as.character(rownames(pca_sc$x)[1:24]),'_')),ncol=7,byrow=T)
tmp1[,1]=gsub("X","",tmp1[,1])
samples_in_pca=paste(tmp1[,1],"_",tmp1[,2],".",tmp1[,3],sep="")
samples_list_ordered=samples_list[match(samples_in_pca, samples_list$cramfile),]
fulllist_ind=c(samples_list_ordered$individual,rownames(pca_sc$x)[25:158])
plot(pca_sc$x[,PC_A],pca_sc$x[,PC_B],type="n", main="PCA on gene counts from Dafni and r23806",col="red",xlab=paste("PC",PC_A,", ", VAR_A, "%"), ylab=paste("PC",PC_B,", ", VAR_B, "%"))
#points(pca_sc$x[,PC_A],pca_sc$x[,PC_B], col="gray", pch=16, cex=1)
idx1=grep("UC", fulllist_ind)
idx2=grep("PSC", fulllist_ind)
idx_rst=25:46
idx_day0=47:146
idx_pma=147:158
points(pca_sc$x[idx1,PC_A], pca_sc$x[idx1,PC_B],col="black", pch=16, cex=1)
points(pca_sc$x[idx2,PC_A], pca_sc$x[idx2,PC_B],col="orange", pch=16, cex=1)
points(pca_sc$x[idx_rst,PC_A], pca_sc$x[idx_rst,PC_B],col="darkorchid", pch=16, cex=1)
points(pca_sc$x[idx_day0,PC_A], pca_sc$x[idx_day0,PC_B],col="darkslategray1", pch=16, cex=1)
points(pca_sc$x[idx_pma,PC_A], pca_sc$x[idx_pma,PC_B],col="gray", pch=16, cex=1)
idx3=grep("CD4+", samples_list_ordered$celltype,fixed=T)
idx3a=grep("CD25+", samples_list_ordered$celltype,fixed=T)
idx3b=grep("CD25-", samples_list_ordered$celltype,fixed=T)
idx4=grep("CD8", samples_list_ordered$celltype,fixed=T)
#points(pca_sc$x[idx4,PC_A], pca_sc$x[idx4,PC_B],col="red", pch=16, cex=1)
legend("top",c("UC","PSC","dafni_rest","dafni_day0","dafni_PMA"),pch=19,col=c("black","orange","darkorchid","darkslategray1","gray"))



plot(pca_sc$x[,1],pca_sc$x[,2],type="n", main="PCA on gene counts from run23806", xlab=paste("PC1, ", round(pca_sc.proportionvar[1], 2), "%"), ylab=paste("PC2, ", round(pca_sc.proportionvar[2], 2), "%"))
points(pca_sc$x[,1],pca_sc$x[,2], col="black", pch=16, cex=1)
idx3=grep("CD4+", samples_list_ordered$celltype,fixed=T)
idx3a=grep("CD25+", samples_list_ordered$celltype,fixed=T)
idx3b=grep("CD25-", samples_list_ordered$celltype,fixed=T)
idx4=grep("CD8", samples_list_ordered$celltype,fixed=T)
points(pca_sc$x[idx3,1], pca_sc$x[idx3,2],col="darkblue", pch=19, cex=1)
#points(pca_sc$x[idx3a,1], pca_sc$x[idx3a,2],col="cyan", pch=19, cex=1)
#points(pca_sc$x[idx3b,1], pca_sc$x[idx3b,2],col="gray", pch=19, cex=1)
points(pca_sc$x[idx4,1], pca_sc$x[idx4,2],col="gold", pch=19, cex=1)
legend("bottom",c("CD4+","CD8+"),pch=19,col=c("darkblue","gold"))

high_expr_day0=which(apply(scaled_rnaseq[,idx_day0],1,mean)>5)


#END