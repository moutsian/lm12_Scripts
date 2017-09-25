# FOR R 3.3.0 #

##installation from the tar.gz file for 3.3.0. For example:
#install.packages("/software/team152/Rpackages/S4Vectors_0.6.6.tar.gz",repos=NULL,type="source",lib="/software/team152/Rpackages/")

##After installation, for the latest R version (3.3.0) --note that I installed the packages by downloading the tar.gz files and then installing from them.



#Library Loading #################################################################################
library(BiocGenerics,lib="/software/team152/Rpackages/")
library(S4Vectors,lib="/software/team152/Rpackages/")
library(IRanges,lib="/software/team152/Rpackages/")
library(GenomeInfoDb,lib="/software/team152/Rpackages/")
library(GenomicRanges,lib="/software/team152/Rpackages/")
library(XVector,lib="/software/team152/Rpackages/")
library(Biostrings,lib="/software/team152/Rpackages/")
library(Rsamtools,lib="/software/team152/Rpackages/")
library(Biobase,lib="/software/team152/Rpackages/")
library(SummarizedExperiment,lib="/software/team152/Rpackages/")
library(GenomicAlignments,lib="/software/team152/Rpackages/")
library(intervals,lib="/software/team152/Rpackages/")
library(genomeIntervals,lib="/software/team152/Rpackages/")
library(RSQLite)
library(AnnotationDbi,lib="/software/team152/Rpackages/")
library(GenomicFeatures,lib="/software/team152/Rpackages/")
library(limma,lib="/software/team152/Rpackages/")
library(edgeR,lib="/software/team152/Rpackages/")
library(ensembldb,lib="/software/team152/Rpackages/")


coding_genes=read.table("/nfs/users/nfs_l/lm12/coding_genes_and_processed_transcripts_from_biomart_july2016.txt",sep="\t",head=T)
blood=read.table("/lustre/scratch115/projects/pneumo-inf/TranscriptCountFromfastqFiles.kallisto.NewCleanBam/txData.blood.counts",head=T)
nasal=read.table("/lustre/scratch115/projects/pneumo-inf/TranscriptCountFromfastqFiles.kallisto.NewCleanBam/txData.nasal.counts",head=T)
bloodheader=read.table("/lustre/scratch115/projects/pneumo-inf/listSeqFiles.Blood2",head=T,sep="\t")
nasalheader=read.table("/lustre/scratch115/projects/pneumo-inf/listSeqFiles.Nasal2",head=T,sep="\t")
allcounts=cbind(blood,nasal)
#get list of coding genes:
tmpallcounts= gsub("\\..*","",rownames(allcounts))
tmp=which(tmpallcounts%in%as.character(coding_genes[,1]))
codingcounts=allcounts[tmp,]
moreinfo=read.table("/lustre/scratch115/projects/pneumo-inf/19743.bamFiles.New.Sorted/RNAextractions.PneumoRNAseqPilot.txt",head=T,sep="\t")

#now ready to create DGEList object and plot MDS
group <- c(rep("Blood",38),rep("Nasal",29))
d <- DGEList(counts = allcounts, group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d,verbose=T)
d <- estimateTagwiseDisp(d)
MDS=plotMDS(d, method="bcv", col=as.numeric(d$samples$group),ndim=10)
allsamples=c(as.character(bloodheader[,1]),as.character(nasalheader[,1]))
output_matrix=MDS$cmdscale.out
outputmatrix=cbind(allsamples,output_matrix)
colnames(outputmatrix)=c("SampleID","MDS1","MDS2","MDS3","MDS4","MDS5","MDS6","MDS7","MDS8","MDS9","MDS10")
#now output the mds info
write.table(outputmatrix,"/lustre/scratch115/projects/pneumo-inf/19743.bamFiles.New.Sorted/MDS.Javi.clean.txt",row.names=F,col.names=T,sep="\t",quote=F)


#now ready to create DGEList object and plot MDS
group <- c(rep("Blood",38),rep("Nasal",29))
cod <- DGEList(counts = codingcounts, group=group)
cod <- calcNormFactors(cod)
cod <- estimateCommonDisp(cod,verbose=T)
cod <- estimateTagwiseDisp(cod)
MDScod=plotMDS(cod, method="bcv", col=as.numeric(cod$samples$group),ndim=10)
allsamples=c(as.character(bloodheader[,1]),as.character(nasalheader[,1]))
output_matrix=MDScod$cmdscale.out
outputmatrix=cbind(allsamples,output_matrix)
colnames(outputmatrix)=c("SampleID","MDS1","MDS2","MDS3","MDS4","MDS5","MDS6","MDS7","MDS8","MDS9","MDS10")
#now output the mds info
write.table(outputmatrix,"/lustre/scratch115/projects/pneumo-inf/19743.bamFiles.New.Sorted/MDS.Javi.clean.coding_genes_only.txt",row.names=F,col.names=T,sep="\t",quote=F)




#END