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


CHR=21

ens84_txdb_pth <- './Homo_sapiens.GRCh38.84.sqlite'
if (!file.exists(ens84_txdb_pth)) {
    ens84_txdb_pth <- ensDbFromGtf(gtf="/software/team152/Rpackages/Homo_sapiens.GRCh38.84.gtf.original")
}

files19828 <- list.files(path="/lustre/scratch115/projects/pneumo-inf/19828.bamFiles.New.Sorted/", pattern = "\\.bam.sorted$")
files19828=paste("/lustre/scratch115/projects/pneumo-inf/19828.bamFiles.New.Sorted/",files19828,sep="")
files19743 <- list.files(path="/lustre/scratch115/projects/pneumo-inf/19743.bamFiles.New.Sorted/", pattern = "\\.bam.sorted$")
files19743=paste("/lustre/scratch115/projects/pneumo-inf/19743.bamFiles.New.Sorted/",files19743,sep="")
files=c(files19743,files19828)
INFO=matrix(ncol=2,nrow=length(files))
INFO[,1]=files
for(i in 1:dim(INFO)[1]){
which <- GRanges(paste("chr",CHR,":9000001-300000000",sep=""))
param <- ScanBamParam(which=which)
aln2 = readGAlignmentPairs(files[i],param=param,index=files[i]) #index file should have the name file as the bam with added extension "bai"
INFO[i,2]=length(aln2)
print(aln2)

txdb_ens84 <- EnsDb(ens84_txdb_pth)
txdb_ens84  # Preview the metadata
exbygene_ensembl=exonsBy(txdb_ens84,by="gene",filter=SeqnameFilter(CHR))
length(exbygene_ensembl)
# When performing overlap operations the seqlevels and genome of the objects must match. Here were modify the VCF to match the TxDb.
genome(aln2)=unique(genome(exbygene_ensembl))
seqlevels(exbygene_ensembl)=paste0("chr",seqlevels(exbygene_ensembl))
#make sure we know have overlapping levels:
intersect(seqlevels(exbygene_ensembl),seqlevels(aln2))
#quantification by gene using Ensembl as reference
genequant= summarizeOverlaps(exbygene_ensembl, aln2, mode="IntersectionNotEmpty",SingleEnd=FALSE,ignore.strand=TRUE) 
head(table(assays(genequant)$counts))
gene_counts=assays(genequant)$counts
gene_lengths=lengthOf(txdb_ens84,of="gene",filter=SeqnameFilter(CHR))
#now prepare edgeR objects
output_file=cbind(names(genequant),gene_lengths,gene_counts)
geneobj=DGEList(counts=gene_counts,genes=output_file[,1])
TPM=cpm(geneobj,normalised.lib.sizes=FALSE,log=FALSE)
RPKM=rpkm(geneobj,normalised.lib.sizes=FALSE,log=FALSE,gene.length=as.numeric(output_file[,2]))
output_file=cbind(output_file,TPM,RPKM)
colnames(output_file)=c("GENE_ID","gene_length","Raw_counts","TPM","RPKM")
write.table(output_file,paste(files[i],".gene_quant.chr",CHR,".from9Mb.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
}

#END
