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

#Reading Alignments from bam ######################################################################
# Now read alignments in from the bam file:

BAMFILE="19743_5.1.bam.from.samtools" #note that this may have been moved to the QCtesting dir
aln = readGAlignments(BAMFILE) #this will take long to load and will look to take a lot of memory but it will actually be deallocated when done
# (at the peak though is about 12% of the total of gen1b). It contains all reads that were aligned to the genome

#If they are paired-end:
aln2 = readGAlignmentPairs(BAMFILE) #this is the one to use I think, as we have paired-end reads.


#check coverage (by chromosome)
cvg=coverage(aln2)

#Retrieve model (reference) from Ensembl/UCSC ######################################################
#This is required for quantification
#note that there may be ready libraries with the human genome preloaded, so it would definitely save time to use them. Have a look at p15 of : http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf 

txdb = makeTxDbFromBiomart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
library(TxDb.Hsapiens.UCSC.hg38.knownGene,lib="/software/team152/Rpackages/") #instead of making it ourselves, we are using the existing one - make sure it is the right one.

#library(TxDb.Hsapiens.BioMart.igis,lib="/software/team152/Rpackages/") # hopefully this is the equivalent for ensembl #*** this doesn't seem to work... also it is for b37


#For working with Ensembl, I got the advice from: https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/
# xxx_DB in the vignette is just a string to the SQLite db file path
##### ************** NOTE THAT The gene model must match the genome build the reads in the BAM file were aligned to **************************
ens84_txdb_pth <- './Homo_sapiens.GRCh38.84.sqlite'
if (!file.exists(ens84_txdb_pth)) {
    ens84_txdb_pth <- ensDbFromGtf(gtf="/software/team152/Rpackages/Homo_sapiens.GRCh38.84.gtf.original")
}
txdb_ens84 <- EnsDb(ens84_txdb_pth)
txdb_ens84  # Preview the metadata



#library(TxDb.Hsapiens.UCSC.hg19.knownGene,lib="/software/team152/Rpackages/")
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
#note that for the object above the IDs are Entrez Gene IDs
exbygene = exonsBy(txdb, "gene") #use.names=T cannot be used when grouping by gene
length(exbygene) #total gene count in exbygene
transcripts = exonsBy(txdb, by="tx", use.names=TRUE)
length(transcripts)

#exbygene is a GRangesList object. Here is some info on it : http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/GRangesList-class.html
exbygene_ensembl=exonsBy(txdb_ens84,by="gene")
length(exbygene_ensembl)
transcripts_ensembl= exonsBy(txdb_ens84, by="tx", use.names=TRUE)
transcripts_ensembl_bycds= cdsBy(txdb_ens84, by="tx", use.names=TRUE)
length(transcripts_ensembl)

# Quantification step (table of counts) ##############################################################

###UCSC####
#next, do quantification (read counts) with summarizeOverlaps --- The gene model must match the genome build the reads in the BAM file were aligned to **check with Javi
overlaps_INE= summarizeOverlaps(exbygene, aln2, mode="IntersectionNotEmpty",SingleEnd=FALSE) #the object it returns is a SummarizedExperiment object
head(table(assays(overlaps_INE)$counts))

#quantification *by transcript* using UCSC 
transcripts_overlap_INE= summarizeOverlaps(transcripts, aln2, mode="IntersectionNotEmpty",SingleEnd=FALSE) #the object it returns is a SummarizedExperiment object
head(table(assays(transcripts_overlap_INE)$counts))

####ENSEMBL ####
#Before continuing with the quantification I need to make sure that the seqlevels are the same before using summarizeOverlaps. Otherwise, the two objects will have
#no sequence levels in common and everything will be counted as zero.
#for instance check:
seqlevels(aln2)[1:10]
seqlevels(exbygene)[1:10]
seqlevels(exbygene_ensembl)[1:10]
#the latter has chromosomes coded as "1","2","3", etc..., whilst the others as "chr1","chr2","chr3", etc,....

#Changing seq levels of the ensembl objects: (make sure this corresponds to the original ensembl alignment, the rank is funny sometimes, e.g. chr21 comes after chr22 for some reason)
# seqlengths() can help in the matching process!

#for ensembl transcripts
newnames_transcripts_ensembl = seqlevels(transcripts_ensembl)
newnames_transcripts_ensembl[1:25]=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX","chr8","chr9","chr11","chr10","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chr19","chrY","chr22","chr21","chrM")
updated_transcripts_ensembl=renameSeqlevels(transcripts_ensembl, newnames_transcripts_ensembl)

#for ensembl genes
newnames_genes_ensembl=seqlevels(exbygene_ensembl)
newnames_genes_ensembl[1:24]=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX","chr8","chr9","chr11","chr10","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chr19","chr22","chr21","chrY")
newnames_genes_ensembl[48]="chrM"
updated_exbygene_ensembl=renameSeqlevels(exbygene_ensembl, newnames_genes_ensembl)

#same for the transcripts by cds
nn = seqlevels(transcripts_ensembl_bycds)
nn[1:25]=c("chr7","chr12","chr11","chr2","chr6","chr16","chr4","chr3","chr1","chr8","chr19","chr17","chr22","chr5","chr14","chrX","chr10","chr18","chr20","chr13","chr15","chrY","chr9","chr21","chrM")
updated_transcripts_ensembl_bycds=renameSeqlevels(transcripts_ensembl_bycds, nn)

#quantification by gene using Ensembl as reference
overlaps_ensembl_INE= summarizeOverlaps(updated_exbygene_ensembl, aln2, mode="IntersectionNotEmpty",SingleEnd=FALSE) #the object it returns is a SummarizedExperiment object
overlaps_ensembl_INE_ignorestrand= summarizeOverlaps(updated_exbygene_ensembl, aln2, mode="IntersectionNotEmpty",SingleEnd=FALSE,ignore.strand=TRUE) 

head(table(assays(overlaps_ensembl_INE)$counts))

#quantification by transcript using Ensembl transcript IDs to compare with Kallisto/Salmon from Javi
transcripts_overlap_ensembl_INE= summarizeOverlaps(updated_transcripts_ensembl, aln2, mode="IntersectionNotEmpty",SingleEnd=FALSE) #the object it returns is a SummarizedExperiment object
head(table(assays(transcripts_overlap_ensembl_INE)$counts))

#quantification by transcript using Ensembl transcript IDs to compare with Kallisto/Salmon from Javi --but with cdsBy function
transcripts_overlap_ensembl_bycds_INE= summarizeOverlaps(updated_transcripts_ensembl_bycds, aln2, mode="IntersectionNotEmpty",SingleEnd=FALSE) #the object it returns is a SummarizedExperiment object
head(table(assays(transcripts_overlap_ensembl_bycds_INE)$counts))


##now also get gene and transcripts lengths from ensembl
 gene_lengths=lengthOf(txdb_ens84,of="gene")
 transcript_lengths=lengthOf(txdb_ens84,of="tx")


# get output to compare with Javi's - transcripts by gene
#make an edgeR object
raw_counts_genes_ensembl_INE=cbind(names(overlaps_ensembl_INE),gene_lengths,assays(overlaps_ensembl_INE)$counts,TPM,RPKM)
gene1=DGEList(counts=as.numeric(raw_counts_genes_ensembl_INE[,3]),genes=raw_counts_genes_ensembl_INE[,1])
TPM=cpm(gene1,normalised.lib.sizes=FALSE,log=FALSE)
RPKM=rpkm(gene1,normalised.lib.sizes=FALSE,log=FALSE,gene.length=as.numeric(raw_counts_genes_ensembl_INE[,2]))
colnames(raw_counts_genes_ensembl_INE)=c("GENE_ID","gene_length","Raw_counts","TPM","RPKM")
write.table(raw_counts_genes_ensembl_INE,paste("/lustre/scratch115/projects/pneumo-inf/quantification_output_from_aligned_bams/",BAMFILE,".gene_quantification",sep=""),quote=F,row.names=F,col.names=T,sep="\t")


# get output to compare with Javi's - transcripts by gene
#make an edgeR object
raw_counts_genes_ensembl_INE_ignorestrand=cbind(names(overlaps_ensembl_INE_ignorestrand),gene_lengths,assays(overlaps_ensembl_INE_ignorestrand)$counts,TPM,RPKM)
gene1is=DGEList(counts=as.numeric(raw_counts_genes_ensembl_INE_ignorestrand[,3]),genes=raw_counts_genes_ensembl_INE_ignorestrand[,1])
TPM=cpm(gene1is,normalised.lib.sizes=FALSE,log=FALSE)
RPKM=rpkm(gene1is,normalised.lib.sizes=FALSE,log=FALSE,gene.length=as.numeric(raw_counts_genes_ensembl_INE_ignorestrand[,2]))
colnames(raw_counts_genes_ensembl_INE_ignorestrand)=c("GENE_ID","gene_length","Raw_counts","TPM","RPKM")
write.table(raw_counts_genes_ensembl_INE_ignorestrand,paste("/lustre/scratch115/projects/pneumo-inf/quantification_output_from_aligned_bams/",BAMFILE,".gene_quantification_ignoring_strand",sep=""),quote=F,row.names=F,col.names=T,sep="\t")



# get output to compare with Javi's
#make an edgeR object
raw_counts_transcripts_ensembl_INE=cbind(names(updated_transcripts_ensembl),transcript_lengths,assays(transcripts_overlap_ensembl_INE)$counts,TPM,RPKM)
trans1=DGEList(counts=as.numeric(raw_counts_transcripts_ensembl_INE[,3]),genes=raw_counts_transcripts_ensembl_INE[,1])
TPM=cpm(trans1,normalised.lib.sizes=FALSE,log=FALSE)
RPKM=rpkm(trans1,normalised.lib.sizes=FALSE,log=FALSE,gene.length=as.numeric(raw_counts_transcripts_ensembl_INE[,2]))
colnames(raw_counts_transcripts_ensembl_INE)=c("Transcript_ID","tx_length","Raw_counts","TPM","RPKM")
write.table(raw_counts_transcripts_ensembl_INE,paste("/lustre/scratch115/projects/pneumo-inf/quantification_output_from_aligned_bams/",BAMFILE,".tx_quantification",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

raw_counts_transcripts_bycds_ensembl_INE=cbind(names(updated_transcripts_ensembl_bycds),transcript_lengths,assays(transcripts_overlap_ensembl_bycds_INE)$counts,TPM,RPKM)
trans2=DGEList(counts=as.numeric(raw_counts_transcripts_bycds_ensembl_INE[,3]),genes=raw_counts_transcripts_bycds_ensembl_INE[,1])
TPM=cpm(trans2,normalised.lib.sizes=FALSE,log=FALSE)
RPKM=rpkm(trans2,normalised.lib.sizes=FALSE,log=FALSE,gene.length=as.numeric(raw_counts_transcripts_bycds_ensembl_INE[,2]))
colnames(raw_counts_transcripts_bycds_ensembl_INE)=c("Transcript_ID","tx_length","Raw_counts","TPM","RPKM")
write.table(raw_counts_transcripts_bycds_ensembl_INE,paste("/lustre/scratch115/projects/pneumo-inf/quantification_output_from_aligned_bams/",BAMFILE,".cdsBy.tx_quantification",sep=""),quote=F,row.names=F,col.names=T,sep="\t")


#output gene coordinates for Javi
tmp=genes(txdb_ens84,columns=c("gene_id","gene_seq_start","gene_seq_end"))
tmp_tbl=cbind(rownames(as.matrix(ranges(tmp))),as.character(seqnames(tmp)),start(ranges(tmp)),end(ranges(tmp)),width(ranges(tmp)))
colnames(tmp_tbl)=c("GeneID","Chrom","Start","End","Width")
write.table(tmp_tbl,"/lustre/scratch115/projects/pneumo-inf/ensembl_genes_with_coordinates_and_width_v2.txt",sep="\t",quote=F,row.names=F)

#now, if you have more than one samples, you can use edgeR for differential expression checks:
edger <- DGEList(assays(overlaps_INE)$counts, group=rownames(colData(overlaps_INE)))



#get info for specific gene, introns and exons ######################################################
ERAP1 <- "51752" #note that as mentioned earlier we are using Gene IDs.
ERAP1_txs = transcriptsBy(txdb, by="gene")[[ERAP1]]
#ERAP1_txs has one range per transcript. 
ERAP1_tx_names <- mcols(ERAP1_txs)$tx_name
#extracting exon regions of the gene:
ERAP1_exbytx <- exonsBy(txdb, by="tx", use.names=TRUE)[ERAP1_tx_names]
elementNROWS(ERAP1_exbytx)
#now intronic regions:
ERAP1_inbytx <- intronsByTranscript(txdb, use.names=TRUE)[ERAP1_tx_names]
elementNROWS(ERAP1_inbytx)

ERAP1_exons_overlaps_INE=  summarizeOverlaps(ERAP1_exbytx,aln2,mode="IntersectionNotEmpty")














#END