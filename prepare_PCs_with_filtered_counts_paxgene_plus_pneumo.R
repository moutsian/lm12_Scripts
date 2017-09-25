#First bit in unix R to obtain PCs, second in Windows for plotting
# This is a nice site for statistics and visualisation using PCs in R: http://www.sthda.com/english/wiki/principal-component-analysis-in-r-prcomp-vs-princomp-r-software-and-data-mining


run21121=read.table("/lustre/scratch115/projects/paxgene/SalmonOutput/IU.21121.transcript_to_gene.counts",sep="\t",head=T)
run21364=read.table("/lustre/scratch115/projects/paxgene/SalmonOutput/IU.21364.transcript_to_gene.counts",sep="\t",head=T)
nasal=read.table("/lustre/scratch115/projects/pneumo-inf/txData.nasal.countsMatrix",head=T)
blood=read.table("/lustre/scratch115/projects/pneumo-inf/txData.blood.countsMatrix",head=T)

#load gene lists
mitogenes=read.table("/lustre/scratch115/projects/paxgene/gene_lists/mitochondrial_genes.txt",skip=1)
ribogenes=read.table("/lustre/scratch115/projects/paxgene/gene_lists/ribogenes.txt")
refseqgenes=read.table("/lustre/scratch115/projects/paxgene/gene_lists/gene_list_with_REFSEQ_protein_IDst.txt",skip=1)

###############################################################################################################
#first, run PCs after keeping only refseq genes and removing mito and ribo genes #
###############################################################################################################
overlap1=which(rownames(run21121)%in%rownames(nasal) & rownames(run21121)%in%rownames(blood) )
overlap2=which(rownames(nasal)%in%rownames(run21121) )

fulldata=cbind(run21121[overlap1,],run21364[overlap1,],blood[overlap2,],nasal[overlap2,])
tokeep=which(rownames(fulldata)%in%refseqgenes[,1] & (!rownames(fulldata)%in%ribogenes[,1]) & (!rownames(fulldata)%in%mitogenes[,1]) ) 
logfulldata=log2(fulldata[tokeep,]+1)
#now calculate principal components
# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
logfulldata.pca <- prcomp(logfulldata,
                 center = TRUE,
                 scale. = TRUE) 
				 
png("pcaplot.counts.log2.scaled.paxgene_with_pneumo.png")
plot(logfulldata.pca, type='l')
dev.off()				 


head(unclass(logfulldata.pca$rotation)[,1:4])
outfile="/lustre/scratch115/projects/paxgene/SalmonOutput/PAXGENE.refseq.logcounts.10pcs.paxgene_plus_pneumo.txt"
write.table(logfulldata.pca$rotation[,1:10],outfile,quote=F,row.names=T,col.names=T,sep="\t")

### The variance retained by each principal component can be obtained as follow :
# Eigenvalues
eig <- (logfulldata.pca$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.logfulldata.active <- data.frame(eig = eig, variance = variance,
                     cumvariance = cumvar)
head(eig.logfulldata.active)

png("pcaplot.counts.log2.scaled.var_by_pc.paxgene_with_pneumo.png")
barplot(eig.logfulldata.active[, 2], names.arg=1:nrow(eig.logfulldata.active), 
       main = "Variances",
       xlab = "Principal Components",
       ylab = "Percentage of variances",
       col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.logfulldata.active), 
      eig.logfulldata.active[, 2], 
      type="b", pch=19, col = "red")
dev.off()				 
	  
####################################################################################################################################################################################
####the below I am doing in Windows for plotting, reading data from /lustre/scratch115/projects/paxgene/SalmonOutput/PAXGENE.refseq.logcounts.10pcs.txt in the SANGER/PAXGENE dir. #
####################################################################################################################################################################################

setwd("C:/Academic/SANGER/PAXGENE/")
pcs=read.table("PAXGENE.refseq.logcounts.10pcs.txt",sep="\t")

plot(pcs[,1],pcs[,2],col=c(rep("darkblue",16),rep("orange",16)),main=list("PCA, Refseq Genes w/o mito and ribo genes, log2counts",cex=1.1),pch=c(rep(15:18,8),pt.cex=c(rep(c(1.5,1.5,1.2,1.1),8))))
legend(y=-0.1,x=-0.174,pch=19,col=c("darkblue","orange"),legend=c("Run21121","Run21364"))


addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visible
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.

  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

COL1=rgb(red=0,green=33,blue=71,max=255)
COL2=rgb(red=255,green=69,blue=0,max=255)
NEWCOL=addTrans(COL,80)
NEWCOL2=addTrans(COL2,80)
with(pcs, symbols(x=pcs[,1],y=pcs[,2],inches=1/5,circles=rep(c(3,3,2,1),8),bg=c(rep(NEWCOL,16),rep(NEWCOL2,16)),main=list("PCA, Refseq Genes w/o mito and ribo genes, log2counts",cex=1.1),ylab="PC2",xlab="PC1"))
legend(y=-0.1,x=-0.174,pch=19,col=c("darkblue","orange"),legend=c("Run21121","Run21364"))

#now for PC3 vs PC4
with(pcs, symbols(x=pcs[,1],y=pcs[,4],inches=1/5,circles=rep(c(3,3,2,1),8),bg=c(rep(NEWCOL,16),rep(NEWCOL2,16)),main=list("PCA, Refseq Genes w/o mito and ribo genes, log2counts",cex=1.1),ylab="PC4",xlab="PC3"))
legend(y=-0.1,x=-0.174,pch=19,col=c("darkblue","orange"),legend=c("Run21121","Run21364"))




#END
