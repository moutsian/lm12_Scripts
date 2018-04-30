## Make PCA plot for Liz's samples

setwd("C:/Academic/SANGER/PSC2020/RNASEQ/RNASEQ_RUN23806/")
rnaseq=read.table("23806_ISR.transcript_to_gene.counts",head=T)


#a) when you don't scale the rnaseq data:
pca <- prcomp(t(rnaseq))
pca.proportionvar <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100
#the above is alrady provided through summary:
pca_prop=summary(pca)
barplot(pca.proportionvar, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot (Gene counts, run 23806)", ylim=c(0,100))


#b) standardising the rnaseq data first:
scaled_rnaseq=scale(rnaseq)
pca_sc <- prcomp(t(scaled_rnaseq))
pca_sc_prop=summary(pca_sc)
pca_sc.proportionvar <- ((pca_sc$sdev^2) / (sum(pca_sc$sdev^2)))*100

barplot(pca_sc.proportionvar, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca_sc$sdev)), ylab="Proportion of variation (%)", main="Scree plot (Gene counts, run 23806)", ylim=c(0,100))

#in order to label the PCA points, I load the sample list and then re-order it so it matches the order of the samples in the PC data.
samples_list=read.table("samples_with_corresponding_Cram_files.txt",sep="\t",head=T,stringsAsFactors=F)
tmp1=matrix(unlist(strsplit(as.character(rownames(pca_sc$x)),'_')),ncol=7,byrow=T)
tmp1[,1]=gsub("X","",tmp1[,1])
samples_in_pca=paste(tmp1[,1],"_",tmp1[,2],".",tmp1[,3],sep="")
samples_list_ordered=samples_list[match(samples_in_pca, samples_list$cramfile),]
plot(pca_sc$x[,1],pca_sc$x[,2],type="n", main="PCA on gene counts from run23806", xlab=paste("PC1, ", round(pca_sc.proportionvar[1], 2), "%"), ylab=paste("PC2, ", round(pca_sc.proportionvar[2], 2), "%"))
points(pca_sc$x[,1],pca_sc$x[,2], col="black", pch=16, cex=1)
idx1=grep("UC", samples_list_ordered$individual)
idx2=grep("PSC", samples_list_ordered$individual)
points(pca_sc$x[idx1,1], pca_sc$x[idx1,2],col="black", pch=16, cex=1)
points(pca_sc$x[idx2,1], pca_sc$x[idx2,2],col="orange", pch=16, cex=1)
legend("bottom",c("UC","PSC"),pch=19,col=c("black","orange"))

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



##END