run21121=read.table("/lustre/scratch115/realdata/mdt3/projects/paxgene/data_from_other_projects/txData.21121.rpkm.intersect",sep="\t",head=T)
run21364=read.table("/lustre/scratch115/realdata/mdt3/projects/paxgene/data_from_other_projects/txData.21364.rpkm.intersect",sep="\t",head=T)
pneumo=read.table("/lustre/scratch115/realdata/mdt3/projects/paxgene/data_from_other_projects/txData.Allvaccines.rpkm.modified.intersect",sep="\t",head=T)

#I will add one to everything because I will then use the log2 values for the PCs
alldata=cbind(run21121[,-1],run21364[,-1],pneumo[,-1])
logalldata=log(alldata+1)

#now calculate principal components
# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
logalldata.pca <- prcomp(logalldata,
                 center = TRUE,
                 scale. = TRUE) 
				 
				 
png("pcaplot.rpkm.log.scaled.png")
plot(logalldata.pca, type='l')
dev.off()				 

####the below I am doing in Windows for plotting, reading data from first10pcs.mine_and_javi.rpkm.txt in the PAXGENE dir.

#rank for 21121 is: 21121.4,21121.12,21121.7, 21121.3, 21121.16, 21121.11, 21121.8, 21121.15, 21121.6, 21121.14, 21121.13,21121.5, 21121.2, 21121.9, 21121.10, 21121.1
#rank for 213634 is: 21364.4,21364.2,21364.1,21364.10,21364.6, 21364.5,21364.8,21364.9,21364.3,21364.15,21364.13,21364.7,21364.14,21364.11,21364.16,21364.12

#based on the above ranks, reorder the two runs, change symbols based on individuals.
data_reordered=data[c(16,13,4,1,12,9,3,7,14,15,6,2,11,10,8,5,19,18,25,17,22,21,28,23,24,20,30,32,27,29,26,31,33:dim(data)[1]),]
plot(data_reordered[,3],data_reordered[,4],col=c(rep("darkblue",16),rep("lightblue",16),rep("darkred",29),rep("red",38)),main=list("PCA, Genes RPKM",cex=1.1),pch=c(rep(c(rep(15,4),rep(16,4),rep(17,4),rep(18,4)),2),rep(19,67)))


#now change symbols by blood volume
plot(data_reordered[,3],data_reordered[,4],col=c(rep("darkblue",16),rep("lightblue",16),rep("darkred",29),rep("red",38)),main=list("PCA, Genes RPKM",cex=1.1),pch=c(rep(15:18,8),rep(19,67)))


#END
