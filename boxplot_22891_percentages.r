#Sept 2017


data=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/22891_runs.txt",head=T,sep="\t")
mike_stats=read.table("../PAXGENE/KIT_COMPARISON_APRIL_2017/Run_22891_NEB and Kapa Bake Off RNAseq QC_Full.txt",head=T,sep="\t")
#mapping rate by kit
boxplot(data[c(64:69,76:81),11],col="darkblue",at=1,xlim=c(0,3),ylim=c(0,100),ylab=list("Mapping Rate (%)"),main=list("Mapping Rate, as obtained from Salmon",cex=1.1))
boxplot(data[c(70:75,82:87),11],col="darkred",add=T,at=2)
axis(side=1,at=1:2,labels=c("NEB","kapa"))

#mapping rate by RNA input
boxplot(data[c(64:69,76:81),11],col="darkblue",at=1,xlim=c(0,3),ylim=c(0,100),ylab=list("Mapping Rate (%)"),main=list("Mapping Rate, as obtained from Salmon",cex=1.1))
boxplot(data[c(70:75,82:87),11],col="darkred",add=T,at=2)
axis(side=1,at=1:2,labels=c("NEB","kapa"))

#mapping rate by RNA input
boxplot(mike_stats[c(1:6),4]*100,col="darkblue",at=1,xlim=c(0,3),ylim=c(0,100),ylab=list("Exonic Rate (%)"),main=list("Exonic Rate, as obtained from Mike's stats",cex=1.1))
boxplot(mike_stats[c(7:12),4]*100,col="darkred",add=T,at=2)
axis(side=1,at=1:2,labels=c("NEB","kapa"))


boxplot(data[c(67:69,73:75,79:81,85:87),11],col="pink",at=1,xlim=c(0,3),ylim=c(0,100),ylab=list("Mapping Rate (%)"),main=list("Mapping Rate, as obtained from Salmon, by RNA input",cex=1.1))
boxplot(data[c(64:66,70:73,76:78,82:84),11],col="purple",add=T,at=2)
axis(side=1,at=1:2,labels=c("200ng","500ng"))

#END