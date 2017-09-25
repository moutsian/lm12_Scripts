#Sept 2017


data=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/22891_runs.txt",head=T,sep="\t")
mike_stats=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/Run_22891_NEB and Kapa Bake Off RNAseq QC_Full.txt",head=T,sep="\t")
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


totalseq_22891=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/fastQC_total_sequences_22891.txt")
totalseq_22219=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/fastQC_total_sequences_22219.txt")
totalseq_22226=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/fastQC_total_sequences_22226.txt")


dedup_22891=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/fastQC_deduplication_rate_22891.txt")
dedup_22219=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/fastQC_deduplication_rate_22219.txt")
dedup_22226=read.table("../PAXGENE/KIT COMPARISON APRIL 2017/fastQC_deduplication_rate_22226.txt")

boxplot(totalseq_22891[,1]*8.25,col="darkblue",at=1,xlim=c(0,4),ylim=c(0,max(c(totalseq_22219[,1],totalseq_22226[,1]))),ylab=list("Total Sequences"),main=list("Total Sequences, normalised*",cex=1.1))
boxplot(totalseq_22219[,1],col="darkred",add=T,at=2)
boxplot(totalseq_22226[,1],col="yellow",add=T,at=3)
axis(side=1,at=1:3,labels=c("22891*","22219","22226"))



boxplot(totalseq_22891[c(1:3,7:9,13:15,19:21),1]*8.25,col="darkblue",at=1,xlim=c(0,5),ylim=c(0,max(c(totalseq_22219[,1],totalseq_22226[,1]))),ylab=list("Total Sequences"),main=list("Total Sequences, normalised*",cex=1.1))
boxplot(totalseq_22891[c(4:6,10:12,16:18,22:24),1]*8.25,col="lightblue",add=T,at=2)
boxplot(totalseq_22219[,1],col="darkred",add=T,at=3)
boxplot(totalseq_22226[,1],col="yellow",add=T,at=4)
axis(side=1,at=1:4,labels=c("22891*\n500ng\n(Kapa&NEB)","22891*\n200ng\n(Kapa&NEB)","22219\n100ng\n(Kapa&NEB)","22226\n100ng\n(TruSeq)"),tick=F,pos=-4.5e+06)


boxplot(dedup_22891[c(1:3,7:9,13:15,19:21),1],col="darkblue",at=1,xlim=c(0,5),ylim=c(0,100),ylab=list("Total Sequences"),main=list("Deduplication Rate from fastQC",cex=1.1))
boxplot(dedup_22891[c(4:6,10:12,16:18,22:24),1],col="lightblue",add=T,at=2)
boxplot(dedup_22219[,1],col="darkred",add=T,at=3)
boxplot(dedup_22226[,1],col="yellow",add=T,at=4)
axis(side=1,at=1:4,labels=c("22891\n500ng\n(Kapa&NEB)","22891\n200ng\n(Kapa&NEB)","22219\n100ng\n(Kapa&NEB)","22226\n100ng\n(TruSeq)"),tick=F,pos=-8)




#END