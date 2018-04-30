##

setwd("C:/Academic/SANGER/PSC2020/RNASEQ_RUN23806/")
r22891=read.table("run22891.sorted.deduplicated_percentage.txt")
r21121=read.table("run21121.sorted.deduplicated_percentage.txt")
r22294=read.table("run22294.sorted.deduplicated_percentage.txt")
r21364=read.table("run21364.sorted.deduplicated_percentage.txt")
r22219=read.table("run22219.sorted.deduplicated_percentage.txt")
r22226=read.table("run22226.sorted.deduplicated_percentage.txt")

#and now Liz's samples:
r23806=read.table("run23806_dedup_percentage.txt",stringsAsFactors=F)

r22891_500=r22891[which(r22891[,3]%in%c(30,31,32)|r22891[,3]%in%c(61,62,63)),]
r22891_200=r22891[which(r22891[,3]%in%c(33,34,35)|r22891[,3]%in%c(64,65,66)),]
library(wesanderson)
COLS=wes_palette(10, name = "Zissou", type = "continuous")
boxplot(r23806[,2],col=wes_palette(5, name = "Chevalier", type = "continuous")[1],at=0.75,xlim=c(0,9),ylim=c(0,90), ylab=list("% deduplicated",cex=1.1),xaxt='n')
boxplot(r22891_500[,4],col=COLS[2],add=T,at=2)
boxplot(r22891_200[,4],col=COLS[3],add=T,at=3)
boxplot(r21364[,4],col=COLS[4],add=T,at=4)
boxplot(r21121[,4],col=COLS[5],add=T,at=5)
boxplot(r22219[,4],col=COLS[9],add=T,at=6)
boxplot(r22226[,4],col=COLS[9],add=T,at=7)
boxplot(r22294[,4],col=COLS[10],add=T,at=8)
abline(v=1.5,lty=2,lwd=2,col="orange")
axis(1,at=c(0.75,2:8),labels=c("r23806\n?ng\nLibraryPrep?","r22891\n500ng\nKapa","r22891\n200ng\nKapa","r21364\n170-500ng\nKapa","r21121\n34-410ng\nTruSeq","r22219\n100ng\nKapa/NEB","r22226\n100ng\nTruSeq","r22294\n10ng\nSmartSeq2"),tick=F)

fragment_ratio=read.table("23806_compatible_fragment_ratio.txt",head=F,sep="\t")
boxplot(fragment_ratio[,2],col=wes_palette(5, name = "Chevalier", type = "continuous")[1],at=0.75,ylim=c(0.5,1), ylab=list("compatible fragment ratio",cex=1.1),xaxt='n')

assigned_fragments=read.table("23806_ISR.num_assigned_fragments.txt",head=F,sep="\t")
boxplot(assigned_fragments[,2]/1e+06,col=wes_palette(5, name = "Chevalier", type = "continuous")[1],at=0.75, ylab=list("number of assigned fragments by Salmon (in millions)",cex=1.1),xaxt='n')

percent_mapped=read.table("23806_percent_mapped.txt",head=F,sep="\t")
boxplot(percent_mapped[,2],col=wes_palette(5, name = "Chevalier", type = "continuous")[1],at=0.75, ylim=c(50,100),ylab=list("% mapped to transcriptome",cex=1.1),main=list("Percentage of reads mapped to transcriptome\n(num_mapped / num_procesed)",cex=1.1),xaxt='n')

#check these by samples
samples_list=read.table("samples_with_corresponding_Cram_files.txt",sep="\t",head=T,stringsAsFactors=F)
#make sure sample names are in the same format:
tmp1=matrix(unlist(strsplit(as.character(assigned_fragments[,1]),'_')),ncol=4,byrow=T)
sample_names1=paste(tmp1[,1],"_",tmp1[,2],".",tmp1[,3],sep="")
assigned_fragments_23806=cbind(sample_names1,assigned_fragments[,2])
colnames(assigned_fragments_23806)=c("sample_names","assigned_fragments")
tmp=matrix(unlist(strsplit(r23806[,1],'\\.')),ncol=3,byrow=T)
sample_names=paste(tmp[,1],tmp[,2],sep=".")
r23806_dedup=cbind(sample_names,r23806[,2])
colnames(r23806_dedup)=c("sample_names","dedup_perc")

mrg=merge(samples_list,assigned_fragments_23806,by.x="cramfile",by.y="sample_names",keep.all="x")
mrg_all=merge(mrg,r23806_dedup,by.x="cramfile",by.y="sample_names",keep.all="x")

#now plot for dedup percentage
plot(as.numeric(as.character(mrg_all$dedup_perc)),as.numeric(as.character(mrg_all$RIN)),xlab=list("% deduplicated",cex=1.1),ylab=list("RIN",cex=1.1),pch=19,col=wes_palette(5, name = "Chevalier", type = "continuous")[1])
plot(as.numeric(as.character(mrg_all$dedup_perc)),as.numeric(as.character(mrg_all$RNAconc)),xlab=list("% deduplicated",cex=1.1),ylab=list("RNA concentration (pg/ul)",cex=1.1),pch=19,col=COLS[10])
plot(as.numeric(as.character(mrg_all$dedup_perc)),as.numeric(as.character(mrg_all$cell_number)),xlab=list("% deduplicated",cex=1.1),ylab=list("Number of cells",cex=1.1),pch=19,col=COLS[2])


cor(as.numeric(as.character(mrg_all$dedup_perc)),as.numeric(as.character(mrg_all$cell_number)),use="complete.obs")
cor(as.numeric(as.character(mrg_all$dedup_perc)),as.numeric(as.character(mrg_all$RNAconc)),use="complete.obs")

#now plot for number of assigned fragments
plot(as.numeric(as.character(mrg_all$assigned_fragments))/1e+06,as.numeric(as.character(mrg_all$RIN)),xlab=list("assigned_fragments (Millions)",cex=1.1),ylab=list("RIN",cex=1.1),pch=19,col=wes_palette(5, name = "Chevalier", type = "continuous")[1])
plot(as.numeric(as.character(mrg_all$assigned_fragments))/1e+06,as.numeric(as.character(mrg_all$RNAconc)),xlab=list("assigned_fragments (Millions)",cex=1.1),ylab=list("RNA concentration (pg/ul)",cex=1.1),pch=19,col=COLS[10])
plot(as.numeric(as.character(mrg_all$assigned_fragments))/1e+06,as.numeric(as.character(mrg_all$cell_number)),xlab=list("assigned_fragments (Millions)",cex=1.1),ylab=list("Number of cells",cex=1.1),pch=19,col=COLS[2])

cor(as.numeric(as.character(mrg_all$assigned_fragments))/1e+06,as.numeric(as.character(mrg_all$cell_number)),use="complete.obs")
cor(as.numeric(as.character(mrg_all$assigned_fragments))/1e+06,as.numeric(as.character(mrg_all$RNAconc)),use="complete.obs")
cor(as.numeric(as.character(mrg_all$assigned_fragments))/1e+06,as.numeric(as.character(mrg_all$cell_number)),use="complete.obs")


## END
