
#total fragments (from Salmon)
total21121=read.table("Stats across all runs/run21121.sorted.num_processed_fragments.txt",stringsAsFactors=F)
total21364=read.table("Stats across all runs/run21364.sorted.num_processed_fragments.txt",stringsAsFactors=F)
total22219=read.table("Stats across all runs/run22219.sorted.num_processed_fragments.txt",stringsAsFactors=F)
total22226=read.table("Stats across all runs/run22226.sorted.num_processed_fragments.txt",stringsAsFactors=F)
total22294=read.table("Stats across all runs/run22294.sorted.num_processed_fragments.txt",stringsAsFactors=F)
total22891=read.table("Stats across all runs/run22891.sorted.num_processed_fragments.txt",stringsAsFactors=F)

#assigned fragments (from Salmon)
assigned21121=read.table("Stats across all runs/run21121.sorted.num_assigned_fragments.txt",stringsAsFactors=F)
assigned21364=read.table("Stats across all runs/run21364.sorted.num_assigned_fragments.txt",stringsAsFactors=F)
assigned22219=read.table("Stats across all runs/run22219.sorted.num_assigned_fragments.txt",stringsAsFactors=F)
assigned22226=read.table("Stats across all runs/run22226.sorted.num_assigned_fragments.txt",stringsAsFactors=F)
assigned22294=read.table("Stats across all runs/run22294.sorted.num_assigned_fragments.txt",stringsAsFactors=F)
assigned22891=read.table("Stats across all runs/run22891.sorted.num_assigned_fragments.txt",stringsAsFactors=F)

#percentage mapped (from Salmon)
perc_mapped21121=read.table("Stats across all runs/run21121.sorted.percent_mapped_fragments.txt",stringsAsFactors=F)
perc_mapped21364=read.table("Stats across all runs/run21364.sorted.percent_mapped_fragments.txt",stringsAsFactors=F)
perc_mapped22219=read.table("Stats across all runs/run22219.sorted.percent_mapped_fragments.txt",stringsAsFactors=F)
perc_mapped22226=read.table("Stats across all runs/run22226.sorted.percent_mapped_fragments.txt",stringsAsFactors=F)
perc_mapped22294=read.table("Stats across all runs/run22294.sorted.percent_mapped_fragments.txt",stringsAsFactors=F)
perc_mapped22891=read.table("Stats across all runs/run22891.sorted.percent_mapped_fragments.txt",stringsAsFactors=F)

#now prepare three boxplots, but without the spike ins for run 22219
tokeep22219=which(total22219[,2]<= 6)
tokeep22226=which(total22226[,2]==7 | (total22226[,2]==8 & total22226[,3]<=4)) 

# i ) Total 
MAX=max(c(total21121[,4],total21364[,4],total22219[tokeep22219,4],total22226[tokeep22226,4],total22294[,4],total22891[,4]))
boxplot(total21121[,4],col="pink",at=1,xlim=c(0,7),ylim=c(0,MAX),ylab=list("Total Fragments"),main=list("Total Fragments, as obtained from Salmon",cex=1.1))
boxplot(total21364[,4],col="purple",at=2,add=T)
boxplot(total22219[tokeep22219,4],col="orange",at=3,add=T)
boxplot(total22226[tokeep22226,4],col="darkblue",at=4,add=T)
boxplot(total22294[,4],col="khaki",at=5,add=T)
boxplot(total22891[,4],col="dark gray",at=6,add=T)
axis(1,at=c(1:6),labels=c("21121\n16spl","21364\n16spl","22219\n8spl","22226\n8spl","22294\n7spl","22891\n66spl"),tick=F)


# ii) Assigned
tokeep22219=which(perc_mapped22219[,2]<= 6)
tokeep22226=which(perc_mapped22226[,2]==7 | (perc_mapped22226[,2]==8 & perc_mapped22226[,3]<=4)) 
MAX=max(c(perc_mapped21121[,4],perc_mapped21364[,4],perc_mapped22219[tokeep22219,4],perc_mapped22226[tokeep22226,4],perc_mapped22294[,4],perc_mapped22891[,4]))
boxplot(perc_mapped21121[,4],col="pink",at=1,xlim=c(0,7),ylim=c(0,MAX),ylab=list("Assigned Fragments"),main=list("Assigned Fragments, as obtained from Salmon",cex=1.1))
boxplot(perc_mapped21364[,4],col="purple",at=2,add=T)
boxplot(perc_mapped22219[tokeep22219,4],col="orange",at=3,add=T)
boxplot(perc_mapped22226[tokeep22226,4],col="darkblue",at=4,add=T)
boxplot(perc_mapped22294[,4],col="khaki",at=5,add=T)
boxplot(perc_mapped22891[,4],col="dark gray",at=6,add=T)
axis(1,at=c(1:6),labels=c("21121\n16spl","21364\n16spl","22219\n8spl","22226\n8spl","22294\n7spl","22891\n66spl"),tick=F)

# iii) Percentage 
tokeep22219=which(perc_mapped22219[,2]<= 6)
tokeep22226=which(perc_mapped22226[,2]==7 | (perc_mapped22226[,2]==8 & perc_mapped22226[,3]<=4)) 
MAX=max(c(perc_mapped21121[,4],perc_mapped21364[,4],perc_mapped22219[tokeep22219,4],perc_mapped22226[tokeep22226,4],perc_mapped22294[,4],perc_mapped22891[,4]))
boxplot(perc_mapped21121[,4],col="pink",at=1,xlim=c(0,7),ylim=c(0,MAX),ylab=list("% Assigned Fragments"),main=list("% Assigned Fragments, as obtained from Salmon",cex=1.1))
boxplot(perc_mapped21364[,4],col="purple",at=2,add=T)
boxplot(perc_mapped22219[tokeep22219,4],col="orange",at=3,add=T)
boxplot(perc_mapped22226[tokeep22226,4],col="darkblue",at=4,add=T)
boxplot(perc_mapped22294[,4],col="khaki",at=5,add=T)
boxplot(perc_mapped22891[,4],col="dark gray",at=6,add=T)
axis(1,at=c(1:6),labels=c("21121\n16spl","21364\n16spl","22219\n8spl","22226\n8spl","22294\n7spl","22891\n66spl"),tick=F)


# iv) Deduplicated percentage from fastQC
perc_dedup21121=read.table("Stats across all runs/run21121.sorted.deduplicated_percentage.txt",stringsAsFactors=F)
perc_dedup21364=read.table("Stats across all runs/run21364.sorted.deduplicated_percentage.txt",stringsAsFactors=F)
perc_dedup22219=read.table("Stats across all runs/run22219.sorted.deduplicated_percentage.txt",stringsAsFactors=F)
perc_dedup22226=read.table("Stats across all runs/run22226.sorted.deduplicated_percentage.txt",stringsAsFactors=F)
perc_dedup22294=read.table("Stats across all runs/run22294.sorted.deduplicated_percentage.txt",stringsAsFactors=F)
perc_dedup22891=read.table("Stats across all runs/run22891.sorted.deduplicated_percentage.txt",stringsAsFactors=F)
#for these we don't need to filter out spike ins, because we've done so already



#END