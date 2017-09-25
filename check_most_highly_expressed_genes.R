dataset="21364"
library="IU"
data=read.table(paste("/lustre/scratch115/projects/paxgene/SalmonOutput/",library,".",dataset,".transcript_to_gene.counts",sep=""),head=T)
#MAKE SURE THAT DATA COLUMNS ARE ORGANIZED BY SAMPLE (from 1 to 16) - OTHERWISE REORDER
orig2_5 =apply(data[,c(1,5,9,13)],1,sum)
mod2_5  =apply(data[,c(2,6,10,14)],1,sum)
mod1    =apply(data[,c(3,7,11,15)],1,sum)
mod0_5  =apply(data[,c(4,8,12,16)],1,sum)

write.table(  mod0_5[which(mod0_5>511000)],paste("~lm12/",library,dataset,"_0_5mod.toptranscripts.csv",sep=""),sep=",",col.names=F,row.names=T,quote=F)

#END