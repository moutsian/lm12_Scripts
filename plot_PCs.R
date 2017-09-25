setwd("C:/Academic/SANGER/LAIV/PC plots/")
pc=read.table("allEdgeRmds.Condition.txt",head=T)
plot(pc[,1],pc[,2])
info_cond=read.table("listSeqFiles.Condition",head=T)
info_tissue=read.table("listSeqFiles.Tissue",head=T)
moreinfo=read.table("../RNAextractions.PneumoRNAseqPilot.txt",head=T,sep="\t")

nasal=which(info_tissue[,2]=="Nasal")
blood=which(info_tissue[,2]=="Blood")
nasal_minus5=which(info_tissue[,2]=="Nasal" & info_cond[,2]=="Minus5")
nasal_plus2=which(info_tissue[,2]=="Nasal" & info_cond[,2]=="Plus2")
blood_plus2=which(info_tissue[,2]=="Blood" & info_cond[,2]=="Plus2")
blood_minus5=which(info_tissue[,2]=="Blood" & info_cond[,2]=="Minus5")

par(mfrow=c(2,4))
#1
plot(pc[,1],pc[,2],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC2",cex=1.2))
points(pc[nasal_minus5,1],pc[nasal_minus5,2],col="darkblue",pch=19)
points(pc[nasal_plus2,1],pc[nasal_plus2,2],col="cyan",pch=19)
points(pc[blood_minus5,1],pc[blood_minus5,2],col="darkred",pch=19)
points(pc[blood_plus2,1],pc[blood_plus2,2],col="red",pch=19)
#2
PCA=3
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC3",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)
#lines(x=c(pc[blood_plus2[1],PCA],pc[blood_plus2[2],PCA]), y=c(pc[blood_plus2[1],PCB],pc[blood_plus2[2],PCB]))
#lines(x=c(pc[blood_plus2[3],PCA],pc[blood_plus2[4],PCA]), y=c(pc[blood_plus2[3],PCB],pc[blood_plus2[4],PCB]))
#lines(x=c(pc[blood_plus2[5],PCA],pc[blood_plus2[6],PCA]), y=c(pc[blood_plus2[5],PCB],pc[blood_plus2[6],PCB]))
#lines(x=c(pc[blood_plus2[7],PCA],pc[blood_plus2[8],PCA]), y=c(pc[blood_plus2[7],PCB],pc[blood_plus2[8],PCB]))
#lines(x=c(pc[blood_plus2[9],PCA],pc[blood_plus2[10],PCA]), y=c(pc[blood_plus2[9],PCB],pc[blood_plus2[10],PCB]))
#lines(x=c(pc[blood_plus2[11],PCA],pc[blood_plus2[12],PCA]), y=c(pc[blood_plus2[11],PCB],pc[blood_plus2[12],PCB]))
#3
PCA=1
PCB=3
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC3",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)
#4
PCA=1
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)
#4b
PCA=2
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC2",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)

#5
PCA=5
PCB=6
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC5",cex=1.2), ylab=list("PC6",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)

#6
PCA=7
PCB=8
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC7",cex=1.2), ylab=list("PC8",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)

#7
PCA=9
PCB=10
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC9",cex=1.2), ylab=list("PC10",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)




par(mfrow=c(2,2))
#1
plot(pc[,1],pc[,2],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC2",cex=1.2))
points(pc[nasal_minus5,1],pc[nasal_minus5,2],col="darkblue",pch=19)
points(pc[nasal_plus2,1],pc[nasal_plus2,2],col="cyan",pch=19)
points(pc[blood_minus5,1],pc[blood_minus5,2],col="darkred",pch=19)
points(pc[blood_plus2,1],pc[blood_plus2,2],col="red",pch=19)
#2
PCA=3
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC3",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)

#3
PCA=1
PCB=3
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC3",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)
#4
PCA=1
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)


###plot by well


well1=which(moreinfo[,2]==1)
well2=which(moreinfo[,2]==2)
well3=which(moreinfo[,2]==3)
well4=which(moreinfo[,2]==4)
well5=which(moreinfo[,2]==5)
well6=which(moreinfo[,2]==6)
well7=which(moreinfo[,2]==7)
well8=which(moreinfo[,2]==8)
well9=which(moreinfo[,2]==9)

par(mfrow=c(2,2),oma=c(0,0,1,0))
#1
PCA=1
PCB=2
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC2",cex=1.2))
points(pc[well1,PCA],pc[well1,PCB],col="darkblue",pch=19)
points(pc[well2,PCA],pc[well2,PCB],col="cyan",pch=19)
points(pc[well3,PCA],pc[well3,PCB],col="darkred",pch=19)
points(pc[well4,PCA],pc[well4,PCB],col="red",pch=19)
points(pc[well5,PCA],pc[well5,PCB],col="purple",pch=19)
points(pc[well6,PCA],pc[well6,PCB],col="gold",pch=19)
points(pc[well7,PCA],pc[well7,PCB],col="gray",pch=19)
points(pc[well8,PCA],pc[well8,PCB],col="chocolate",pch=19)
points(pc[well9,PCA],pc[well9,PCB],col="pink",pch=19)

#2
PCA=3
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC3",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[well1,PCA],pc[well1,PCB],col="darkblue",pch=19)
points(pc[well2,PCA],pc[well2,PCB],col="cyan",pch=19)
points(pc[well3,PCA],pc[well3,PCB],col="darkred",pch=19)
points(pc[well4,PCA],pc[well4,PCB],col="red",pch=19)
points(pc[well5,PCA],pc[well5,PCB],col="purple",pch=19)
points(pc[well6,PCA],pc[well6,PCB],col="gold",pch=19)
points(pc[well7,PCA],pc[well7,PCB],col="gray",pch=19)
points(pc[well8,PCA],pc[well8,PCB],col="chocolate",pch=19)
points(pc[well9,PCA],pc[well9,PCB],col="pink",pch=19)

#3
PCA=1
PCB=3
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC3",cex=1.2))
points(pc[well1,PCA],pc[well1,PCB],col="darkblue",pch=19)
points(pc[well2,PCA],pc[well2,PCB],col="cyan",pch=19)
points(pc[well3,PCA],pc[well3,PCB],col="darkred",pch=19)
points(pc[well4,PCA],pc[well4,PCB],col="red",pch=19)
points(pc[well5,PCA],pc[well5,PCB],col="purple",pch=19)
points(pc[well6,PCA],pc[well6,PCB],col="gold",pch=19)
points(pc[well7,PCA],pc[well7,PCB],col="gray",pch=19)
points(pc[well8,PCA],pc[well8,PCB],col="chocolate",pch=19)
points(pc[well9,PCA],pc[well9,PCB],col="pink",pch=19)
#4
PCA=1
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[well1,PCA],pc[well1,PCB],col="darkblue",pch=19)
points(pc[well2,PCA],pc[well2,PCB],col="cyan",pch=19)
points(pc[well3,PCA],pc[well3,PCB],col="darkred",pch=19)
points(pc[well4,PCA],pc[well4,PCB],col="red",pch=19)
points(pc[well5,PCA],pc[well5,PCB],col="purple",pch=19)
points(pc[well6,PCA],pc[well6,PCB],col="gold",pch=19)
points(pc[well7,PCA],pc[well7,PCB],col="gray",pch=19)
points(pc[well8,PCA],pc[well8,PCB],col="chocolate",pch=19)
points(pc[well9,PCA],pc[well9,PCB],col="pink",pch=19)
title("PCs, samples coloured by well",outer=T)



###plot by sample
library(RColorBrewer)
tmp=table(moreinfo[,7])
allcols=c(brewer.pal(dim(tmp), "Spectral"),brewer.pal(dim(tmp),"Set3"))
allpchs=sample(15:20,dim(tmp),replace=TRUE)
par(mfrow=c(2,2),oma=c(0,0,1,0))
#1
PCA=1
PCB=2
plot(pc[,PCA],pc[,PCB],pch=19,col="white",xlab=list("PC1",cex=1.2), ylab=list("PC2",cex=1.2))
for(i in 1:dim(tmp)[1]){
idx=which(moreinfo[,7]==names(tmp)[i])
print(idx)
points(pc[idx,PCA],pc[idx,PCB],col=allcols[i],pch=allpchs[i])
}

#2
PCA=3
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="white",xlab=list("PC3",cex=1.2), ylab=list("PC4",cex=1.2))
for(i in 1:dim(tmp)[1]){
idx=which(moreinfo[,7]==names(tmp)[i])
print(idx)
points(pc[idx,PCA],pc[idx,PCB],col=allcols[i],pch=allpchs[i])
}

#3
PCA=1
PCB=3
plot(pc[,PCA],pc[,PCB],pch=19,col="white",xlab=list("PC1",cex=1.2), ylab=list("PC3",cex=1.2))
for(i in 1:dim(tmp)[1]){
idx=which(moreinfo[,7]==names(tmp)[i])
print(idx)
points(pc[idx,PCA],pc[idx,PCB],col=allcols[i],pch=allpchs[i])
}

#4
PCA=1
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="white",xlab=list("PC1",cex=1.2), ylab=list("PC4",cex=1.2))
for(i in 1:dim(tmp)[1]){
idx=which(moreinfo[,7]==names(tmp)[i])
print(idx)
points(pc[idx,PCA],pc[idx,PCB],col=allcols[i],pch=allpchs[i])
}
title("PCs, samples coloured by individual ID",outer=T)


par(mfrow=c(1,3))
PCA=3
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="white",xlab=list("PC3",cex=1.2), ylab=list("PC4",cex=1.2))
for(i in 1:dim(tmp)[1]){
idx=which(moreinfo[,7]==names(tmp)[i])
print(idx)
points(pc[idx,PCA],pc[idx,PCB],col=allcols[i],pch=allpchs[i],bg=BG,lwd=2,cex=1.2)
}
title("PCs, kallisto quantification, samples coloured by individual ID",outer=F)
#2
PCA=2
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC2",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[nasal_minus5,PCA],pc[nasal_minus5,PCB],col="darkblue",pch=19)
points(pc[nasal_plus2,PCA],pc[nasal_plus2,PCB],col="cyan",pch=19)
points(pc[blood_minus5,PCA],pc[blood_minus5,PCB],col="darkred",pch=19)
points(pc[blood_plus2,PCA],pc[blood_plus2,PCB],col="red",pch=19)
title("PCs, samples coloured by condition and tissue",outer=F)
#3
PCA=2
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC2",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[well1,PCA],pc[well1,PCB],col="darkblue",pch=19)
points(pc[well2,PCA],pc[well2,PCB],col="cyan",pch=19)
points(pc[well3,PCA],pc[well3,PCB],col="darkred",pch=19)
points(pc[well4,PCA],pc[well4,PCB],col="red",pch=19)
points(pc[well5,PCA],pc[well5,PCB],col="purple",pch=19)
points(pc[well6,PCA],pc[well6,PCB],col="gold",pch=19)
points(pc[well7,PCA],pc[well7,PCB],col="gray",pch=19)
points(pc[well8,PCA],pc[well8,PCB],col="chocolate",pch=19)
points(pc[well9,PCA],pc[well9,PCB],col="pink",pch=19)
title("PCs, samples coloured by well",outer=F)


###plot by carrier status
COL1="navy"
COL2="red"
COL3="orange"
natural_carrier=which(moreinfo$sample_ID2=="sample3" | moreinfo$sample_ID2 =="sample14")
exp_samples=c("sample9","sample11","sample12","sample13","sample17","sample18","sample21")
exp_carriers=which(as.character(moreinfo$sample_ID2)%in%exp_samples & as.character(moreinfo$Condition)=="Plus2")
exp_carriers_before=which(as.character(moreinfo$sample_ID2)%in%exp_samples & as.character(moreinfo$Condition)=="Minus5")
par(mfrow=c(2,2))
#1
PCA=1
PCB=2
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC2",cex=1.2))
points(pc[natural_carrier,PCA],pc[natural_carrier,PCB],col=COL1,pch=19)
points(pc[exp_carriers,1],pc[exp_carriers,PCB],col=COL2,pch=19)
points(pc[exp_carriers_before,1],pc[exp_carriers_before,2],col=COL3,pch=19)
#2
PCA=3
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC3",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[natural_carrier,PCA],pc[natural_carrier,PCB],col=COL1,pch=19)
points(pc[exp_carriers,PCA],pc[exp_carriers,PCB],col=COL2,pch=19)
points(pc[exp_carriers_before,PCA],pc[exp_carriers_before,PCB],col=COL3,pch=19)

#3
PCA=1
PCB=3
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC3",cex=1.2))
points(pc[natural_carrier,PCA],pc[natural_carrier,PCB],col=COL1,pch=19)
points(pc[exp_carriers,PCA],pc[exp_carriers,PCB],col=COL2,pch=19)
points(pc[exp_carriers_before,PCA],pc[exp_carriers_before,2],col=COL3,pch=19)
#4
PCA=1
PCB=4
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[natural_carrier,PCA],pc[natural_carrier,PCB],col=COL1,pch=19)
points(pc[exp_carriers,PCA],pc[exp_carriers,PCB],col=COL2,pch=19)
points(pc[exp_carriers_before,PCA],pc[exp_carriers_before,PCB],col=COL3,pch=19)
title("PCs, kallisto quantification, samples coloured by carrier status",outer=T)
legend(-1,1.3,col=c(COL1,COL2,COL3),c("natural carriers","exp.carriers on +2","exp.carriers on -5"),pch=19,cex=0.75,pt.cex=1,ncol=1,title.adj=0.5,y.intersp=0.65,inset=c(0,-0.2))





#just testing
par(mfrow=c(1,1))
PCA=1
PCB=2
plot(pc[,PCA],pc[,PCB],pch=19,col="gray",xlab=list("PC1",cex=1.2), ylab=list("PC4",cex=1.2))
points(pc[natural_carrier,PCA],pc[natural_carrier,PCB],col=COL1,pch=19)
points(pc[exp_carriers,PCA],pc[exp_carriers,PCB],col=COL2,pch=19)
points(pc[exp_carriers_before,PCA],pc[exp_carriers_before,PCB],col=COL3,pch=19)
title("PCs, kallisto quantification, samples coloured by carrier status",outer=T)
legend(-1,1.3,col=c(COL1,COL2,COL3),c("natural carriers","exp.carriers on +2","exp.carriers on -5"),pch=19,cex=0.75,pt.cex=1,ncol=1,title.adj=0.5,y.intersp=0.65,inset=c(0,-0.2))






#END
