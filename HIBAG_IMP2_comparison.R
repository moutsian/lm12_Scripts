
setwd("C:/Academic/SANGER/HIBAG comparison of MS cohorts/")
AIMP2=read.table("A.ctrls.ms2012.Gil.SWE")
AHIBAG=read.table("A.ctrls.ms2012.HIBAG.SWE",head=T)
al1=gsub(":", "",AHIBAG[,2] )
al2=gsub(":", "",AHIBAG[,3] )
AHIBAG=cbind(AHIBAG,al1,al2)
Aconcordance=matrix(ncol=1,nrow=dim(AHIBAG)[1],0);
for(i in 1:dim(HIBAG)[1]){
hibag=c(as.numeric(as.character(AHIBAG[i,5])),as.numeric(as.character(AHIBAG[i,6])))
hibag=hibag[order(hibag)]
imp2=c(as.numeric(as.character(AIMP2[i,2])),as.numeric(as.character(AIMP2[i,3])))
imp2=imp2[order(imp2)]
if(hibag[1]==imp2[1] && hibag[2]==imp2[2]){
Aconcordance[i,1]=2;
}else if(hibag[1]==imp2[1] || hibag[1]==imp2[2] || hibag[2]==imp2[1] || hibag[2]==imp2[2]){
Aconcordance[i,1]=1
}
}

CIMP2=read.table("C.ctrls.ms2012.Gil.SWE")
CHIBAG=read.table("C.ctrls.ms2012.HIBAG.SWE",head=T)
al1=gsub(":", "",CHIBAG[,2] )
al2=gsub(":", "",CHIBAG[,3] )
CHIBAG=cbind(CHIBAG,al1,al2)
Cconcordance=matrix(ncol=1,nrow=dim(CHIBAG)[1],0);
for(i in 1:dim(HIBAG)[1]){
hibag=c(as.numeric(as.character(CHIBAG[i,5])),as.numeric(as.character(CHIBAG[i,6])))
hibag=hibag[order(hibag)]
imp2=c(as.numeric(as.character(CIMP2[i,2])),as.numeric(as.character(CIMP2[i,3])))
imp2=imp2[order(imp2)]
if(hibag[1]==imp2[1] && hibag[2]==imp2[2]){
Cconcordance[i,1]=2;
}else if(hibag[1]==imp2[1] || hibag[1]==imp2[2] || hibag[2]==imp2[1] || hibag[2]==imp2[2]){
Cconcordance[i,1]=1
}
}

BIMP2=read.table("B.ctrls.ms2012.Gil.SWE")
BHIBAG=read.table("B.ctrls.ms2012.HIBAG.SWE",head=T)
al1=gsub(":", "",BHIBAG[,2] )
al2=gsub(":", "",BHIBAG[,3] )
BHIBAG=cbind(BHIBAG,al1,al2)
Bconcordance=matrix(ncol=1,nrow=dim(BHIBAG)[1],0);
for(i in 1:dim(HIBAG)[1]){
hibag=c(as.numeric(as.character(BHIBAG[i,5])),as.numeric(as.character(BHIBAG[i,6])))
hibag=hibag[order(hibag)]
imp2=c(as.numeric(as.character(BIMP2[i,2])),as.numeric(as.character(BIMP2[i,3])))
imp2=imp2[order(imp2)]
if(hibag[1]==imp2[1] && hibag[2]==imp2[2]){
Bconcordance[i,1]=2;
}else if(hibag[1]==imp2[1] || hibag[1]==imp2[2] || hibag[2]==imp2[1] || hibag[2]==imp2[2]){
Bconcordance[i,1]=1
}
}

DRB1IMP2=read.table("DRB.ctrls.ms2012.Gil.SWE")
DRB1HIBAG=read.table("DRB1.ctrls.ms2012.HIBAG.SWE",head=T)
al1=gsub(":", "",DRB1HIBAG[,2] )
al2=gsub(":", "",DRB1HIBAG[,3] )
DRB1HIBAG=cbind(DRB1HIBAG,al1,al2)
DRB1concordance=matrix(ncol=1,nrow=dim(DRB1HIBAG)[1],0);
for(i in 1:dim(HIBAG)[1]){
hibag=c(as.numeric(as.character(DRB1HIBAG[i,5])),as.numeric(as.character(DRB1HIBAG[i,6])))
hibag=hibag[order(hibag)]
imp2=c(as.numeric(as.character(DRB1IMP2[i,2])),as.numeric(as.character(DRB1IMP2[i,3])))
imp2=imp2[order(imp2)]
if(hibag[1]==imp2[1] && hibag[2]==imp2[2]){
DRB1concordance[i,1]=2;
}else if(hibag[1]==imp2[1] || hibag[1]==imp2[2] || hibag[2]==imp2[1] || hibag[2]==imp2[2]){
DRB1concordance[i,1]=1
}
}

DQA1IMP2=read.table("DQA.ctrls.ms2012.Gil.SWE")
DQA1HIBAG=read.table("DQA1.ctrls.ms2012.HIBAG.SWE",head=T)
al1=gsub(":", "",DQA1HIBAG[,2] )
al2=gsub(":", "",DQA1HIBAG[,3] )
DQA1HIBAG=cbind(DQA1HIBAG,al1,al2)
DQA1concordance=matrix(ncol=1,nrow=dim(DQA1HIBAG)[1],0);
for(i in 1:dim(HIBAG)[1]){
hibag=c(as.numeric(as.character(DQA1HIBAG[i,5])),as.numeric(as.character(DQA1HIBAG[i,6])))
hibag=hibag[order(hibag)]
imp2=c(as.numeric(as.character(DQA1IMP2[i,2])),as.numeric(as.character(DQA1IMP2[i,3])))
imp2=imp2[order(imp2)]
if(hibag[1]==imp2[1] && hibag[2]==imp2[2]){
DQA1concordance[i,1]=2;
}else if(hibag[1]==imp2[1] || hibag[1]==imp2[2] || hibag[2]==imp2[1] || hibag[2]==imp2[2]){
DQA1concordance[i,1]=1
}
}


DQB1IMP2=read.table("DQB.ctrls.ms2012.Gil.SWE")
DQB1HIBAG=read.table("DQB1.ctrls.ms2012.HIBAG.SWE",head=T)
al1=gsub(":", "",DQB1HIBAG[,2] )
al2=gsub(":", "",DQB1HIBAG[,3] )
DQB1HIBAG=cbind(DQB1HIBAG,al1,al2)
DQB1concordance=matrix(ncol=1,nrow=dim(DQB1HIBAG)[1],0);
for(i in 1:dim(HIBAG)[1]){
hibag=c(as.numeric(as.character(DQB1HIBAG[i,5])),as.numeric(as.character(DQB1HIBAG[i,6])))
hibag=hibag[order(hibag)]
imp2=c(as.numeric(as.character(DQB1IMP2[i,2])),as.numeric(as.character(DQB1IMP2[i,3])))
imp2=imp2[order(imp2)]
if(hibag[1]==imp2[1] && hibag[2]==imp2[2]){
DQB1concordance[i,1]=2;
}else if(hibag[1]==imp2[1] || hibag[1]==imp2[2] || hibag[2]==imp2[1] || hibag[2]==imp2[2]){
DQB1concordance[i,1]=1
}
}

DPB1IMP2=read.table("DPB.ctrls.ms2012.Gil.SWE")
DPB1HIBAG=read.table("DPB1.ctrls.ms2012.HIBAG.SWE",head=T)
al1=gsub(":", "",DPB1HIBAG[,2] )
al2=gsub(":", "",DPB1HIBAG[,3] )
DPB1HIBAG=cbind(DPB1HIBAG,al1,al2)
DPB1concordance=matrix(ncol=1,nrow=dim(DPB1HIBAG)[1],0);
for(i in 1:dim(HIBAG)[1]){
hibag=c(as.numeric(as.character(DPB1HIBAG[i,5])),as.numeric(as.character(DPB1HIBAG[i,6])))
hibag=hibag[order(hibag)]
imp2=c(as.numeric(as.character(DPB1IMP2[i,2])),as.numeric(as.character(DPB1IMP2[i,3])))
imp2=imp2[order(imp2)]
if(hibag[1]==imp2[1] && hibag[2]==imp2[2]){
DPB1concordance[i,1]=2;
}else if(hibag[1]==imp2[1] || hibag[1]==imp2[2] || hibag[2]==imp2[1] || hibag[2]==imp2[2]){
DPB1concordance[i,1]=1
}
}

CONCO=matrix(nrow=7,ncol=8,0)
rownames(CONCO)=c("A","C","B","DRB1","DQA1","DQB1","DPB1")
CONCO[1,1:3]=table(Aconcordance)
CONCO[2,1:3]=table(Cconcordance)
CONCO[3,1:3]=table(Bconcordance)
CONCO[4,1:3]=table(DRB1concordance)
CONCO[5,1:3]=table(DQA1concordance)
CONCO[6,1:3]=c(0,table(DQB1concordance))
CONCO[7,1:3]=table(DPB1concordance)
Aidx=which(AHIBAG[,4]>0.9 & AIMP2[,5]>0.9)
Cidx=which(CHIBAG[,4]>0.9 & CIMP2[,5]>0.9)
Bidx=which(BHIBAG[,4]>0.9 & BIMP2[,5]>0.9)
DRB1idx=which(DRB1HIBAG[,4]>0.9 & DRB1IMP2[,5]>0.9)
DQA1idx=which(DQA1HIBAG[,4]>0.9 & DQA1IMP2[,5]>0.9)
DQB1idx=which(DQB1HIBAG[,4]>0.9 & DQB1IMP2[,5]>0.9)
DPB1idx=which(DRB1HIBAG[,4]>0.9 & DPB1IMP2[,5]>0.9)
CONCO[1,4:6]=c(0,table(Aconcordance[Aidx]))
CONCO[2,4:6]=c(0,table(Cconcordance[Cidx]))
CONCO[3,4:6]=c(0,0,table(Bconcordance[Bidx]))
CONCO[4,4:6]=c(0,0,table(DRB1concordance[DRB1idx]))
CONCO[5,4:6]=table(DQA1concordance[DQA1idx])
CONCO[6,4:6]=c(0,table(DQB1concordance[DQB1idx]))
CONCO[7,4:6]=c(0,table(DPB1concordance[DPB1idx]))
CONCO[1,7]=length(Aidx)
CONCO[2,7]=length(Cidx)
CONCO[3,7]=length(Bidx)
CONCO[4,7]=length(DRB1idx)
CONCO[5,7]=length(DQA1idx)
CONCO[6,7]=length(DQB1idx)
CONCO[7,7]=length(DPB1idx)
CONCO[1,8]=length(Aidx)/length(Aconcordance)
CONCO[2,8]=length(Cidx)/length(Cconcordance)
CONCO[3,8]=length(Bidx)/length(Bconcordance)
CONCO[4,8]=length(DRB1idx)/length(DRB1concordance)
CONCO[5,8]=length(DQA1idx)/length(DQA1concordance)
CONCO[6,8]=length(DQB1idx)/length(DQB1concordance)
CONCO[7,8]=length(DPB1idx)/length(DPB1concordance)


#plot concordance
plotconco=rbind(CONCO[,3]/length(Aconcordance),CONCO[,6]/CONCO[,7])
library(wesanderson)
SEQ=seq(0,1,0.05)
barplot(plotconco,ylim=c(0,1),col=wes_palette("Chevalier",2),width=4,beside=TRUE,ylab=list("Concordance",cex=1.1),panel.first={
segments(x0=rep(-0.04,20),x1=rep(5,20),y0=SEQ,y1=SEQ,col="gray",lty=2)
})

#END