#this is prior to any filtering
#pca=read.table("TILI_1KG_forPCA.30.evec",head=F,stringsAsFactors=F) 
#fam=read.table("TILI_1KG_forPCA.corrected.fam",stringsAsFactors=F)

#this is after filtering (first done Aug30)
#pca=read.table("TILI_1KG_v2.30.evec",head=F,stringsAsFactors=F)
#fam=read.table("TILI_1KG_forPCA_v2.fam",stringsAsFactors=F)

#this is prior to any filtering, for the re-organised version
pca=read.table("TILI_1KG_qc1d.30.evec",head=F,stringsAsFactors=F) 
fam=read.table("TILI_1KG_forPCA.fam",stringsAsFactors=F)


pcafam=cbind(pca,fam[,c(1,2,6)])
idx1=grepl("HG",pcafam[,1])
idx2=grepl("NA",pcafam[,1])
 
 pca_1kg_idx=which(idx1==T | idx2==T)
 pca_cases_idx=which(pcafam[,dim(pca)[2]+3]==2)
 pca_ctrls_idx=which(pcafam[,dim(pca)[2]+3]==1 & idx1==F & idx2==F)

plot(pca[,2],pca[,3],pch=19,xlab="PC1",ylab="PC2",main=list("PCA 1 vs 2\nTILI with 1kg",cex=1.1),col="gray")
points(pca[pca_cases_idx,2],pca[pca_cases_idx,3],pch=19,col="darkred")
points(pca[pca_ctrls_idx,2],pca[pca_ctrls_idx,3],pch=19,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)



 pc1=1
 pc2=2
plot(pcafam[,pc1+1],pcafam[,pc2+1],pch=19,xlab=paste("PC",pc1,sep=""),ylab=paste("PC",pc2,sep=""),main=list(paste("PC",pc1," vs PC",pc2,sep=""),cex=1.1),col="gray")
points(pcafam[pca_cases_idx,pc1+1],pcafam[pca_cases_idx,pc2+1],pch=19,col="darkred")
points(pcafam[pca_ctrls_idx,pc1+1],pcafam[pca_ctrls_idx,pc2+1],pch=19,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)
abline(h=-0.15,col="darkgrey",lty=2)
abline(v=-0.15,col="darkgrey",lty=2)
abline(h=0.075,col="darkgrey",lty=2)
abline(v=0.075,col="darkgrey",lty=2)
abline(h=-0.075,col="darkgrey",lty=2)
abline(v=-0.075,col="darkgrey",lty=2)
#following will just colour the sanger samples in a different colour, to check pc5-pc6 plot.
#sanger=grepl("_",pcafilt[,12])
#points(pcafilt[sanger,pc1+1],pcafilt[sanger,pc2+1],pch=19,col="orange")


##the following has more colours for 1KG
info=read.table("20130606_sample_info.txt",head=T,sep="\t",stringsAsFactors=F)
hapinfo=read.table("hapmap3_popinfo.txt",head=T,stringsAsFactors=F)
eur_samples_idx=which(info$Population=="CEU"|info$Population=="TSI"|info$Population=="IBS")
eur_samples=info[eur_samples_idx,1]
eur_pca_idx=which(pcafam[,dim(pca)[2]+1]%in%eur_samples)
gbr_samples_idx=which(info$Population=="GBR")
gbr_samples=info[gbr_samples_idx,1]
gbr_pca_idx=which(pcafam[,dim(pca)[2]+1]%in%gbr_samples)

fin_samples_idx=which(info$Population=="FIN")
fin_samples=info[fin_samples_idx,1]
fin_pca_idx=which(pcafam[,dim(pca)[2]+1]%in%fin_samples)

sas_samples_idx=which(info$Population=="GIH"|info$Population=="PJL"|info$Population=="BEB"|info$Population=="STU"|info$Population=="ITU")
sas_samples=info[sas_samples_idx,1]
sas_pca_idx=which(pcafam[,dim(pca)[2]+1]%in%sas_samples)
eas_samples_idx=which(info$Population=="CHB"|info$Population=="JPT"|info$Population=="CHS"|info$Population=="CDX"|info$Population=="KHV")
eas_samples=info[eas_samples_idx,1]
eas_pca_idx=which(pcafam[,dim(pca)[2]+1]%in%eas_samples)
afr_samples_idx=which(info$Population==""|info$Population=="YRI"|info$Population=="LWK"|info$Population=="GWD"|info$Population=="MSL"|info$Population=="ESN"|info$Population=="ASW"|info$Population=="ACB")
afr_samples=info[afr_samples_idx,1]
afr_pca_idx=which(pcafam[,dim(pca)[2]+1]%in%afr_samples)
amr_samples_idx=which(info$Population=="CLM"|info$Population=="PUR"|info$Population=="MEX"|info$Population=="PEL")
amr_samples=info[amr_samples_idx,1]
amr_pca_idx=which(pcafam[,dim(pca)[2]+1]%in%amr_samples)
all_idx=unique(c(amr_pca_idx,afr_pca_idx,eas_pca_idx,sas_pca_idx,gbr_pca_idx,eur_pca_idx,pca_cases_idx,pca_ctrls_idx,fin_pca_idx))
tmp=1:dim(pcafam)[1]
idx_miss=which(!tmp%in%all_idx)
hapinfo[which(hapinfo[,2]%in%pcafam[idx_miss,13]),] #from this you can see that all ones with missing info are mexicans

 pc1=9
 pc2=10
plot(pcafam[,pc1+1],pcafam[,pc2+1],pch=19,xlab=paste("PC",pc1,sep=""),ylab=paste("PC",pc2,sep=""),main=list(paste("PC",pc1," vs PC",pc2,sep=""),cex=1.1),col="white")
abline(h=0.028,lty=2)
points(pcafam[amr_pca_idx,pc1+1],pcafam[amr_pca_idx,pc2+1],pch=19,col="pink")
points(pcafam[idx_miss,pc1+1],pcafam[idx_miss,pc2+1],pch=19,col="pink")

points(pcafam[sas_pca_idx,pc1+1],pcafam[sas_pca_idx,pc2+1],pch=19,col="gray")
points(pcafam[eas_pca_idx,pc1+1],pcafam[eas_pca_idx,pc2+1],pch=19,col="purple")
points(pcafam[afr_pca_idx,pc1+1],pcafam[afr_pca_idx,pc2+1],pch=19,col="green")

points(pcafam[gbr_pca_idx,pc1+1],pcafam[gbr_pca_idx,pc2+1],pch=19,col="yellow")
points(pcafam[eur_pca_idx,pc1+1],pcafam[eur_pca_idx,pc2+1],pch=19,col="yellow")
points(pcafam[fin_pca_idx,pc1+1],pcafam[fin_pca_idx,pc2+1],pch=19,col="yellow2")

points(pcafam[pca_ctrls_idx,pc1+1],pcafam[pca_ctrls_idx,pc2+1],pch=4,col="darkblue")
points(pcafam[pca_cases_idx,pc1+1],pcafam[pca_cases_idx,pc2+1],pch=4,col="darkred")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)
abline(h=-0.15,col="darkgrey",lty=2)
abline(v=-0.15,col="darkgrey",lty=2)
abline(h=0.075,col="darkgrey",lty=2)
abline(v=0.075,col="darkgrey",lty=2)
abline(h=-0.0075,col="darkgrey",lty=2)
abline(v=-0.075,col="darkgrey",lty=2)
#following will just colour the sanger samples in a different colour, to check pc5-pc6 plot.
sanger=grepl("_",pcafam[,34])
#points(pcafam[sanger,pc1+1],pcafam[sanger,pc2+1],pch=4,col="orange")


##### 30 AUG: REMOVAL PROCESS:

toremove_pc1_cases=which(pcafam[pca_cases_idx,1+1]<0.01)
toremove_pc1_ctrls=which(pcafam[pca_ctrls_idx,1+1]<0.01)
#lenient:
#toremove_pc2_cases=which(pcafam[pca_cases_idx,2+1]<0.0075)
#toremove_pc2_ctrls=which(pcafam[pca_ctrls_idx,2+1]<0.0075)
#stricter:
toremove_pc2_cases=which(pcafam[pca_cases_idx,2+1]<0.0125)
toremove_pc2_ctrls=which(pcafam[pca_ctrls_idx,2+1]<0.0125)

#toremove_pc3_cases=which(pcafam[pca_cases_idx,3+1]< -0.015)
toremove_pc3_ctrls=which(pcafam[pca_ctrls_idx,3+1]>0.001)
toremove_pc3_cases=which(pcafam[pca_cases_idx,3+1]>0.001)
toremove_pc4_ctrls=which(pcafam[pca_ctrls_idx,4+1]> -0.0035)
toremove_pc4_cases=which(pcafam[pca_cases_idx,4+1]> -0.0035)
toremove_pc6_ctrls=which(pcafam[pca_ctrls_idx,6+1]< -0.0075)
toremove_pc6_cases=which(pcafam[pca_cases_idx,6+1]< -0.0075)
toremove_pc7_ctrls=which(pcafam[pca_ctrls_idx,7+1]>0.005)
toremove_pc7_cases=which(pcafam[pca_cases_idx,7+1]>0.005)
toremove_pc8_ctrls=which(pcafam[pca_ctrls_idx,8+1]>0.015)
toremove_pc8_cases=which(pcafam[pca_cases_idx,8+1]>0.015)
#toremove_pc7_ctrls2=which(pcafam[pca_ctrls_idx,7+1]< -0.025)
#toremove_pc7_cases2=which(pcafam[pca_cases_idx,7+1]< -0.025)

#***Des tin PC9, stin opoia tha mporousa na efarmosw ena akoma threshold, kai me to opoio tha exana 36 samples (3 cases, 33 ctrls):
toremove_pc9_cases=which(pcafam[pca_cases_idx,9+1]< -0.002)
toremove_pc9_ctrls=which(pcafam[pca_ctrls_idx,9+1]< -0.002)

toremove_cases=unique(c(toremove_pc1_cases,toremove_pc2_cases,toremove_pc3_cases,toremove_pc4_cases,toremove_pc6_cases,toremove_pc7_cases,toremove_pc8_cases,toremove_pc9_cases))
toremove_ctrls=unique(c(toremove_pc1_ctrls,toremove_pc2_ctrls,toremove_pc3_ctrls,toremove_pc4_ctrls,toremove_pc6_ctrls,toremove_pc7_ctrls,toremove_pc8_ctrls,toremove_pc9_ctrls))

#these are the cases and controls to be removed:
pcafam[pca_cases_idx[toremove_cases],33:34]
pcafam[pca_ctrls_idx[toremove_ctrls],33:34]
write.table(pcafam[pca_ctrls_idx[toremove_ctrls],33:34],"ctrls_to_remove_due_to_PCA_qc1d.txt",quote=F,col.names=F,row.names=F)
write.table(pcafam[pca_cases_idx[toremove_cases],33:34],"cases_to_remove_due_to_PCA_qc1d.txt",quote=F,col.names=F,row.names=F)
##replot


 pc1=5
 pc2=9
plot(pcafam[,pc1+1],pcafam[,pc2+1],pch=19,xlab=paste("PC",pc1,sep=""),ylab=paste("PC",pc2,sep=""),main=list(paste("PC",pc1," vs PC",pc2,sep=""),cex=1.1),col="white")
points(pcafam[amr_pca_idx,pc1+1],pcafam[amr_pca_idx,pc2+1],pch=19,col="purple")
points(pcafam[idx_miss,pc1+1],pcafam[idx_miss,pc2+1],pch=19,col="purple")
points(pcafam[sas_pca_idx,pc1+1],pcafam[sas_pca_idx,pc2+1],pch=19,col="purple")
points(pcafam[eas_pca_idx,pc1+1],pcafam[eas_pca_idx,pc2+1],pch=19,col="purple")
points(pcafam[afr_pca_idx,pc1+1],pcafam[afr_pca_idx,pc2+1],pch=19,col="purple")
points(pcafam[gbr_pca_idx,pc1+1],pcafam[gbr_pca_idx,pc2+1],pch=19,col="yellow")
points(pcafam[eur_pca_idx,pc1+1],pcafam[eur_pca_idx,pc2+1],pch=19,col="yellow")
points(pcafam[fin_pca_idx,pc1+1],pcafam[fin_pca_idx,pc2+1],pch=19,col="yellow")
points(pcafam[pca_cases_idx[-toremove_cases],pc1+1],pcafam[pca_cases_idx[-toremove_cases],pc2+1],pch=4,col="darkred")
points(pcafam[pca_ctrls_idx[-toremove_ctrls],pc1+1],pcafam[pca_ctrls_idx[-toremove_ctrls],pc2+1],pch=4,col="darkblue")
abline(h=0.15,col="darkgrey",lty=2)
abline(v=0.15,col="darkgrey",lty=2)
abline(h=-0.15,col="darkgrey",lty=2)
abline(v=-0.15,col="darkgrey",lty=2)
abline(h=0.075,col="darkgrey",lty=2)
abline(v=0.075,col="darkgrey",lty=2)
abline(h=-0.075,col="darkgrey",lty=2)
abline(v=-0.075,col="darkgrey",lty=2)
#following will just colour the sanger samples in a different colour, to check pc5-pc6 plot.
#sanger=grepl("_",pcafilt[,12])
#points(pcafilt[sanger,pc1+1],pcafilt[sanger,pc2+1],pch=19,col="orange")


#END