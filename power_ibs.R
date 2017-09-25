
#from medical records
ibsMED_1_1=read.table("Z:/ibs10_MED_1.1.res",head=F)
ibsMED_1_2=read.table("Z:/ibs10_MED_1.2.res",head=F)
ibsMED_1_3=read.table("Z:/ibs10_MED_1.3.res",head=F)

#from self-reported
ibsSELF_1_1=read.table("Z:/ibs10_SELF_1.1.res",head=F)
ibsSELF_1_2=read.table("Z:/ibs10_SELF_1.2.res",head=F)
ibsSELF_1_3=read.table("Z:/ibs10_SELF_1.3.res",head=F)

plot(ibsSELF_1_1[,1],ibsSELF_1_1[,2],xlim=c(0.,0.3),col="darkblue",xlab=list("Risk Allele Frequency",cex=1.1),ylab=list("Power",cex=1.1),
main=list("Power to detect associations in UKB individuals with IBS\nand >450K population controls",cex=1.1),pch=19)
points(ibsSELF_1_2[,1],ibsSELF_1_2[,2],col="darkgreen",pch=19)
points(ibsMED_1_1[,1],ibsMED_1_1[,2],col="darkblue",pch=0)
points(ibsMED_1_2[,1],ibsMED_1_2[,2],col="darkgreen",pch=0)
abline(h=seq(0,1,0.2),lty=2,col="gray")
abline(v=seq(0.05,0.25,0.05),lty=2,col="gray")
legend(x=0.18,y=0.4,title="Genotype Relative Risk",c("1.1, self-reported","1.2, self-reported","1.1, from med records","1.2, from med records"),pch=c(19,19,0,0),bg="white",col=c("darkblue","darkgreen","darkblue","darkgreen"),lwd=2)

#END