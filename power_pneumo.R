#power plot for pneumococcal carriage status - oct2015
setwd("Z:/Scripts/")
res1_2=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_1.2.res")
res1_1=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_1.1.res")
res1_3=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_1.3.res")
res1_5=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_1.5.res")
res1_1=rbind(c(0,0),res1_1)
res1_2=rbind(c(0,0),res1_2)
res1_3=rbind(c(0,0),res1_3)
res1_5=rbind(c(0,0),res1_5)



plot(res1_2[,1],res1_2[,2],type='l',lwd=2,xlab=list("Risk Allele Frequency", cex=1.3),xlim=c(0,0.2),ylab=list("Power to detect association",cex=1.3),
main=list("Pneumococcus carriage status for healthy Nepalese children",cex=1.3))
points(res1_3[,1],res1_3[,2],type='l',lwd=2,col="blue")
points(res1_5[,1],res1_5[,2],type='l',lwd=2,col="red")
abline(h=0.8,col="gray",lwd=2,lty=2)
legend(x=0.15,y=0.4,title="Genotype Relative Risk",c("1.2","1.3","1.5"),pch=19,col=c("black","blue","red"),lwd=2)


#now for a dominant serotype

dres1_3=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_1.3.dom.res")
dres1_5=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_1.5.dom.res")
dres1_7=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_1.7.dom.res")
dres2_0=read.table("Z:/PNEUMOCOCCUS/healthy_GRR_2.dom.res")
dres1_3=rbind(c(0,0),dres1_3)
dres1_5=rbind(c(0,0),dres1_5)
dres1_7=rbind(c(0,0),dres1_7)
dres2_0=rbind(c(0,0),dres2_0)

plot(dres1_5[,1],dres1_5[,2],type='l',lwd=2,col="red",xlab=list("Risk Allele Frequency", cex=1.3),xlim=c(0,0.5),ylim=c(0,1),ylab=list("Power to detect association",cex=1.3),
main=list("Specific serotype carriage status\nSerotype Prevalence of 10% ",cex=1.3))
points(dres1_7[,1],dres1_7[,2],type='l',lwd=2,col="blue")
points(dres2_0[,1],dres2_0[,2],type='l',lwd=2,col="orange")
abline(h=0.8,col="gray",lwd=2,lty=2)
legend(x=0.3,y=0.4,title="Genotype Relative Risk",c("1.5","1.7","2.0"),pch=19,col=c("red","blue","orange"),lwd=2)


#END