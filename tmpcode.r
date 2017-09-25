

#javi's qc file:
/lustre/scratch113/teams/anderson/users/jga/001_Projects/Project3_Vaccines/pneumoPilot/README.QCfromFASTQ

d1=read.table("21121.1.postfilter.newidx.ISR.raw/quant.sf",head=T)
d2=read.table("21121.2.postfilter.newidx.ISR.raw/quant.sf",head=T)
d3=read.table("21121.3.postfilter.newidx.ISR.raw/quant.sf",head=T)
d4=read.table("21121.4.postfilter.newidx.ISR.raw/quant.sf",head=T)
d5=read.table("21121.5.postfilter.newidx.ISR.raw/quant.sf",head=T)
d6=read.table("21121.6.postfilter.newidx.ISR.raw/quant.sf",head=T)
d7=read.table("21121.7.postfilter.newidx.ISR.raw/quant.sf",head=T)
d8=read.table("21121.8.postfilter.newidx.ISR.raw/quant.sf",head=T)
d9=read.table("21121.9.postfilter.newidx.ISR.raw/quant.sf",head=T)
d10=read.table("21121.10.postfilter.newidx.ISR.raw/quant.sf",head=T)
d11=read.table("21121.11.postfilter.newidx.ISR.raw/quant.sf",head=T)
d12=read.table("21121.12.postfilter.newidx.ISR.raw/quant.sf",head=T)
d13=read.table("21121.13.postfilter.newidx.ISR.raw/quant.sf",head=T)
d14=read.table("21121.14.postfilter.newidx.ISR.raw/quant.sf",head=T)
d15=read.table("21121.15.postfilter.newidx.ISR.raw/quant.sf",head=T)
d16=read.table("21121.16.postfilter.newidx.ISR.raw/quant.sf",head=T)


 TPM=cbind(d1$TPM,d2$TPM,d3$TPM,d4$TPM,d5$TPM,d6$TPM,d7$TPM,d8$TPM,d9$TPM,d10$TPM,d11$TPM,d12$TPM,d13$TPM,d14$TPM,d15$TPM,d16$TPM) 
 NumReads=cbind(d1$NumReads,d2$NumReads,d3$NumReads,d4$NumReads,d5$NumReads,d6$NumReads,d7$NumReads,d8$NumReads,d9$NumReads,d10$NumReads,d11$NumReads,d12$NumReads,d13$NumReads,d14$NumReads,d15$NumReads,d16$NumReads) 

cor_pear_NumReads=matrix(ncol=16,nrow=16,0)
for(i in 1:16){
	for(j in 1:16){
	cor_pear_NumReads[i,j]=cor(NumReads[,i],NumReads[,j])
	}
}
