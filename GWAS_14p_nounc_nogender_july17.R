#################################################################
# The function takes as input a file in the form of *_details.txt, that 
# Alex has produced, with posteriors for all HLA alleles. Note that 
# the first two columns will not be included (as in the original file they
# are sampleID and chromosome). Also note that,as it is, the function takes
# as in input a combined case/controls prob file, although these are separate in the 
# directory.
# It is assumed that the haplotypes are transmitted independently of each other, which
# makes the probabilities for genotypes being hom1,het,hom2  straightforward to calculate.
# Variables:
# status is the binary case/control status (0/1)
# probs is the matrix with haplotypes as rows, alleles as columns.
# allele is the column of the allele of interest.
# model1 is, as before:
# Model: 1=Additive, 2=Dominant, 3=Recessive, 4=General (not yet implemented) and 5=Heterozygote

SCRIPTDIR="/nfs/users/nfs_l/lm12/Scripts/"  
RELTOL=1e-04

#the following is without the SNP-SNP interaction we are testing (to be used to test interaction between DDX39B and IL7R) .
gwas14_nounc_noG_july17_INT = function(interaction_term,status,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,PC){ ## where INT refers to H8 here, which is included as an interaction with DRB*15:01 (presence/absence)
  pH1=H1;
  pH2=H2;
  pH3=H3;
  pH4=H4;
  pH5=H5;
  pH6=H6;
  pH7=H7;
  pH8=H8;
  pH9=H9;
  pH10=H10;
  pH11=H11;
  pH12=H12;
  pH13=H13;
  pH14=H14;

  pH1[pH1>1]=1
  pH2[pH2>1]=1 # in case we get probs > 1 because of rounding
  pH3[pH3>1]=1
  pH4[pH4>1]=1
  pH5[pH5>1]=1
  pH6[pH6>1]=1
  pH7[pH7>1]=1
  pH8[pH8>1]=1
  pH9[pH9>1]=1
  pH10[pH10>1]=1
  pH11[pH11>1]=1
  pH12[pH12>1]=1
  pH13[pH13>1]=1
  pH14[pH14>1]=1



G1=abs(apply(pH1,1,which.max)-3);
G2=abs(apply(pH2,1,which.max)-3);
G3=abs(apply(pH3,1,which.max)-3);
G4=abs(apply(pH4,1,which.max)-3);
G5=abs(apply(pH5,1,which.max)-3);
G6=abs(apply(pH6,1,which.max)-3);
G7=abs(apply(pH7,1,which.max)-3);
G8=abs(apply(pH8,1,which.max)-3);
G9=abs(apply(pH9,1,which.max)-3);
G10=abs(apply(pH10,1,which.max)-3);
G11=abs(apply(pH11,1,which.max)-3);
G12=abs(apply(pH12,1,which.max)-3);
G13=abs(apply(pH13,1,which.max)-3);
G14=abs(apply(pH14,1,which.max)-3);


G1HOM=G1
G1HOM[G1HOM==1]=0
G1HOM[G1HOM==2]=1

I1501=G1
I1501[I1501>1]=1

G2HOM=G2
G2HOM[G2HOM==1]=0
G2HOM[G2HOM==2]=1

G4HOM=G4
G4HOM[G4HOM==1]=0
G4HOM[G4HOM==2]=1

G9DOM=G9
G9DOM[G9DOM>1]=1 #dominant effect

I1=NULL
if(interaction_term==1){
I1=G1
}else if(interaction_term==2){
I1=G2
}else if(interaction_term==3){
I1=G3
}else if(interaction_term==4){
I1=G4HOM
}else if(interaction_term==5){
I1=G5
}else if(interaction_term==6){
I1=G6;
}else if(interaction_term==7){
I1=G7;
}else if(interaction_term==8){
I1=G8;
}else if(interaction_term==9){ #this is what was used for the paper 
I1=G9DOM;
}
I1[I1==2]=1
G10=G10*I1501
G11=G11*I1 # Note that this is the effect of the last HLA allele (not SNP!) included in the data, in the presence/absence of the interactor (DQB*03:02 in the paper) 
# G10 is the new HLA/ SNP to test (where we test all in turn)

# note that in the G1,G2,...,G5 genotype arrays the number of copies are that of the actual allele. So G[i]=2 means that the SNP or HLA file has two copies of the HLA allele
# or 2 copies of the minor SNP allele.

M14= glm(status~G1+G1HOM+G2+G2HOM+G3+G4+G4HOM+G5+G6+G7+G8+G9DOM+G10+G11+G12+G13+G14+PC,family=binomial) #independent model for DRB*15:01 (G1) and A*02:01 (G2), recessive for DRB*03:01 (4G), additive model for rest

if(dim(summary(M14)$coeff)[1]==23){ 
fit_results = data.frame(summary(M14)$coeff[1,1],
summary(M14)$coeff[2,1],summary(M14)$coeff[3,1],summary(M14)$coeff[4,1],
summary(M14)$coeff[5,1],summary(M14)$coeff[6,1],summary(M14)$coeff[7,1],
summary(M14)$coeff[8,1],summary(M14)$coeff[9,1],summary(M14)$coeff[10,1],
summary(M14)$coeff[11,1],summary(M14)$coeff[12,1],summary(M14)$coeff[13,1],
summary(M14)$coeff[14,1],summary(M14)$coeff[15,1],summary(M14)$coeff[16,1],
summary(M14)$coeff[17,1],summary(M14)$coeff[18,1],
summary(M14)$coeff[2,2],summary(M14)$coeff[3,2],summary(M14)$coeff[4,2],
summary(M14)$coeff[5,2],summary(M14)$coeff[6,2], summary(M14)$coeff[7,2],
summary(M14)$coeff[8,2],summary(M14)$coeff[9,2],summary(M14)$coeff[10,2],
summary(M14)$coeff[11,2],summary(M14)$coeff[12,2],summary(M14)$coeff[13,2],
summary(M14)$coeff[14,2],summary(M14)$coeff[15,2],summary(M14)$coeff[16,2],
summary(M14)$coeff[17,2],summary(M14)$coeff[18,2],
summary(M14)$coeff[2,4],summary(M14)$coeff[3,4],summary(M14)$coeff[4,4],
summary(M14)$coeff[5,4],summary(M14)$coeff[6,4], summary(M14)$coeff[7,4],
summary(M14)$coeff[8,4],summary(M14)$coeff[9,4],summary(M14)$coeff[10,4],
summary(M14)$coeff[11,4],summary(M14)$coeff[12,4],summary(M14)$coeff[13,4],
summary(M14)$coeff[14,4],summary(M14)$coeff[15,4],summary(M14)$coeff[16,4],
summary(M14)$coeff[17,4],summary(M14)$coeff[18,4],
M14$null, M14$dev)
#print(fit_results)
#fit_results=as.table(fit.results)
colnames(fit_results)=c("mu","gamma1","gamma1HOM","gamma2","gamma2HOM","gamma3","gamma4","gamma4HOM","gamma5","gamma6","gamma7","gamma8","gamma9","gamma10INT","gamma11INT","gamma12","gamma13NONHLASNP","gamma14",
"se1","se1HOM","se2","se2HOM","se3","se4","se4HOM","se5","se6","se7","se8","se9","se10","se11","se12","se143NONHLASNP","se14",
"pval1","pval1HOM","pval2","pval2HOM","pval3","pval4","pval4HOM","pvalue5","pvalue6","pvalue7","pvalue8","pvalue9","pvalue10INT","pvalue11INT","pvalue12SNP","pvalue13NONHLASNP","pvalue14",
"dev_NULL","dev_M")
}else{
fit_results=rep(NA,54)
}
 return (fit_results);

}



# the following is as above, but where instead of a single term for the last SNP, we have three terms based on the copies of the IL7R SNP.
#with an added term for the SNP-SNP interaction we are testing (to be used to test interaction between DDX39B and IL7R) .
gwas14_nounc_noG_july17_WITH_3INT = function(interaction_term,status,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,PC){ ## where INT refers to H8 here, which is included as an interaction with DRB*15:01 (presence/absence)
  pH1=H1;
  pH2=H2;
  pH3=H3;
  pH4=H4;
  pH5=H5;
  pH6=H6;
  pH7=H7;
  pH8=H8;
  pH9=H9;
  pH10=H10;
  pH11=H11;
  pH12=H12;
  pH13=H13;
  pH14=H14;

  pH1[pH1>1]=1
  pH2[pH2>1]=1 # in case we get probs > 1 because of rounding
  pH3[pH3>1]=1
  pH4[pH4>1]=1
  pH5[pH5>1]=1
  pH6[pH6>1]=1
  pH7[pH7>1]=1
  pH8[pH8>1]=1
  pH9[pH9>1]=1
  pH10[pH10>1]=1
  pH11[pH11>1]=1
  pH12[pH12>1]=1
  pH13[pH13>1]=1
  pH14[pH14>1]=1



G1=abs(apply(pH1,1,which.max)-3);
G2=abs(apply(pH2,1,which.max)-3);
G3=abs(apply(pH3,1,which.max)-3);
G4=abs(apply(pH4,1,which.max)-3);
G5=abs(apply(pH5,1,which.max)-3);
G6=abs(apply(pH6,1,which.max)-3);
G7=abs(apply(pH7,1,which.max)-3);
G8=abs(apply(pH8,1,which.max)-3);
G9=abs(apply(pH9,1,which.max)-3);
G10=abs(apply(pH10,1,which.max)-3);
G11=abs(apply(pH11,1,which.max)-3);
G12=abs(apply(pH12,1,which.max)-3);
G13=abs(apply(pH13,1,which.max)-3);
G14=abs(apply(pH14,1,which.max)-3);

G1HOM=G1
G1HOM[G1HOM==1]=0
G1HOM[G1HOM==2]=1

I1501=G1
I1501[I1501>1]=1

G2HOM=G2
G2HOM[G2HOM==1]=0
G2HOM[G2HOM==2]=1

G4HOM=G4
G4HOM[G4HOM==1]=0
G4HOM[G4HOM==2]=1

G9DOM=G9
G9DOM[G9DOM>1]=1 #dominant effect

I1=NULL
if(interaction_term==1){
I1=G1
}else if(interaction_term==2){
I1=G2
}else if(interaction_term==3){
I1=G3
}else if(interaction_term==4){
I1=G4HOM
}else if(interaction_term==5){
I1=G5
}else if(interaction_term==6){
I1=G6;
}else if(interaction_term==7){
I1=G7;
}else if(interaction_term==8){
I1=G8;
}else if(interaction_term==9){ #this is what was used for the paper 
I1=G9DOM;
}
I1[I1==2]=1
G10=G10*I1501
G11=G11*I1 # Note that this is the effect of the last HLA allele (not SNP!) included in the data, in the presence/absence of the interactor (DQB*03:02 in the paper) 
# G10 is the new HLA/ SNP to test (where we test all in turn)


#the added bit here: the three interaction terms for the DDX39B effect based on the IL7R genotype

#0 copies
IL7R_0=G13
IL7R_0[IL7R_0==0]=-1
IL7R_0[IL7R_0>0]=0
IL7R_0[IL7R_0==-1]=1

#1 copy
IL7R_1=G13
IL7R_1[IL7R_1==2]=0

#2 copies
IL7R_2=G13
IL7R_2[IL7R_2==1]=0
IL7R_2[IL7R_2==2]=1

#now prepare the three interaction terms
G14INT0=G14*IL7R_0
G14INT1=G14*IL7R_1
G14INT2=G14*IL7R_2

# note that in the G1,G2,...,G5 genotype arrays the number of copies are that of the actual allele. So G[i]=2 means that the SNP or HLA file has two copies of the HLA allele
# or 2 copies of the minor SNP allele.

M14= glm(status~G1+G1HOM+G2+G2HOM+G3+G4+G4HOM+G5+G6+G7+G8+G9DOM+G10+G11+G12+G13+G14INT0+G14INT1+G14INT2+PC,family=binomial) #independent model for DRB*15:01 (G1) and A*02:01 (G2), recessive for DRB*03:01 (4G), additive model for rest

if(dim(summary(M14)$coeff)[1]==25){ 
fit_results = data.frame(summary(M14)$coeff[1,1],
summary(M14)$coeff[2,1],summary(M14)$coeff[3,1],summary(M14)$coeff[4,1],
summary(M14)$coeff[5,1],summary(M14)$coeff[6,1],summary(M14)$coeff[7,1],
summary(M14)$coeff[8,1],summary(M14)$coeff[9,1],summary(M14)$coeff[10,1],
summary(M14)$coeff[11,1],summary(M14)$coeff[12,1],summary(M14)$coeff[13,1],
summary(M14)$coeff[14,1],summary(M14)$coeff[15,1],summary(M14)$coeff[16,1],
summary(M14)$coeff[17,1],summary(M14)$coeff[18,1],summary(M14)$coeff[19,1],summary(M14)$coeff[20,1],
summary(M14)$coeff[2,2],summary(M14)$coeff[3,2],summary(M14)$coeff[4,2],
summary(M14)$coeff[5,2],summary(M14)$coeff[6,2], summary(M14)$coeff[7,2],
summary(M14)$coeff[8,2],summary(M14)$coeff[9,2],summary(M14)$coeff[10,2],
summary(M14)$coeff[11,2],summary(M14)$coeff[12,2],summary(M14)$coeff[13,2],
summary(M14)$coeff[14,2],summary(M14)$coeff[15,2],summary(M14)$coeff[16,2],
summary(M14)$coeff[17,2],summary(M14)$coeff[18,2],summary(M14)$coeff[19,2],summary(M14)$coeff[20,2],
summary(M14)$coeff[2,4],summary(M14)$coeff[3,4],summary(M14)$coeff[4,4],
summary(M14)$coeff[5,4],summary(M14)$coeff[6,4], summary(M14)$coeff[7,4],
summary(M14)$coeff[8,4],summary(M14)$coeff[9,4],summary(M14)$coeff[10,4],
summary(M14)$coeff[11,4],summary(M14)$coeff[12,4],summary(M14)$coeff[13,4],
summary(M14)$coeff[14,4],summary(M14)$coeff[15,4],summary(M14)$coeff[16,4],
summary(M14)$coeff[17,4],summary(M14)$coeff[18,4],summary(M14)$coeff[19,4],summary(M14)$coeff[20,4],
M14$null, M14$dev)
#print(fit_results)
#fit_results=as.table(fit.results)
colnames(fit_results)=c("mu","gamma1","gamma1HOM","gamma2","gamma2HOM","gamma3","gamma4","gamma4HOM","gamma5","gamma6","gamma7","gamma8","gamma9","gamma10INT","gamma11INT","gamma12","gamma13NONHLASNP","g14INT0","g14INT1","g14INT2",
"se1","se1HOM","se2","se2HOM","se3","se4","se4HOM","se5","se6","se7","se8","se9","se10","se11","se12","se143NONHLASNP","se14INT0","seINT1","seINT2",
"pval1","pval1HOM","pval2","pval2HOM","pval3","pval4","pval4HOM","pvalue5","pvalue6","pvalue7","pvalue8","pvalue9","pvalue10INT","pvalue11INT","pvalue12SNP","pvalue13NONHLASNP","pval14INT0","pval14INT1","pval14INT2",
"dev_NULL","dev_M")
}else{
fit_results=rep(NA,60)
}
 return (fit_results);

}



# the following is as above, but where instead of stratifying by the IL7R variant, we have three terms based on the copies of the DDX39B SNP.
gwas14_nounc_noG_july17_WITH_3INT_stratified_by_DDX389B = function(interaction_term,status,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,PC){ ## where INT refers to H8 here, which is included as an interaction with DRB*15:01 (presence/absence)
  pH1=H1;
  pH2=H2;
  pH3=H3;
  pH4=H4;
  pH5=H5;
  pH6=H6;
  pH7=H7;
  pH8=H8;
  pH9=H9;
  pH10=H10;
  pH11=H11;
  pH12=H12;
  pH13=H13;
  pH14=H14;

  pH1[pH1>1]=1
  pH2[pH2>1]=1 # in case we get probs > 1 because of rounding
  pH3[pH3>1]=1
  pH4[pH4>1]=1
  pH5[pH5>1]=1
  pH6[pH6>1]=1
  pH7[pH7>1]=1
  pH8[pH8>1]=1
  pH9[pH9>1]=1
  pH10[pH10>1]=1
  pH11[pH11>1]=1
  pH12[pH12>1]=1
  pH13[pH13>1]=1
  pH14[pH14>1]=1



G1=abs(apply(pH1,1,which.max)-3);
G2=abs(apply(pH2,1,which.max)-3);
G3=abs(apply(pH3,1,which.max)-3);
G4=abs(apply(pH4,1,which.max)-3);
G5=abs(apply(pH5,1,which.max)-3);
G6=abs(apply(pH6,1,which.max)-3);
G7=abs(apply(pH7,1,which.max)-3);
G8=abs(apply(pH8,1,which.max)-3);
G9=abs(apply(pH9,1,which.max)-3);
G10=abs(apply(pH10,1,which.max)-3);
G11=abs(apply(pH11,1,which.max)-3);
G12=abs(apply(pH12,1,which.max)-3);
G13=abs(apply(pH13,1,which.max)-3);
G14=abs(apply(pH14,1,which.max)-3);

G1HOM=G1
G1HOM[G1HOM==1]=0
G1HOM[G1HOM==2]=1

I1501=G1
I1501[I1501>1]=1

G2HOM=G2
G2HOM[G2HOM==1]=0
G2HOM[G2HOM==2]=1

G4HOM=G4
G4HOM[G4HOM==1]=0
G4HOM[G4HOM==2]=1

G9DOM=G9
G9DOM[G9DOM>1]=1 #dominant effect

I1=NULL
if(interaction_term==1){
I1=G1
}else if(interaction_term==2){
I1=G2
}else if(interaction_term==3){
I1=G3
}else if(interaction_term==4){
I1=G4HOM
}else if(interaction_term==5){
I1=G5
}else if(interaction_term==6){
I1=G6;
}else if(interaction_term==7){
I1=G7;
}else if(interaction_term==8){
I1=G8;
}else if(interaction_term==9){ #this is what was used for the paper 
I1=G9DOM;
}
I1[I1==2]=1
G10=G10*I1501
G11=G11*I1 # Note that this is the effect of the last HLA allele (not SNP!) included in the data, in the presence/absence of the interactor (DQB*03:02 in the paper) 
# G10 is the new HLA/ SNP to test (where we test all in turn)


#the added bit here: the three interaction terms for the IL7R effect based on the DDX39B genotype

#0 copies
DDX39B_0=G14
DDX39B_0[DDX39B_0==0]=-1
DDX39B_0[DDX39B_0>0]=0
DDX39B_0[DDX39B_0==-1]=1

#1 copy
DDX39B_1=G14
DDX39B_1[DDX39B_1==2]=0

#2 copies
DDX39B_2=G14
DDX39B_2[DDX39B_2==1]=0
DDX39B_2[DDX39B_2==2]=1

#now prepare the three interaction terms
IL7R_INT0=G13*DDX39B_0
IL7R_INT1=G13*DDX39B_1
IL7R_INT2=G13*DDX39B_2

# note that in the G1,G2,...,G5 genotype arrays the number of copies are that of the actual allele. So G[i]=2 means that the SNP or HLA file has two copies of the HLA allele
# or 2 copies of the minor SNP allele.

M14= glm(status~G1+G1HOM+G2+G2HOM+G3+G4+G4HOM+G5+G6+G7+G8+G9DOM+G10+G11+G12+G14+IL7R_INT0+IL7R_INT1+IL7R_INT2+PC,family=binomial) #independent model for DRB*15:01 (G1) and A*02:01 (G2), recessive for DRB*03:01 (4G), additive model for rest

if(dim(summary(M14)$coeff)[1]==25){ 
fit_results = data.frame(summary(M14)$coeff[1,1],
summary(M14)$coeff[2,1],summary(M14)$coeff[3,1],summary(M14)$coeff[4,1],
summary(M14)$coeff[5,1],summary(M14)$coeff[6,1],summary(M14)$coeff[7,1],
summary(M14)$coeff[8,1],summary(M14)$coeff[9,1],summary(M14)$coeff[10,1],
summary(M14)$coeff[11,1],summary(M14)$coeff[12,1],summary(M14)$coeff[13,1],
summary(M14)$coeff[14,1],summary(M14)$coeff[15,1],summary(M14)$coeff[16,1],
summary(M14)$coeff[17,1],summary(M14)$coeff[18,1],summary(M14)$coeff[19,1],summary(M14)$coeff[20,1],
summary(M14)$coeff[2,2],summary(M14)$coeff[3,2],summary(M14)$coeff[4,2],
summary(M14)$coeff[5,2],summary(M14)$coeff[6,2], summary(M14)$coeff[7,2],
summary(M14)$coeff[8,2],summary(M14)$coeff[9,2],summary(M14)$coeff[10,2],
summary(M14)$coeff[11,2],summary(M14)$coeff[12,2],summary(M14)$coeff[13,2],
summary(M14)$coeff[14,2],summary(M14)$coeff[15,2],summary(M14)$coeff[16,2],
summary(M14)$coeff[17,2],summary(M14)$coeff[18,2],summary(M14)$coeff[19,2],summary(M14)$coeff[20,2],
summary(M14)$coeff[2,4],summary(M14)$coeff[3,4],summary(M14)$coeff[4,4],
summary(M14)$coeff[5,4],summary(M14)$coeff[6,4], summary(M14)$coeff[7,4],
summary(M14)$coeff[8,4],summary(M14)$coeff[9,4],summary(M14)$coeff[10,4],
summary(M14)$coeff[11,4],summary(M14)$coeff[12,4],summary(M14)$coeff[13,4],
summary(M14)$coeff[14,4],summary(M14)$coeff[15,4],summary(M14)$coeff[16,4],
summary(M14)$coeff[17,4],summary(M14)$coeff[18,4],summary(M14)$coeff[19,4],summary(M14)$coeff[20,4],
M14$null, M14$dev)
#print(fit_results)
#fit_results=as.table(fit.results)
colnames(fit_results)=c("mu","gamma1","gamma1HOM","gamma2","gamma2HOM","gamma3","gamma4","gamma4HOM","gamma5","gamma6","gamma7","gamma8","gamma9","gamma10INT","gamma11INT","gamma12","gamma14HLASNP","gIL7R_INT0","gIL7R_INT1","gIL7R_INT2",
"se1","se1HOM","se2","se2HOM","se3","se4","se4HOM","se5","se6","se7","se8","se9","se10","se11","se12","se14HLASNP","se13_IL7R_INT0","se13_IL7R_INT1","se13_IL7R_INT2",
"pval1","pval1HOM","pval2","pval2HOM","pval3","pval4","pval4HOM","pvalue5","pvalue6","pvalue7","pvalue8","pvalue9","pvalue10INT","pvalue11INT","pvalue12SNP","pvalue14HLASNP","pval13IL7R_INT0","pval13IL7R_INT1","pval13IL7R_INT2",
"dev_NULL","dev_M")
}else{
fit_results=rep(NA,60)
}
 return (fit_results);

}

