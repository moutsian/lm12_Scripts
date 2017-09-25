# June 2012
## This function is for 5-parameter logistic regression incorporating uncertainty on posteriors  
## It also plots 3x1 plots of results.


#example:
#case_file="snps.cases.ms2012.UK"
#control_file="snps.cases.ms2012.UK"
# locus="DPB"
# allele_res="4D"
# cond_case_file="DRB.cases.ms2012.UK"
# cond_control_file="DRB.ctrls.ms2012.UK"
# cond_locus="DRB"
# cond_type="1501"
# cond_allele_res = "4D"
# cond_case_file2="A.cases.ms2012.UK"
# cond_control_file2="A.ctrls.ms2012.UK"
# cond_locus2="A"
# cond_type2="0201"
# cond_allele_res2 = "4D"
# cond_case_file3="DRB.cases.ms2012.UK"
# cond_control_file3="DRB.ctrls.ms2012.UK"
# cond_locus3="DRB"
# cond_type3="1303"
# cond_allele_res3 = "4D"
# cond_case_file4="DRB.cases.ms2012.UK"
# cond_control_file4="DRB.ctrls.ms2012.UK"
# cond_locus4="DRB"
# cond_type4="0301"
# cond_allele_res4 = "4D"
# results_output ="test.out"
#alelle="0301"
logreg_14p_nounc_nogender_sept16_SNP_INT = function(case_file,control_file,
						cond_case_file,cond_control_file,cond_locus,cond_type,cond_allele_res,
						cond_case_file2,cond_control_file2,cond_locus2,cond_type2,cond_allele_res2,
						cond_case_file3,cond_control_file3,cond_locus3,cond_type3,cond_allele_res3,
						cond_case_file4,cond_control_file4,cond_locus4,cond_type4,cond_allele_res4,
						cond_case_file5,cond_control_file5,cond_locus5,cond_type5,cond_allele_res5,
						cond_case_file6,cond_control_file6,cond_locus6,cond_type6,cond_allele_res6,
						cond_case_file7,cond_control_file7,cond_locus7,cond_type7,cond_allele_res7,
						cond_case_file8,cond_control_file8,cond_locus8,cond_type8,cond_allele_res8,
						cond_case_file9,cond_control_file9,cond_locus9,cond_type9,cond_allele_res9,
						cond_case_file10,cond_control_file10,cond_locus10,cond_type10,cond_allele_res10,
						cond_case_file11,cond_control_file11,cond_locus11,cond_type11,cond_allele_res11,
						cond_snp_casefile,cond_snp_ctrlfile,cond_SNP,cond_SNP2,interaction_term,
						results_output,pc_cases,pc_ctrls,from,to,DATADIR,RESDIR,SCRIPTDIR){


#first, load functions for 12 par log reg. Main function: gwas8_nounc_noG_oct14_INT
source(paste(SCRIPTDIR,"GWAS_14p_nounc_nogender_sept16.R",sep=""))
source(paste(SCRIPTDIR,"prep_alleleprobs2012.R",sep=""))
#setwd(DATADIR)
	

#load PCs
PC_CASES=read.table(pc_cases,sep=" ",head=T);
PC_CTRLS=read.table(pc_ctrls,sep=" ",head=T);
PCMATRIX=as.matrix(rbind(PC_CTRLS[,-1],PC_CASES[,-1])); # first column contains sample ID.
print(paste("PC_CASES:",dim(PC_CASES)," PC_CTRLS:",dim(PC_CTRLS)," PC: ",dim(PCMATRIX),sep=""))
	
						
CASES = read.table(case_file,head=T);
CTRLS = read.table(control_file,head=T);
if(dim(CASES)[2]!=dim(CTRLS)[2]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
CASE_SAMPLES=CASES[,1]
CTRL_SAMPLES=CTRLS[,1]
CASES=CASES[,-1] #to remove sample IDs and keep the SNP data
CTRLS=CTRLS[,-1]
SNP_IDs=colnames(CASES);
RISKALLELE=matrix(ncol=1,nrow=length(SNP_IDs),NA);
OTHERALLELE=matrix(ncol=1,nrow=length(SNP_IDs),NA);


#for the SNP we are conditioning on:
CCASES = read.table(cond_snp_casefile,head=T);
CCTRLS = read.table(cond_snp_ctrlfile,head=T);
if(dim(CCASES)[2]!=dim(CCTRLS)[2]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
CCASE_SAMPLES=CCASES[,1]
CCTRL_SAMPLES=CCTRLS[,1]
CCASES=CCASES[,-1] #to remove sample IDs and keep the SNP data
CCTRLS=CCTRLS[,-1]
CSNP_IDs=colnames(CCASES);
snp_idx=which(CSNP_IDs==cond_SNP);
if(length(snp_idx)!=1){
stop(paste("Either NO or multiple SNPs with ID ",cond_SNP," were found - Aborting.",sep=""));
}
CDATA=cbind(t(as.character(CCTRLS[,snp_idx])),t(as.character(CCASES[,snp_idx])));
CSNP=t(CDATA)
CSNP=gsub('\\?','NA',CSNP);
CALLELES=unique(unlist(strsplit(CSNP,"/")))
CSNPDATA=matrix(ncol=3,nrow=length(CSNP),0)
if(length(CALLELES)==1)#monomorphic
{
stop(paste("SNP ",cond_SNP," is monomorphic and we cannot condition upon it. Aborting.",sep=""));
}
else{
myallele="NA";
mylabels=c(1,0);
if(CALLELES[1]=="A" || CALLELES[2]=="A"){
	myallele="A";
	if(CALLELES[2]=="A"){
		mylabels=c(0,1);
	}
}else 
if(CALLELES[1]=="C" || CALLELES[2]=="C"){
	myallele="C";
	if(CALLELES[2]=="C"){
		mylabels=c(0,1);
	}

}else{
myallele="G";
if(CALLELES[2]=="G"){
		mylabels=c(0,1);
	}

}
if(length(CALLELES)>2){
which(CALLELES=="NA");
if(which(CALLELES=="NA")==3){
mylabels=c(mylabels,"NA");
}else
if(which(CALLELES=="NA")==2){
mylabels=c(mylabels[1],"NA",mylabels[2]);
}else
if(which(CALLELES=="NA")==1){
mylabels=c("NA",mylabels);
}
}
FULL=unlist(strsplit(CSNP,"/"))
FULL=factor(FULL,labels = mylabels)
FULL1=abs(as.numeric(FULL[seq(from=1,to=length(FULL),by=2)])-2)
FULL2=abs(as.numeric(FULL[seq(from=2,to=length(FULL),by=2)])-2)
GENO=FULL1+FULL2
CSNPDATA[GENO==1,2]=1
CSNPDATA[GENO==2,3]=1
CSNPDATA[GENO==0,1]=1
}


#for the 2nd SNP we are conditioning on:
snp_idx2=which(CSNP_IDs==cond_SNP2);
if(length(snp_idx2)!=1){
stop(paste("Either NO or multiple SNPs with ID ",cond_SNP2," were found - Aborting.",sep=""));
}
CDATA=cbind(t(as.character(CCTRLS[,snp_idx2])),t(as.character(CCASES[,snp_idx2])));
CSNP=t(CDATA)
CSNP=gsub('\\?','NA',CSNP);
CALLELES=unique(unlist(strsplit(CSNP,"/")))
CSNPDATA2=matrix(ncol=3,nrow=length(CSNP),0)
if(length(CALLELES)==1)#monomorphic
{
stop(paste("SNP ",cond_SNP2," is monomorphic and we cannot condition upon it. Aborting.",sep=""));
}
else{
myallele="NA";
mylabels=c(1,0);
if(CALLELES[1]=="A" || CALLELES[2]=="A"){
	myallele="A";
	if(CALLELES[2]=="A"){
		mylabels=c(0,1);
	}
}else 
if(CALLELES[1]=="C" || CALLELES[2]=="C"){
	myallele="C";
	if(CALLELES[2]=="C"){
		mylabels=c(0,1);
	}

}else{
myallele="G";
if(CALLELES[2]=="G"){
		mylabels=c(0,1);
	}

}
if(length(CALLELES)>2){
which(CALLELES=="NA");
if(which(CALLELES=="NA")==3){
mylabels=c(mylabels,"NA");
}else
if(which(CALLELES=="NA")==2){
mylabels=c(mylabels[1],"NA",mylabels[2]);
}else
if(which(CALLELES=="NA")==1){
mylabels=c("NA",mylabels);
}
}
FULL=unlist(strsplit(CSNP,"/"))
FULL=factor(FULL,labels = mylabels)
FULL1=abs(as.numeric(FULL[seq(from=1,to=length(FULL),by=2)])-2)
FULL2=abs(as.numeric(FULL[seq(from=2,to=length(FULL),by=2)])-2)
GENO=FULL1+FULL2
CSNPDATA2[GENO==1,2]=1
CSNPDATA2[GENO==2,3]=1
CSNPDATA2[GENO==0,1]=1
}

if(to==-1){
to = length(SNP_IDs);
}
# the status refers to genotypes, not haplotypes
status = c(rep(0,dim(CTRLS)[1]),rep(1,dim(CASES)[1]))

#now for the allele to condition on:
CNDCASES = prep_alleleprobs2012(cond_case_file,cond_allele_res);
CNDCTRLS = prep_alleleprobs2012(cond_control_file,cond_allele_res);
if(dim(CNDCASES)[3]!=dim(CNDCTRLS)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type)
if(length(coltopick)==0){
stop("There is no allele with this type in the conditional locus. Aborting...\n\n");
}

COND_DATA=rbind(CNDCTRLS[,,coltopick],CNDCASES[,,coltopick])


#now for the 2nd allele to condition on:
CNDCASES2 = prep_alleleprobs2012(cond_case_file2,cond_allele_res2);
CNDCTRLS2 = prep_alleleprobs2012(cond_control_file2,cond_allele_res2);
if(dim(CNDCASES2)[3]!=dim(CNDCTRLS2)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES2)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type2)
if(length(coltopick)==0){
stop("There is no allele with this type in the 2nd conditional locus. Aborting...\n\n");
}
COND_DATA2=rbind(CNDCTRLS2[,,coltopick],CNDCASES2[,,coltopick])


#now for the 3rd allele to condition on:
CNDCASES3 = prep_alleleprobs2012(cond_case_file3,cond_allele_res3);
CNDCTRLS3 = prep_alleleprobs2012(cond_control_file3,cond_allele_res3);
if(dim(CNDCASES3)[3]!=dim(CNDCTRLS3)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES3)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type3)
if(length(coltopick)==0){
stop("There is no allele with this type in the 3rd conditional locus. Aborting...\n\n");
}
COND_DATA3=rbind(CNDCTRLS3[,,coltopick],CNDCASES3[,,coltopick])

#now for the 4th allele to condition on:
CNDCASES4 = prep_alleleprobs2012(cond_case_file4,cond_allele_res4);
CNDCTRLS4 = prep_alleleprobs2012(cond_control_file4,cond_allele_res4);
if(dim(CNDCASES4)[3]!=dim(CNDCTRLS4)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES4)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type4)
if(length(coltopick)==0){
stop("There is no allele with this type in the 4th conditional locus. Aborting...\n\n");
}
COND_DATA4=rbind(CNDCTRLS4[,,coltopick],CNDCASES4[,,coltopick])



#now for the 5th allele to condition on:
CNDCASES5 = prep_alleleprobs2012(cond_case_file5,cond_allele_res5);
CNDCTRLS5 = prep_alleleprobs2012(cond_control_file5,cond_allele_res5);
if(dim(CNDCASES5)[3]!=dim(CNDCTRLS5)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES5)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type5)
if(length(coltopick)==0){
stop("There is no allele with this type in the 5th conditional locus. Aborting...\n\n");
}
COND_DATA5=rbind(CNDCTRLS5[,,coltopick],CNDCASES5[,,coltopick])


#now for the 6th allele to condition on:
CNDCASES6 = prep_alleleprobs2012(cond_case_file6,cond_allele_res6);
CNDCTRLS6 = prep_alleleprobs2012(cond_control_file6,cond_allele_res6);
if(dim(CNDCASES6)[3]!=dim(CNDCTRLS6)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES6)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type6)
if(length(coltopick)==0){
stop("There is no allele with this type in the 6th conditional locus. Aborting...\n\n");
}
COND_DATA6=rbind(CNDCTRLS6[,,coltopick],CNDCASES6[,,coltopick])

#now for the 7th allele to condition on:
CNDCASES7 = prep_alleleprobs2012(cond_case_file7,cond_allele_res7);
CNDCTRLS7 = prep_alleleprobs2012(cond_control_file7,cond_allele_res7);
if(dim(CNDCASES7)[3]!=dim(CNDCTRLS7)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES7)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type7)
if(length(coltopick)==0){
stop("There is no allele with this type in the 7th conditional locus. Aborting...\n\n");
}
COND_DATA7=rbind(CNDCTRLS7[,,coltopick],CNDCASES7[,,coltopick])

#now for the 8th allele to condition on:
CNDCASES8 = prep_alleleprobs2012(cond_case_file8,cond_allele_res8);
CNDCTRLS8 = prep_alleleprobs2012(cond_control_file8,cond_allele_res8);
if(dim(CNDCASES8)[3]!=dim(CNDCTRLS8)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES8)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type8)
if(length(coltopick)==0){
stop("There is no allele with this type in the 8th conditional locus. Aborting...\n\n");
}
COND_DATA8=rbind(CNDCTRLS8[,,coltopick],CNDCASES8[,,coltopick])

#now for the 9th allele to condition on:
CNDCASES9 = prep_alleleprobs2012(cond_case_file9,cond_allele_res9);
CNDCTRLS9 = prep_alleleprobs2012(cond_control_file9,cond_allele_res9);
if(dim(CNDCASES9)[3]!=dim(CNDCTRLS9)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES9)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type9)
if(length(coltopick)==0){
stop("There is no allele with this type in the 9th conditional locus. Aborting...\n\n");
}
COND_DATA9=rbind(CNDCTRLS9[,,coltopick],CNDCASES9[,,coltopick])


#now for the 10nth allele to condition on:
CNDCASES10 = prep_alleleprobs2012(cond_case_file10,cond_allele_res10);
CNDCTRLS10 = prep_alleleprobs2012(cond_control_file10,cond_allele_res10);
if(dim(CNDCASES10)[3]!=dim(CNDCTRLS10)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES10)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type10)
if(length(coltopick)==0){
stop("There is no allele with this type in the 10nth conditional locus. Aborting...\n\n");
}
COND_DATA10=rbind(CNDCTRLS10[,,coltopick],CNDCASES10[,,coltopick])


#now for the 11nth allele to condition on:
CNDCASES11 = prep_alleleprobs2012(cond_case_file11,cond_allele_res11);
CNDCTRLS11 = prep_alleleprobs2012(cond_control_file11,cond_allele_res11);
if(dim(CNDCASES11)[3]!=dim(CNDCTRLS11)[3]){
stop("Different number of alleles in cases and controls. Aborting...\n\n");
}
cond_unique_alleles=unlist(dimnames(CNDCASES11)[3]) #all unique alleles in the dataset.
coltopick = which(cond_unique_alleles==cond_type11)
if(length(coltopick)==0){
stop("There is no allele with this type in the 11nth conditional locus. Aborting...\n\n");
}
COND_DATA11=rbind(CNDCTRLS11[,,coltopick],CNDCASES11[,,coltopick])


#time to perform the logistic regression test
log_reg_results = matrix(nrow=to-from+1,ncol=51)
for(i in from:to){
cat(i," ");
	DATA=cbind(t(as.character(CTRLS[,i])),t(as.character(CASES[,i])));
	SNP=t(DATA)
	SNP=gsub('\\?','NA',SNP);
	ALLELES=unique(unlist(strsplit(SNP,"/")))
	if(is.na(ALLELES[1])){
	ALLELES[1]=ALLELES[2]
	ALLELES[2]=ALLELES[3]
	ALLELES[3]=NA
	}else if(is.na(ALLELES[2])){
	ALLELES[2]=ALLELES[3]
	ALLELES[3]=NA	
	}else{}
	if(length(ALLELES)==1)#monomorphic
	{
		cat(i," mono ");
		log_reg_results[i-from+1,]=rep(-1,dim(log_reg_results)[2])		
	}
	else{
	cat(i," poly  ");

	# I count the copies of A,C,G,T in that order of priority.
	# E.g. if the two alleles are A and C, the genotype contains the counts for A,
	# whereas if the two alleles are C and G, the genotype contains the counts for C.
	# This is for consistency across cohorts without using info on the reference allele
	#I assume no differences in strand between cohorts.
	myallele="NA";
	mylabels=c(1,0);
	if(ALLELES[1]=="A" || ALLELES[2]=="A"){
		myallele="A";
		if(ALLELES[2]=="A"){
			mylabels=c(0,1);
		}
	}else 
	if(ALLELES[1]=="C" || ALLELES[2]=="C"){
		myallele="C";
		if(ALLELES[2]=="C"){
			mylabels=c(0,1);
		}

	}else{
	myallele="G";
	if(ALLELES[2]=="G"){
			mylabels=c(0,1);
		}

	}
	# if(length(ALLELES)>2){
	# cat(i," NA entries ");
	# which(is.na(ALLELES));
	# if(which(is.na(ALLELES))==3){
	# mylabels=c(mylabels,"NA");
	# }else
	# if(which(is.na(ALLELES))==2){
	# mylabels=c(mylabels[1],"NA",mylabels[2]);
	# }else
	# if(which(is.na(ALLELES))==1){
	# mylabels=c("NA",mylabels);
	# }
	# }
	RISKALLELE[i]=ALLELES[mylabels[1]+1]
	OTHERALLELE[i]=ALLELES[mylabels[2]+1]
	SNP[is.na(SNP)]="NA/NA"
	FULL=unlist(strsplit(SNP,"/"))
	FULL[FULL=="NA"]=NA
	FULL=factor(FULL,labels = mylabels,exclude=NA)
	FULL1=abs(as.numeric(FULL[seq(from=1,to=length(FULL),by=2)])-2)
	FULL2=abs(as.numeric(FULL[seq(from=2,to=length(FULL),by=2)])-2)
	GENO=FULL1+FULL2
	DATA=matrix(ncol=3,nrow=length(GENO),0)
	DATA[GENO==1,2]=1
	DATA[GENO==2,3]=1
	DATA[GENO==0,1]=1
	log_reg_results[i-from+1,]=as.matrix(gwas14_nounc_noG_sept16_INT(interaction_term,status,COND_DATA,COND_DATA2,COND_DATA3,COND_DATA4,COND_DATA5,COND_DATA6,COND_DATA7,COND_DATA8,COND_DATA9,COND_DATA10,COND_DATA11,CSNPDATA,CSNPDATA2,DATA,PCMATRIX))
	} 
}
cat("\n");

rownames(log_reg_results)= SNP_IDs[from:to]
colnames(log_reg_results)=c("mu","gamma1","gamma1HOM","gamma2","gamma2HOM","gamma3","gamma4REC","gamma5","gamma6","gamma7","gamma8","gamma9","gamma10INT","gamma11INT","gamma12","gamma13","gamma14",
"se1","se1HOM","se2","se2HOM","se3","se4REC","se5","se6","se7","se8","se9","se10","se11","se12","se13","se14",
"pval1","pval1HOM","pval2","pval2HOM","pval3","pval4REC","pvalue5","pvalue6","pvalue7","pvalue8","pvalue9","pvalue10INT","pvalue11INT","pvalue12SNP","pvalue13SNP","pvalue14",
"dev_NULL","dev_M")
##no plotting inherent to this function. Plot separately.						
									
 RES = data.frame(allele=SNP_IDs[from:to],RISKALLELE[from:to],OTHERALLELE[from:to],log_reg_results)
 write.table(RES,paste(RESDIR,results_output,sep=""),quote=F,row.names=F,sep="\t",col.names=T)


}#end of function



