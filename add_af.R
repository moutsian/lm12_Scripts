#MAR2017

complement <- function(allele)
{
reval=NULL
if(allele=="C"){
reval="G"}else 
if(allele=="G"){
reval="C"}else
if(allele=="T"){
reval="A"}else
if(allele=="A"){
reval="T"
}else{
reval="IGNORE" #
}
return(reval)
}

#chr=16
for(chr in 2:22){
print(paste0("Running for chrom: ",chr))
infile=paste0("t1dmeta_all_summary.Chr_",chr,".sorted.txt",sep="")
t1d=read.table(infile,head=T,sep=" ")
psc=read.table(paste0("../pscgwas/psc.all_for_reference.chr",chr,".txt"),head=T)
idx_shared=which(t1d$Position%in%psc$coor)
t1dfull=t1d #just to keep if we need it, but work downstream with the overlap only
t1d=t1dfull[idx_shared,]
idx_shared=which(t1d$Position%in%psc$coor)
freq=matrix(ncol=2,nrow=length(idx_shared),NA)
for(i in 1:length(idx_shared)){
idx_psc=which(psc$coor==t1d$Position[idx_shared[i]])
#we want the frequency of the MinorAllele, since this is the allele we are reporting the OR for.
	if( as.character(psc$a1[idx_psc])==as.character(t1d$MinorAllele[idx_shared[i]]) && as.character(psc$a0[idx_psc])==as.character(t1d$MajorAllele[idx_shared[i]]) ){
	freq[i,1]=psc$freq_a1[idx_psc]
	freq[i,2]=1
	}else if( as.character(psc$a1[idx_psc])==as.character(t1d$MajorAllele[idx_shared[i]]) && as.character(psc$a0[idx_psc])==as.character(t1d$MinorAllele[idx_shared[i]]) ){
	freq[i,1]=1-psc$freq_a1[idx_psc]
	freq[i,2]=2
	}else if (complement(as.character(psc$a1[idx_psc]))==as.character(t1d$MinorAllele[idx_shared[i]]) && complement(as.character(psc$a0[idx_psc]))==as.character(t1d$MajorAllele[idx_shared[i]]) ){ #complement
	freq[i,1]=psc$freq_a1[idx_psc]
	freq[i,2]=3	
	}else if ( complement(as.character(psc$a1[idx_psc]))==as.character(t1d$MajorAllele[idx_shared[i]]) && complement(as.character(psc$a0[idx_psc]))==as.character(t1d$MinorAllele[idx_shared[i]]) ){ #reverse complement
	freq[i,1]=1-psc$freq_a1[idx_psc]
	freq[i,2]=4
	}
	else{
	#if alleles don't match, just leave out --this happens because the two alleles are missing for some of the T1D entries.
	}
}

#now update T1D file with frequencies:
t1d=cbind(t1d,freq[,1])
colnames(t1d)[5]="OR_MinorAllele"
colnames(t1d)[9]="Freq_MinorAllele"
outfile=paste0(infile,".ready")
write.table(t1d,outfile,quote=F,col.names=T,row.names=F,sep="\t")
}
#END