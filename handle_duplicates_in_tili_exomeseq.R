#deal with duplicates in the exome seq data (bcf file  - tili project)
#allvariants=read.table("tili.poly.biallelic.nomiss.recode.rs1.variants",stringsAsFactors=F) #these are all the variants in the bcf file
allvariants=read.table("tili.poly.biallelic.qc1.T.nodup.variants",stringsAsFactors=F)  #this is to see if it has worked

varstokeep=read.table("snps_to_keep_after_hwe_check.txt",stringsAsFactors=F) #this is the list of the variants to keep, after checking for HWE deviations.

#for now I will remove ALL duplicates, since this means that for the same SNP we have more than two alleles (sometimes just a triallelic SNPs, sometimes copy number variations.)
n_occur=data.frame(table(allvariants[,3]))
dupvars=n_occur[n_occur$Freq>1,]
TOREMOVE=matrix(ncol=1,nrow=dim(varstokeep)[1],0)
for(i in 1:dim(dupvars)[1])
{
idx=which(allvariants[,3]==dupvars[i,1])
#now match to the varstokeep 
v_idx=which(varstokeep[,1]==allvariants[idx[1],1] & varstokeep[,2]==allvariants[idx[1],2] )
TOREMOVE[v_idx,1]=1 #flag for removal
}

#second round where we remove variants in the same position but where the two alleles have been reversed
n_occur=data.frame(table(allvariants[,2])) 
dupvars=n_occur[n_occur$Freq>1,] #note the dups by position here and not by ID.

#this should be faster
v_idx=NULL
for(i in 1:dim(dupvars)[1])
{
if(i%%100==0){
print(paste(i," out of ",dim(dupvars)[1],sep=""))}
idx=which(allvariants[,2]==dupvars[i,1])
#now match to the varstokeep 
for (cnt in 1:length(idx)){
v_idx=c(v_idx,which(varstokeep[,1]==allvariants[idx[cnt],1] & varstokeep[,2]==allvariants[idx[cnt],2] ))
}}
TOREMOVE[as.numeric(names(which(table(v_idx)>1))),1]=1#flag for removal


varstokeep_nodup=varstokeep[which(TOREMOVE[,1]==0),]
varstokeep_nodup_unique=unique(varstokeep_nodup)

write.table(varstokeep_nodup_unique,"snps_to_keep_after_hwe_check_nodup_updated.txt",sep="\t",quote=F,row.names=F,col.names=F)

#END
