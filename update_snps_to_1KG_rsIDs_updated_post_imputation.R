#updating the post-imputation variants to rsID variants
# this script should be fine with handling triplicates and more (i.e. not just duplicates like the earlier version)
# However, it will still keep entries for which there are more than one potential rsIDs with the same length of alleles, e.g. a case where one rsID has A,C as alleles and the other A,G or A,. 
options(stringsAsFactors = FALSE)
dataset="ichip" #e.g. gwas1
DATASET="ichip" #e.g. GWAS1
for(CHROM in 22:1){
vcfCHROM=read.table(paste("/lustre/scratch115/projects/ibdgwas/new_imputation/",DATASET,"M25.vcfs/",CHROM,".",dataset,".tokeep.txt",sep=""), stringsAsFactors=FALSE) #the alleles from the vcf file
vcfCHROMout=paste("/lustre/scratch115/projects/ibdgwas/new_imputation/",DATASET,"M25.vcfs/",CHROM,".",dataset,".tokeep.updatedIDs.txt",sep="") 
kgCHROM=read.table(paste("/lustre/scratch115/projects/ibdgwas/aux_files/chr",CHROM,".1KG_b37.alleles",sep=""), stringsAsFactors=FALSE) #chr22 of the 1KG file (loading the whole variant set from 1KG could take long)
colnames(kgCHROM)=c("chr","pos","rsID","al1","al2")
colnames(vcfCHROM)=c("chr","pos","id","al1vcf","al2vcf","AF","other_info")
togetherCHROM=merge(vcfCHROM,kgCHROM,by.x="pos",by.y="pos",all.x=T)
#only replace the IDs for the ones who don't already have an rsID
toreplace=which((substr(as.character(togetherCHROM[,3]),start=1,stop=2)=="rs")==F)
togetherCHROMs=togetherCHROM[toreplace,]
#for the positions where the two alleles are the same, we can immediately update the SNP IDs.
allele_match=which(as.character(togetherCHROMs$al1vcf)==as.character(togetherCHROMs$al1) & as.character(togetherCHROMs$al2vcf)==as.character(togetherCHROMs$al2))
#for many of the positions  for which this is not the case, this is because the variant is triallelic in 1KG. In these cases, we check if the two variants are 
#part of the variants in the vcf file and, if so, we update the SNP IDs too.
allele_mismatch=which(!(as.character(togetherCHROMs$al1vcf)==as.character(togetherCHROMs$al1) & as.character(togetherCHROMs$al2vcf)==as.character(togetherCHROMs$al2)))
allele_match2=NULL
for(i in 1:length(allele_mismatch)){
AL1=grepl(togetherCHROMs[allele_mismatch[i],4],togetherCHROMs[allele_mismatch[i],10])
AL2=grepl(togetherCHROMs[allele_mismatch[i],5],togetherCHROMs[allele_mismatch[i],11]) 
#This can create issues - for instance I had a variant with alleles C/T which was wrongly assigned the same RSID as a variant with alleles CCCC/CTCT
#It could be solved by matching the first allele exactly, unless we have REF/ALT switch. For now I've left it as is.
AL1b=grepl(togetherCHROMs[allele_mismatch[i],4],togetherCHROMs[allele_mismatch[i],11])
AL2b=grepl(togetherCHROMs[allele_mismatch[i],5],togetherCHROMs[allele_mismatch[i],10])
if((AL1==T & AL2==T)|(AL1b==T & AL2b==T) ){ 
allele_match2=c(allele_match2,i)#note that this is index to allele_mismatch
}
}
#based on the above, add these to the allele_match index
allele_match_full=c(allele_match,allele_mismatch[allele_match2])
togetherCHROMs[allele_match_full,3]=as.character(togetherCHROMs[allele_match_full,9])
togetherCHROMs[-allele_match_full,3]=paste(togetherCHROMs[-allele_match_full,2],togetherCHROMs[-allele_match_full,1],sep="_")
write.table(togetherCHROMs[,c(2,1,3:7)],vcfCHROMout,quote=F,col.names=F,row.names=F,sep="\t") #this is the one that will be used to replace the IDs
}

#END