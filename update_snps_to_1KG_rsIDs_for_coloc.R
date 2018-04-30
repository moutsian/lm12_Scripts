#updating the post-imputation variants to rsID variants
# this script should be fine with handling triplicates and more (i.e. not just duplicates like the earlier version)
# However, it will still keep entries for which there are more than one potential rsIDs with the same length of alleles, e.g. a case where one rsID has A,C as alleles and the other A,G or A,. 
options(stringsAsFactors = FALSE)
trait="CD" #e.g. UC or CD
COLOC_FILE=NULL

for(CHROM in 22:1){
if(trait=="UC"){
COLOCFILE=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/",tolower(trait),"_b37_45975.Chr_",CHROM,".sorted.txt",sep=""), stringsAsFactors=FALSE,head=T)
COLOCFILEout=paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/",tolower(trait),"_b37_45975.Chr_",CHROM,".sorted.updatedIDs.txt",sep="")
}else{ #trait="CD"
COLOCFILE=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/",tolower(trait),"_b37_40266.Chr_",CHROM,".sorted.txt",sep=""), stringsAsFactors=FALSE,head=T)
COLOCFILEout=paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/",tolower(trait),"_b37_40266.Chr_",CHROM,".sorted.updatedIDs.txt",sep="")
}
kgCHROM=read.table(paste("/lustre/scratch115/projects/ibdgwas/aux_files/chr",CHROM,".1KG_b37.alleles",sep=""), stringsAsFactors=FALSE) #chr22 of the 1KG file (loading the whole variant set from 1KG could take long)
colnames(kgCHROM)=c("chr","pos","rsID","al1","al2")
togetherCHROM=merge(COLOCFILE,kgCHROM,by.x="Position",by.y="pos",all.x=T)
togetherCHROMs=togetherCHROM
#for the positions where the two alleles are the same, we can immediately update the SNP IDs.
allele_match=which(tolower(as.character(togetherCHROMs$al1))==tolower(as.character(togetherCHROMs$Al1)) & tolower(as.character(togetherCHROMs$al2))==tolower(as.character(togetherCHROMs$Al2)))
#for many of the positions  for which this is not the case, this is because the variant is triallelic in 1KG. In these cases, we check if the two variants are 
#part of the variants in the vcf file and, if so, we update the SNP IDs too.
allele_mismatch=which(!(tolower(as.character(togetherCHROMs$al1))==tolower(as.character(togetherCHROMs$Al1)) & tolower(as.character(togetherCHROMs$al2))==tolower(as.character(togetherCHROMs$Al2))))
allele_match2=NULL
for(i in 1:length(allele_mismatch)){
AL1=grepl(tolower(togetherCHROMs[allele_mismatch[i],7]),tolower(togetherCHROMs[allele_mismatch[i],12]))
AL2=grepl(tolower(togetherCHROMs[allele_mismatch[i],8]),tolower(togetherCHROMs[allele_mismatch[i],13])) 
#This can create issues - for instance I had a variant with alleles C/T which was wrongly assigned the same RSID as a variant with alleles CCCC/CTCT
#It could be solved by matching the first allele exactly, unless we have REF/ALT switch. For now I've left it as is.
AL1b=grepl(tolower(togetherCHROMs[allele_mismatch[i],7]),tolower(togetherCHROMs[allele_mismatch[i],13]))
AL2b=grepl(tolower(togetherCHROMs[allele_mismatch[i],8]),tolower(togetherCHROMs[allele_mismatch[i],12])) 
if((AL1==T & AL2==T)|(AL1b==T & AL2b==T) ){ 
allele_match2=c(allele_match2,i)#note that this is index to allele_mismatch
}
}
#based on the above, add these to the allele_match index
allele_match_full=c(allele_match,allele_mismatch[allele_match2])
togetherCHROMs[allele_match_full,3]=as.character(togetherCHROMs[allele_match_full,11])
togco=togetherCHROMs[order(togetherCHROMs[,1],decreasing=F),]
#togetherCHROMs[-allele_match_full,3]=paste(togetherCHROMs[-allele_match_full,2],togetherCHROMs[-allele_match_full,1],sep="_")
write.table(togco[,c(2,1,3:9)],COLOCFILEout,quote=F,col.names=T,row.names=F,sep="\t") #this is the one that will be used to replace the IDs
}

#END