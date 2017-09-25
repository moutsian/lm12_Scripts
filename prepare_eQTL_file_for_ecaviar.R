## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"


# Plan: the LD file to be used will come from the same reference for both files. What will change are the two .z files to be produced, which will be aligned to the .ld file
# The result will be a truncated .ld file for the intersection and two corresponding .z files, one for the GWAS dataset and one for the eQTL. 

args = commandArgs(trailingOnly=TRUE)
if (length(args)<5) {
  stop("Four arguments must be supplied", call.=FALSE)
}

chromosome=args[1]
starting_pos=args[2]
ending_pos=args[3]
celltype=args[4]
funcassay=args[5] #e.g. "gene_nor_combat"
out_dir="/lustre/scratch115/projects/psc2020/eCAVIAR/"
ldfile=paste0("/lustre/scratch115/projects/psc2020/finemapping/bcf.out.",chromosome,".",starting_pos,".",ending_pos,".uk10k_1KG.vcf.ld")
snpstatsfile=paste0("/lustre/scratch115/projects/psc2020/finemapping/bcf.out.",chromosome,".",starting_pos,".",ending_pos,".uk10k_1KG.vcf.annotated.concat.snpstats") #-->essentially all we need from this 
#is the alignment of the two SNPs for the LD file
gwas_sumstats=paste0("/lustre/scratch115/projects/psc2020/splitfiles/pscgwas/psc.all_for_reference.chr",chromosome,".txt")
func_sumstats=paste0("/lustre/scratch115/projects/psc2020/splitfiles/",celltype,"/",celltype,"_",funcassay,"_peer_10_all_summary.Chr_",chromosome,".sorted.txt")

#read files in
LD=read.table(ldfile,head=F)
LDsnps=read.table(snpstatsfile,head=T)
GWAS=read.table(gwas_sumstats,head=T)
FUNC=read.table(func_sumstats,head=F)
colnames(FUNC)=c("chrom","pos","ref_alt","snp","phenotypeID","p_value","beta","bonferroni_p_value","FDR","alt_allele_freq","se")

#firstly prepare the GWAS .z and .ld files
LDtokeep=which(LDsnps$position%in%GWAS$coor)
LDsubset=LD[LDtokeep,LDtokeep] #this should be the same (keeping all) in the case of GWAS data, but not in the case of functional data
LDsnps_subset=LDsnps[LDtokeep,]
gwas_vars=which(GWAS$coor%in%LDsnps_subset$position) #subset to keep
GWASsubset=GWAS[gwas_vars,]
refpanplus=merge(LDsnps_subset,GWASsubset,by.y="coor",by.x="position",all.x) #merge and align them
tmp1=which(as.character(refpanplus$A_allele)!=as.character(refpanplus$a0))
tmp2=which(as.character(refpanplus$A_allele)==as.character(refpanplus$a1))
torevert=intersect(tmp1,tmp2)
refpan_z=cbind(as.character(refpanplus$snp),as.numeric(refpanplus$b)/as.numeric(refpanplus$se))
refpan_z[torevert,2]=as.numeric(refpan_z[torevert,2])*(-1)
outfile=paste0(out_dir,"psc_gwas.chr",chromosome,"_",starting_pos,".",ending_pos,".z")
write.table(refpan_z,outfile,quote=F,col.names=F,row.names=F)
outfile=paste0(out_dir,"psc_gwas.chr",chromosome,"_",starting_pos,".",ending_pos,".ld")
write.table(LDsubset,outfile,quote=F,col.names=F,row.names=F)


#now for the functional dataset, we need to do this per phenotype (eg. per transcript for eQTLs)
all_transcripts=unique(as.character(FUNC$phenotypeID))
for i in 1:length(all_transcripts)){
FUNCsubset=FUNC[as.character(FUNC$phenotypeID)==all_transcripts[i],]

LDtokeep=which(LDsnps$position%in%FUNCsubset$pos)
if(length(LDtokeep)[1]>1){ #not including entries with just one SNP
LDsubset=LD[LDtokeep,LDtokeep]
LDsnps_subset=LDsnps[LDtokeep,]
func_vars=which(FUNCsubset$pos%in%LDsnps$position)
FUNCsubset_final=FUNCsubset[func_vars,]
refpanplus=merge(LDsnps_subset,FUNCsubset_final,by.y="pos",by.x="position",all.x) #merge and align them
#prepare the ref and alt alleles in a good format for comparison
tmp=matrix(unlist(strsplit(as.character(refpanplus$ref_alt),"_")))
a0=tmp[seq(from=1,to=length(tmp),by=2)]
a1=tmp[seq(from=2,to=length(tmp),by=2)]
tmp1=which(as.character(refpanplus$A_allele)!=a0)
tmp2=which(as.character(refpanplus$A_allele)==a1)
torevert=intersect(tmp1,tmp2)
refpan_z=cbind(as.character(refpanplus$snp),as.numeric(refpanplus$beta)/as.numeric(refpanplus$se))
if(length(torevert)>0){
refpan_z[torevert,2]=as.numeric(refpan_z[torevert,2])*(-1)
}
}
#now write the z and ld files
outfile=paste0(out_dir,celltype,"_",funcassay,"_",chromosome,".",all_transcripts[i],".z")
write.table(refpan_z,outfile,quote=F,col.names=F,row.names=F)
outfile=paste0(out_dir,celltype,"_",funcassay,"_",chromosome,".",all_transcripts[i],".ld")
write.table(LDsubset,outfile,quote=F,col.names=F,row.names=F)

}



#END
