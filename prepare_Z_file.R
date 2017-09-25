
## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Three arguments must be supplied", call.=FALSE)
}
gwasfile=args[1] #example: rs3748816.variant.full.list.txt
freqfile= args[2] #example bcf.out.4.122999745.123999745.uk10k_1KG.vcf.annotated.concat.snpstats ### i get this from qctool
refpan=args[3] #example: bcf.out.1.2026746.3026746.uk10k_1KG.vcf.annotated.concat.snplist

#GWAS file variants
psc_gwas = read.table("/lustre/scratch115/projects/psc2020/PSC2020_regions.txt",head=T,sep="\t") #here we save this data into a variable we decided to call psc_gwas
#note that sometimes loading data can take time (for instance if our file is 200Mb or more)

gwasvar=read.table(gwasfile,head=F)
refpanvar=read.table(refpan,head=F)
colnames(gwasvar)=c("chr","snp","coor","a0","a1","N","s","MAF","freq_a1","b","se","pvalues")
colnames(refpanvar)=c("chr","pos","rsid","ref","alt")
present_idx=which(gwasvar[,2]%in%refpanvar[,3])
absent_idx=which(!gwasvar[,2]%in%refpanvar[,3])
write.table(gwasvar[absent_idx,],paste(gwasfile,".missing",sep=""),quote=F,col.names=T,row.names=F) #variants not included in fine mapping

refpanplus=merge(refpanvar,gwasvar,by.y="coor",by.x="pos",all.x)
refpan_z=cbind(as.character(refpanplus[,3]),refpanplus[,14]/refpanplus[,15])

#before saving, we need to take into account that the alleles were reversed for some of the entries (and do the same for the Z-scores)

#we want the alleles to match the gwasvar matrix
reverse_alleles=which(!as.character(refpanplus[,4])==as.character(refpanplus[,8]))
refpan_z[reverse_alleles,2]=as.numeric(refpan_z[reverse_alleles,2])*(-1)

#we could perhaps be doing this using the frequencies (for a couple of files I checked, results where almost identical)
#snpstats=read.table(freqfile,head=T)
#freq1=(snpstats$BB*2+snpstats$AB)/(2*(snpstats$AA+snpstats$AB+snpstats$BB))
#dif_a1=abs(freq1-refpanplus$freq_a1)
#dif_a0=abs(freq1-1+refpanplus$freq_a1)
#reverse_alleles=which(dif_a0<dif_a1)
#refpan_z[reverse_alleles,2]=as.numeric(refpan_z[reverse_alleles,2])*(-1)


#ok now write the file to disk
write.table(refpan_z,paste(refpan,".z",sep=""),quote=F,row.names=F,col.names=F)



#END
