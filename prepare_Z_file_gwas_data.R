
## March 15 2017##
## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"
## This is to prepare the .z file for FINEMAP, using the LD info directly from SunGou's PSC binary ped files, located in /lustre/scratch113/projects/pscmeta/psc/finemap/
## The GWAS data it is using are located in : /lustre/scratch115/projects/psc2020/splitfiles/pscgwas/

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("A bim file and an output .z file must be supplied", call.=FALSE)
}
bimfile=args[1] #example: /lustre/scratch113/projects/pscmeta/psc/finemap/IL21.bim
outfile=args[2]


bim=read.table(bimfile,head=F)
chrom=as.character(bim[1,1]) #this assumes all entries are from the same chromosome, which should be the case here.
gwasfile=paste("/lustre/scratch115/projects/psc2020/splitfiles/pscgwas/psc.all_for_reference.chr",chrom,".txt",sep="")
gwasvar=read.table(gwasfile,head=T)
colnames(bim)=c("chr","rsid","nothing","pos","ref","alt")

present_idx=which(gwasvar[,3]%in%bim[,4]) 

refpanplus=merge(bim,gwasvar,by.y="coor",by.x="pos",all.x)
refpan_z=cbind(as.character(refpanplus[,3]),refpanplus[,15]/refpanplus[,16])

#before saving, we need to take into account that the alleles were reversed for some of the entries (and do the same for the Z-scores)
#we want the alleles to match the gwasvar matrix
reverse_alleles=which(!as.character(refpanplus[,5])==as.character(refpanplus[,9])) #NOTE THAT THIS IS NOT A STRAND ISSUE, JUST AN ALLELE RE-ARRANGEMENT ISSUE FREQUENTLY CAUSED BY PLINK
refpan_z[reverse_alleles,2]=as.numeric(refpan_z[reverse_alleles,2])*(-1)

#ok now write the file to disk
write.table(refpan_z,outfile,quote=F,row.names=F,col.names=F)



#END
