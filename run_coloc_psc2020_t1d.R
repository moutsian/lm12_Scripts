
## This is an R file, so these commands should be run in an R environment (just type "R" in your unix terminal)
## In this file we will load the GWAS results for PSC, then load the ichip summary statistics from T1D, and then
## perform colocalisation analysis with coloc 
## Mar 2017

## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("An argument must be supplied for the dataset to be analysed with PSC", call.=FALSE)
}
trait=args[1]
#functional_assay=args[2]

#load coloc library and its dependencies
library(colorspace,lib="/software/team152/Rpackages/")
library(leaps,lib="/software/team152/Rpackages/")
library(robustbase,lib="/software/team152/Rpackages/")
library(inline,lib="/software/team152/Rpackages/")
library(rrcov,lib="/software/team152/Rpackages/")
library(BMA,lib="/software/team152/Rpackages/")
library(coloc,lib="/software/team152/Rpackages/")

case_prop=NULL
if( trait=="t1dcc"){case_prop=0.35}else
if(trait=="t1dmeta"){
case_prop=0.50 #note that I am not sure how to calculate this for family data. If the results are very sensitive to this value we need to decide what exactly to do.
}

##first we need to read in our data from GWAS
## we do this using function read.table

##firstly, let's get in the list of GWAS loci
psc_gwas = read.table("/lustre/scratch115/projects/psc2020/PSC2020_regions.txt",head=T,sep="\t") #here we save this data into a variable we decided to call psc_gwas
#note that sometimes loading data can take time (for instance if our file is 200Mb or more)

#important variables
#snp_of_interest="rs2836883" #in the end, we will do this for all entries but for now let's focus on one SNP


for(j in 1:dim(psc_gwas)[1]){
#for(j in 1:dim(psc_gwas)[1]){snp_of_interest=as.character(psc_gwas[j,1]);print(snp_of_interest)}
snp_of_interest=as.character(psc_gwas[j,1])

#where to save our results:
output_directory=paste("/lustre/scratch115/projects/psc2020/COLOC_RESULTS/",trait,"/",sep="")
output_filename=paste(output_directory,snp_of_interest,".",trait,".coloc.results",sep="")

#now, find the entry (row) where this SNP is on our table
index=which(psc_gwas[,1]==snp_of_interest)

#have a look at this entry:
psc_gwas[index,]
#we get the chromosome from column 2 of the data
chrom=substring(strsplit(as.character(psc_gwas[index,2]),":")[[1]][1],4)

#####################################################################################################################################################
 #now load the data we want from the other GWAS dataset (for the other trait)
psc_func=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/",trait,"/",trait,"_all_summary.Chr_",chrom,".sorted.txt.ready",sep=""),head=T,sep="\t")


#this is big, but we only care about our region around the snp we are interested in
region_start = psc_gwas[index,11] 
region_end=psc_gwas[index,12] #you could also get this using the name of the column, for instance here: psc_gwas$end_position[index]
print(paste0("region: ",region_start,"-",region_end))

#we will use the fact that our data is sorted, to find the starting_position and ending_position
idx_start=which(as.numeric(psc_func[,2])>=region_start)[1]

tmp=which(as.numeric(psc_func[,2])<=region_end)
idx_end=tmp[length(tmp)]
print(paste0("idx_start:",idx_start," idx_end: ",idx_end))

#now that we have that, this is all we care about (from the whole psc_func data)
psc_func_subset = psc_func[idx_start:idx_end,]
#also remove entries where the pvalues are NA
toremove=which(is.na(psc_func_subset[,4]))
toremove_missingMAF=which(is.na(psc_func_subset[,9]))
toremove_all=unique(c(toremove,toremove_missingMAF))
if(length(toremove_all)>0){
psc_func_subset=psc_func_subset[-toremove_all,]
}
#we can now even remove the big data we loaded earlier
rm(psc_func)
######################################################################################################################################################

#now load the data we want from the GWAS data.
psc_gwas_data=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/pscgwas/psc.all_for_reference.chr",chrom,".txt",sep=""),head=T)

#we will use the fact that our data is sorted, to find the starting_position and ending_position
idx_start_gwas=which(as.numeric(psc_gwas_data[,3])>=region_start)[1]
tmp=which(as.numeric(psc_gwas_data[,3])<=region_end)
idx_end_gwas=tmp[length(tmp)]

#now that we have that, this is all we care about (from the whole psc_gwas data)
psc_gwas_data_subset = psc_gwas_data[idx_start_gwas:idx_end_gwas,]

#we can now even remove the big data we loaded earlier
rm(psc_gwas_data)
######################################################################################################################################################

varbeta_gwas=psc_gwas_data_subset[,11]^2


#addition March 1st: Sometimes there is no overlap at all in the set of SNPs (!). In that case just move on.
snp_gwas=psc_gwas_data_subset$snp
snp_func=psc_func_subset[,3]
tmp=which(as.character(snp_func)%in%as.character(snp_gwas))
if(length(tmp)>0){
##this is how SunGou had run coloc:
#results <- coloc.abf(dataset1=list(N=data$N.x, s=data$s.x, pvalues=data$pvalues.x, type="cc"), dataset2=list(N=data$N.y, s=data$s.y, pvalues=data$pvalues.y, type="cc"), MAF=data$MAF, p12=1e-06) # p12 is set as what Fortune et al. 2015 recommended

#this is the one which ran fine:
results = coloc.abf(dataset1=list(snp=as.character(psc_gwas_data_subset$snp),beta=psc_gwas_data_subset[,10],varbeta=varbeta_gwas,type="cc",s=0.25,MAF=psc_gwas_data_subset$MAF),
 dataset2=list(snp=as.character(psc_func_subset[,3]),pvalues=as.numeric(psc_func_subset[,4]),N=psc_func_subset$SampleSize, MAF=psc_func_subset$Freq_MinorAllele, type="cc",s=case_prop)
,p12=1e-06)

}
coloc_results=matrix(ncol=6,nrow=1,0)
coloc_results[1,]=results$summary
colnames(coloc_results)=c("n_snps","H0","H1","H2","H3","H4")
write.table(coloc_results,output_filename,quote=F,col.names=T,row.names=F,sep="\t")
#now we know from them the README file from Louella that the beta is column 5 and the SE column 9.

}#here we close the for loop for the snp_of_interest

#########################################



#END
