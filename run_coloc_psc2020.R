## This is an R file, so these commands should be run in an R environment (just type "R" in your unix terminal)
## In this file we will load the GWAS results for PSC, then load the functional analyses results from Blueprint and then
## perform colocalisation analysis with coloc (at first)
## Feb 2017

## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"

#load coloc library and its dependencies
library(colorspace,lib="/software/team152/Rpackages/")
library(leaps,lib="/software/team152/Rpackages/")
library(robustbase,lib="/software/team152/Rpackages/")
library(inline,lib="/software/team152/Rpackages/")
library(rrcov,lib="/software/team152/Rpackages/")
library(BMA,lib="/software/team152/Rpackages/")
library(coloc,lib="/software/team152/Rpackages/")

##first we need to read in our data from GWAS
## we do this using function read.table

##firstly, let's get in the list of GWAS loci
psc_gwas = read.table("/lustre/scratch115/projects/psc2020/PSC2020_regions.txt",head=T,sep="\t") #here we save this data into a variable we decided to call psc_gwas
#note that sometimes loading data can take time (for instance if our file is 200Mb or more)


#important variables
snp_of_interest="rs2836883" #in the end, we will do this for all entries but for now let's focus on one SNP
celltype_of_interest="mono" #this could also be "tcel" or "mono"
functional_assay="gene_nor_combat" #this could be "K4ME1_log2rpm" for instance


#depending on the celltype and assay, we want to get the number of counts correct:
func_counts=NULL
if(celltype_of_interest=="mono"){
	if(functional_assay=="gene_nor_combat"){
		func_counts=194
	}else if(functional_assay=="meth_M"){
		func_counts=196
	}else if(functional_assay=="K4ME1_log2rpm"){
		func_counts=172
	}else{
		func_counts=162
	}
}else if(celltype_of_interest=="neut"){
	if(functional_assay=="gene_nor_combat"){
		func_counts=192
	}else if(functional_assay=="meth_M"){
		func_counts=197
	}else if(functional_assay=="K4ME1_log2rpm"){
		func_counts=173
	}else{
		func_counts=174
	}
}else{ celltype_of_interest=="tcel"
	if(functional_assay=="gene_nor_combat"){
		func_counts=171
	}else if(functional_assay=="meth_M"){
		func_counts=133
	}else if(functional_assay=="K4ME1_log2rpm"){
		func_counts=104
	}else{
		func_counts=142
	}	
}

#some ways of checking the data we just loaded
head(psc_gwas) # this can give you an idea of what is in there, just as we did "head" in unix
dim(psc_gwas) #this will give you the dimensions (rows x columns) of the table
psc_gwas[,1] # this will give you the contents of the first column of the table
psc_gwas[1,] # this will give you the contents of the first row of the table


#now, find the entry (row) where this SNP is on our table
index=which(psc_gwas[,1]==snp_of_interest)

#have a look at this entry:
psc_gwas[index,]
#we get the chromosome from column 2 of the data
chrom=substring(strsplit(as.character(psc_gwas[index,2]),":")[[1]][1],4)

#where to save our results:
output_directory=paste("/lustre/scratch115/projects/psc2020/COLOC_RESULTS/",celltype_of_interest,"/",functional_assay,"/",sep="")
output_filename=paste(output_directory,snp_of_interest,".",celltype_of_interest,".",functional_assay,".coloc.results",sep="")



 #####################################################################################################################################################
 #now load the data we want from the functional dataset
psc_func=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/",celltype_of_interest,"/",celltype_of_interest,"_",functional_assay,"_peer_10_all_summary.Chr_",chrom,".sorted.txt",sep=""),head=F)

#this is big, but we only care about our region around the snp we are interested in
region_start = psc_gwas[index,11] 
region_end=psc_gwas[index,12] #you could also get this using the name of the column, for instance here: psc_gwas$end_position[index]

#we will use the fact that our data is sorted, to find the starting_position and ending_position
idx_start=which(as.numeric(psc_func[,2])>=region_start)[1]
tmp=which(as.numeric(psc_func[,2])<=region_end)
idx_end=tmp[length(tmp)]

#now that we have that, this is all we care about (from the whole psc_func data)
psc_func_subset = psc_func[idx_start:idx_end,]

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


#we need to split the functional file by transcript
all_transcripts=unique(as.character(psc_func_subset[,5]))
coloc_results=matrix(nrow=length(all_transcripts),ncol=7,NA)
coloc_results[,1]=all_transcripts
##iterate by transcript, and run coloc for each of them.
for(i in 1:length(all_transcripts)){
transcript_index=which(as.character(psc_func_subset[,5])==all_transcripts[i])
##prepare varbetas:
varbeta_func=psc_func_subset[transcript_index,11]^2


##this is how SunGou had run coloc:
#results <- coloc.abf(dataset1=list(N=data$N.x, s=data$s.x, pvalues=data$pvalues.x, type="cc"), dataset2=list(N=data$N.y, s=data$s.y, pvalues=data$pvalues.y, type="cc"), MAF=data$MAF, p12=1e-06) # p12 is set as what Fortune et al. 2015 recommended
#this is the one which ran fine:
results = coloc.abf(dataset1=list(snp=psc_gwas_data_subset$snp,beta=psc_gwas_data_subset[,10], s=0.25, N=psc_gwas_data_subset$N, varbeta=varbeta_gwas,MAF=psc_gwas_data_subset$MAF,type="cc"), dataset2=list(snp=psc_func_subset[transcript_index,4],beta=psc_func_subset[transcript_index,7],N=func_counts, MAF=psc_func_subset[transcript_index,10],varbeta=varbeta_func, type="quant"),p12=1e-06)

coloc_results[i,2:7]=results$summary
}
colnames(coloc_results)=c("Gene/Feature","n_snps","H0","H1","H2","H3","H4")
write.table(coloc_results,output_filename,quote=F,col.names=T,row.names=F,sep="\t")
#now we know from them the README file from Louella that the beta is column 5 and the SE column 9.

#########################################
#########################################
# So what is remaining? To get a paper first and a nobel prize down the line.



#END