
printtopdf=T;
MIN_COUNTS=3
#if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/ibs/data/ukb7725.tab",sep="\t",head=T)}
if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/ibs/data/phenotypes_for_genotyped_samples_only.txt",sep="\t",head=T)} #this only checks the genotyped samples
disease_coding=read.table("/lustre/scratch115/teams/anderson/ukbiobank/project_autoimmunity/uncompressed_files/disease_coding.csv",sep=",",head=T)
Disease="irritable bowel syndrome" #examples: "crohns disease" "irritable bowel syndrome" "sclerosing cholangitis" "malabsorption/coeliac disease"
primafield=disease_coding[which(disease_coding[,2]==Disease),1]
if(length(primafield)==0){print(paste("Disease ",Disease,"could not be found"))}
bridging=read.table("/lustre/scratch115/teams/anderson/ukbiobank/ibs/data/bridging_file.csv",sep=",",head=F) #bridging file for our sample IDs with the central sanger ones
colnames(bridging)=c("ID_ours","ID_sanger")

ioanna_sample_file=read.table("/lustre/scratch115/projects/ukbiobank/it3/data/Ioanna.sample2",head=T)
european_qcpass_kostas=read.table("/lustre/scratch115/projects/ukbiobank/lm12/data/Allchr_dir-typ_sex_97miss_BLVfails_EUR_Het3sd_Relatedness_SNPQC.fam",head=F) #this is from Kostas via Fernando -the samples to keep
PCA=read.table("/lustre/scratch115/projects/ukbiobank/it3/data/Ioanna.cleanphenofileforBOLTLMM_fullwithNAs",head=T)

#The PCA and ioanna_sample_file files are in the same order, and this is the order we need to maintain - specifically, this of Ioanna I presume. 
#disease_info_reordering_idx=matrix(ncol=1,nrow=dim(PCA)[1],NA)
for(i in 1:length(disease_info_reordering_idx)){
print(paste(i,":",sep=""))
bridged=which(bridging[,2]==PCA[i,1])
if(length(bridged)>0){
tmp=which(disease_info[,1]==bridging[bridged,1])
if(length(tmp>0)){
disease_info_reordering_idx[i]=tmp}else{
disease_info_reordering_idx[i]=NA}
}
} 
DIRI=disease_info_reordering_idx
DIRI=DIRI[-length(DIRI)]


idx=which(bridging[,1]%in%disease_info[primary_entries,1])
idx_in_sanger_files=which(PCA[,1]%in%bridging[idx,2])

#these entries are not in the PCA /Ioanna files:
not_in=which(!bridging[idx,2]%in%PCA[,1])
samples_not_in_ioanna=bridging[idx[not_in],]


disease_info_reordered=disease_info[order(DIRI),]
disinfo=cbind(PCA[,1],disease_info_reordered) #we should now be able to use this file for all phenotype extractions without worrying about sample re-ordering


primary_entries=NULL
for(i in 54:140){
primary_entries=c(primary_entries,which(disinfo[,i]==primafield))
}
primary_entries=unique(primary_entries)

ibs_selfrep=matrix(ncol=1,nrow=dim(disinfo)[1],0)
ibs_selfrep[primary_entries]=1
ibs_selfrep[which(is.na(disinfo[,2]))]=NA

#hospital_records=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/Coeliac_malabsorption.status_from_hospital_admission.csv",head=T,sep=",")


hosp_entries_prim=NULL
code1="K580"
code2="K589"
for(i in 673:1049){ #note that the index for disinfo is +2 compared to file:///C:/Academic/SANGER/UK%20Biobank/UK%20Biobank%20-%20IBS%20-%201767/Data/ukb7725.html, because the latter is zero-indexed and 
#also we have added an additional ID column.
hosp_entries_prim=c(hosp_entries_prim,which(disinfo[,i]==code1),which(disinfo[,i]==code2))
}
hosp_entries_prim=unique(hosp_entries_prim)

ibs_hosp_prim=matrix(ncol=1,nrow=dim(disinfo)[1],0)
ibs_hosp_prim[hosp_entries_prim]=1
ibs_hosp_prim[which(is.na(disinfo[,2]))]=NA


hosp_entries_sec=NULL
code1="K580"
code2="K589"
for(i in 1078:1510){ #note that the index for disinfo is +2 compared to file:///C:/Academic/SANGER/UK%20Biobank/UK%20Biobank%20-%20IBS%20-%201767/Data/ukb7725.html, because the latter is zero-indexed and 
#also we have added an additional ID column.
hosp_entries_sec=c(hosp_entries_sec,which(disinfo[,i]==code1),which(disinfo[,i]==code2))
}
hosp_entries_sec=unique(hosp_entries_sec)

ibs_hosp_sec=matrix(ncol=1,nrow=dim(disinfo)[1],0)
ibs_hosp_sec[hosp_entries_sec]=1
ibs_hosp_sec[which(is.na(disinfo[,2]))]=NA

ibs_hosp_all=ibs_hosp_prim+ibs_hosp_sec
ibs_hosp_all[ibs_hosp_all>1]=1

ibs_any=ibs_selfrep+ibs_hosp_all
ibs_any[ibs_any>1]=1
age=2016-disinfo[,4]+9/12-disinfo[,5]/12 #9/12 for Sept

#the following are individuals which have IBS according to the hospital records - primary assessment, but have not been self-reported with it:
idv_disc=record_matching[(which(record_matching[,3]=="41202" & record_matching[,4]!="1" & record_matching[,5]!="1")),1] 
idv_idx=matrix(ncol=1,nrow=length(idv_disc),0)
for(i in 1:length(idv_disc)){
print(i);
idv_idx[i]=which(disease_info[,1]==idv_disc[i])
}
idv_disc_traits=disease_info[idv_idx,11:39]

chip=matrix(ncol=1,nrow=dim(disinfo)[1],0)
chip[disinfo[,539]<0]=1 #BILEVE
chip[disinfo[,539]>0]=2 #AXIOM

#prepare pheno file for output
phenoout=cbind(PCA[,c(1:4)],ibs_selfrep,ibs_hosp_all,ibs_any,chip,age,PCA[,c(6:15)])
colnames(phenoout)[1:3]=c("ID_1","ID_2","missing")
phenofin=rbind(secondline,phenoout)
write.table(phenofin,"/lustre/scratch115/projects/ukbiobank/lm12/data/ibs_pheno.sample",quote=F,row.names=F,col.names=T,sep=" ")

#prepare exclusion list for SNPTEST2 based on Kostas's QC
tmp=which(disinfo[,1]=="3078126")
toremove_idx=c(which(!disinfo[,1]%in%european_qcpass_kostas[,1]),tmp)
write.table(disinfo[toremove_idx,1],"/lustre/scratch115/projects/ukbiobank/lm12/data/ibs_samples_to_exclude.txt",quote=F,row.names=F,col.names=F)




#END