
printtopdf=T;
MIN_COUNTS=3
if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/ukb4940.tab",sep="\t",head=T)}
disease_coding=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/disease_coding.csv",sep=",",head=T)
Disease="irritable bowel syndrome" #examples: "crohns disease" "irritable bowel syndrome" "sclerosing cholangitis" "malabsorption/coeliac disease"
primafield=disease_coding[which(disease_coding[,2]==Disease),1]
if(length(primafield)==0){print(paste("Disease ",Disease,"could not be found"))}

hospital_records=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/IBS.status_from_hospital_admission.csv",head=T,sep=",")

primary_entries=NULL
secondary_entries=NULL
for(i in 11:39){
primary_entries=c(primary_entries,which(disease_info[,i]==primafield))
}
for(i in 40:68){
secondary_entries=c(secondary_entries,which(disease_info[,i]==primafield))
}
newcases=which(!secondary_entries%in%primary_entries)
oldcases=which(secondary_entries%in%primary_entries)

# NOTE that the ranking in the disease_info and the hospital_records datasets are not the same.
idx=which(hospital_records[,2]!=0)
record_matching=matrix(ncol=5,nrow=length(idx),0)
record_matching[,1:3]=as.matrix(hospital_records[idx,])
for(i in 1:length(idx)){
	if(which(disease_info[,1]==hospital_records[idx[i],1])%in%primary_entries){
		record_matching[i,4]=1
	}
	if(which(disease_info[,1]==hospital_records[idx[i],1])%in%secondary_entries){
		record_matching[i,5]=1
	}
}




#the following are individuals which have IBS according to the hospital records - primary assessment, but have not been self-reported with it:
idv_disc=record_matching[(which(record_matching[,3]=="41202" & record_matching[,4]!="1" & record_matching[,5]!="1")),1] #should be 741 individuals
idv_idx=matrix(ncol=1,nrow=length(idv_disc),0)
for(i in 1:length(idv_disc)){
print(i);
idv_idx[i]=which(disease_info[,1]==idv_disc[i])
}
idv_disc_traits=disease_info[idv_idx,11:39]


#we should check the dates, in case most individuals in the above category were entered in the hospital AFTER completing the questionaire (and hence may have not had IBS yet at the time
# of completing the questionaire) 
# The date of assessment is field 53 (column 5 in the disease_info matrix). For the hospital episodes it may be field 41096, but I think (as of jan 2016), that the data is only available in situ,
# so we cannot check that for now.

# In the meantime, I am just checking what are the diseases that these 741 individuals are self-reporting.
counts_idv_disc_traits=table(unlist(idv_disc_traits))
sorted_counts_disc=sort(counts_idv_disc_traits,decreasing=T)

sorted_diseases=matrix(ncol=2,nrow=dim(sorted_counts_disc)[1],NA)
for(i in 1:dim(sorted_diseases)[1]){
idx=which(disease_coding[,1]==as.numeric(names(sorted_counts_disc[i])))
sorted_diseases[i,]=c(disease_coding[idx,1],as.character(disease_coding[idx,2]))
}
prevalence_in_uk_biobank=matrix(ncol=1,nrow=dim(sorted_diseases)[1],0)
for(i in 11:29){
for(j in 1:dim(sorted_diseases)[1]){
prevalence_in_uk_biobank[j]=prevalence_in_uk_biobank[j]+length(which(disease_info[,i]==sorted_diseases[j,1]))
}}
n_individuals=dim(disease_info)[1]
prevalence_in_uk_biobank_pt=prevalence_in_uk_biobank/n_individuals 
prevalence_in_discordant_individuals_pt=sorted_counts_disc/length(idv_disc)

prevalence_comparison=cbind(sorted_diseases,prevalence_in_uk_biobank_pt,prevalence_in_discordant_individuals_pt, sorted_counts_disc)
prevalence_comparison=cbind(prevalence_comparison,as.numeric(prevalence_comparison[,4])/as.numeric(prevalence_comparison[,3]))
colnames(prevalence_comparison)=c("Code","Disease","Prev_in_UKB","Prev_in_discordant_ind","Counts_in_disc_ind","ratio")
#from the prevalence comparison table, we can get the biggest differences in prevalence and make two plots (less/more prevalent)

idx=which(as.numeric(prevalence_comparison[,5])>=MIN_COUNTS)
prevcomp_filtered=prevalence_comparison[idx,]
prevcomp_filt_high=prevcomp_filtered[order(as.numeric(prevcomp_filtered[,6]),decreasing=T),]
prevcomp_filt_low=prevcomp_filtered[order(as.numeric(prevcomp_filtered[,6]),decreasing=F),]

tmp=gsub(' ','_',Disease)
tmp=gsub('/','_',tmp)
tmp=gsub('\\(','_',tmp)
tmp=gsub('\\)','_',tmp)


outfile_contrast=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Self_reported_traits_for_people_with_medical_record_of_IBS.",MIN_COUNTS,"_",tmp,".pdf",sep="")
if(printtopdf){pdf(outfile_contrast)}
par(mar=c(15,4,4,2))
barplot(as.numeric(prevcomp_filt_high[1:10,4]), main=list(paste("Self-reported conditions for individuals\nwith a medical record of ",Disease," but no self-report",sep=""),cex=1.1),col="darkblue",
ylab=list("Counts",cex=1),names.arg=prevcomp_filt_high[1:10,2],las=2,cex.names=0.9,ylim=c(0,max(as.numeric(prevcomp_filt_high[1:10,4]))),)
barplot(as.numeric(prevcomp_filt_high[1:10,3]),xaxt='n',yaxt='n',density=20,col="darkred",add=T)
legend("topright",pch=19,col=c("darkblue","red"),c("Prevalence in the 741 individuals","Prevalence in UKB"),cex=0.8)
mtext(paste("(Min of ",MIN_COUNTS," in the 741 individuals.)",sep=""),3,cex=0.7)
if(printtopdf){dev.off()}

outfile_contrast=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Self_reported_traits_for_people_with_medical_record_of_IBS.",MIN_COUNTS,"_",tmp,".less_prev.pdf",sep="")
if(printtopdf){pdf(outfile_contrast)}
par(mar=c(15,4,4,2))
barplot(as.numeric(prevcomp_filt_low[1:10,4]), main=list(paste("Self-reported conditions for individuals\nwith a medical record of ",Disease," but no self-report",sep=""),cex=1.1),col="darkblue",
ylab=list("Counts",cex=1),names.arg=prevcomp_filt_low[1:10,2],las=2,cex.names=0.9,ylim=c(0,max(as.numeric(prevcomp_filt_low[1:10,3]))),)
barplot(as.numeric(prevcomp_filt_low[1:10,3]),xaxt='n',yaxt='n',density=20,col="darkred",add=T)
legend("topleft",pch=19,col=c("darkblue","red"),c("Prevalence in the 741 individuals","Prevalence in UKB"),cex=0.8)
mtext(paste("(Min of ",MIN_COUNTS," in the 741 individuals.)",sep=""),3,cex=0.7)
if(printtopdf){dev.off()}



#END