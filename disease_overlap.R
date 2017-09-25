
printtopdf=T;

if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/ukb4940.tab",sep="\t",head=T)}
disease_coding=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/disease_coding.csv",sep=",",head=T)
Disease="fibromyalgia" #examples: "crohns disease" "irritable bowel syndrome" "sclerosing cholangitis" "malabsorption/coeliac disease" "chronic fatigue syndrome"
primafield=disease_coding[which(disease_coding[,2]==Disease),1]
if(length(primafield)==0){print(paste("Disease ",Disease,"could not be found"))}

primary_entries=NULL
for(i in 11:39){
primary_entries=c(primary_entries,which(disease_info[,i]==primafield))
}
disease_filtered=disease_info[primary_entries,11:39]
healthy_individuals=sum(is.na(disease_info[,11]))
n_diseases=table(rowSums(!is.na(disease_filtered)))

disease_counts=table(unlist(disease_filtered))
sorted_counts=sort(disease_counts,decreasing=T)

#plot
sorted_diseases=matrix(ncol=2,nrow=10,NA)
for(i in 2:11){
idx=which(disease_coding[,1]==as.numeric(names(sorted_counts[i])))
sorted_diseases[i-1,]=c(disease_coding[idx,1],as.character(disease_coding[idx,2]))
}

tmp=gsub(' ','_',Disease)
tmp=gsub('/','_',tmp)

outfile=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Overlap_of_",tmp,"_w_other_diseases_excl_cancer.pdf",sep="")
if(printtopdf){pdf(outfile,width=10)}
par(mar=c(20,4,4,2))
barplot(sorted_counts[2:11]/sorted_counts[1], main=list(paste("Overlap of ",Disease," with other diseases (excl. cancer)",sep=""),cex=1.1),col="darkblue", ylab=list("Percentage (%) of Overlap",cex=1.2),
names.arg=sorted_diseases[,2],las=2,cex.names=1.2)
mtext("Results from Initial Assessment (Field: 20002)",3)
if(printtopdf){dev.off()}

outfile_ndis=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Total_number_of_diseases_excl_cancer_for_people_with_",tmp,".pdf",sep="")
if(printtopdf){pdf(outfile_ndis,width=10)}
par(mar=c(10,10,10,2))
n11plus=sum(n_diseases[11:length(n_diseases)])
barplot(c(n_diseases[1:10]/dim(disease_filtered)[1],n11plus/dim(disease_filtered)[1]),main=list(paste("Number of diseases (excl. cancer) people with\n ",Disease," have",sep=""),cex=1.1),
ylab=list(paste("Percentage (%) of all individuals with\n",Disease,sep=""),cex=1.1),names.arg=c(1:10,"11+"), xlab=list("Total number of diseases",cex=1.1))
mtext("Results from Initial Assessment (Field: 20002)",3)
if(printtopdf){dev.off()}



prevalence_in_uk_biobank=matrix(ncol=1,nrow=10,0)
for(i in 11:39){
for(j in 1:10){
prevalence_in_uk_biobank[j]=prevalence_in_uk_biobank[j]+length(which(disease_info[,i]==sorted_diseases[j,1]))
}}
prevalence_in_uk_biobank_pt=prevalence_in_uk_biobank/dim(disease_info)[1]
outfile_comp=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Overlap_of_",tmp,"_w_other_diseases_excl_cancer_compared_to_biobank.pdf",sep="")
if(printtopdf){pdf(outfile_comp,width=10)}
par(mar=c(20,4,4,2))
TOPLOT=cbind(sorted_counts[2:11]/sorted_counts[1],prevalence_in_uk_biobank_pt)
barplot(t(TOPLOT),beside=T,space=c(.25,.75), main=list(paste("Overlap of ",Disease," with other diseases (excl. cancer)",sep=""),cex=1.1),col=c("darkblue","red"), ylab=list("Percentage (%) of Overlap",cex=1.1),
names.arg=sorted_diseases[,2],las=2,cex.names=1.2)
mtext("Results from Initial Assessment (Field: 20002)",3)
legend("topright",pch=19,col=c("darkblue","red"),c(paste("Overlap with ",Disease,sep=""),"Prevalence in UK Biobank"))
if(printtopdf){dev.off()}




#END