
printtopdf=T;
MIN_COUNTS=10
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

#the following are individuals which have IBS according to the hospital records but have not been self-reported with it:
idv_disc=record_matching[(which(record_matching[,2]=="K589" & record_matching[,4]!="1" & record_matching[,5]!="1")),1]
idv_idx=matrix(ncol=1,nrow=length(idv_disc),0)
for(i in 1:length(idv_disc)){
print(i);
idv_idx[i]=which(disease_info[,1]==idv_disc[i])
}
idv_disc_traits=disease_info[idv_idx,11:39]


disease_filtered=disease_info[secondary_entries[newcases],11:39]
disease_filtered_secondary=disease_info[secondary_entries[newcases],40:68]

n_diseases=table(rowSums(!is.na(disease_filtered))) #this vector contains the number of individuals with at least x diseases (e.g. n_diseases[4] is the number of ppl with 4 diseases)

disease_counts=table(unlist(disease_filtered))
sorted_counts=sort(disease_counts,decreasing=T)

disease_counts_secondary=table(unlist(disease_filtered_secondary))
disidx=names(sorted_counts[sorted_counts>=MIN_COUNTS])
counts_sec=matrix(ncol=1,nrow=length(disidx),0)
for(i in 1:length(disidx)){
tmpidx=which(names(disease_counts_secondary)==disidx[i])
counts_sec[i]=disease_counts_secondary[tmpidx]
#names(counts_sec)=disease_coding[disease_coding[,1]==disidx[i],2]
}
names(counts_sec)=disidx

sorted_diseases=matrix(ncol=2,nrow=dim(sorted_counts)[1],NA)
for(i in 1:dim(sorted_diseases)[1]){
idx=which(disease_coding[,1]==as.numeric(names(sorted_counts[i])))
sorted_diseases[i,]=c(disease_coding[idx,1],as.character(disease_coding[idx,2]))
}

tmp=gsub(' ','_',Disease)
tmp=gsub('/','_',tmp)
tmp=gsub('\\(','_',tmp)
tmp=gsub('\\)','_',tmp)


outfile_contrast=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Init_vs_Secondary_Assessment_min",MIN_COUNTS,"_",tmp,".pdf",sep="")
if(printtopdf){pdf(outfile_contrast)}
par(mar=c(15,4,4,2))
BOTH1=cbind(sorted_counts[sorted_counts>=MIN_COUNTS],counts_sec)
barplot(BOTH1[,1], main=list(paste("Secondary conditions for individuals\nreporting ",tmp," ONLY in the secondary assessment",sep=""),cex=1.1),col="darkblue",
ylab=list("Counts",cex=1),names.arg=sorted_diseases[sorted_counts>=MIN_COUNTS,2],las=2,cex.names=0.9,ylim=c(0,max(BOTH1)),)
barplot(t(BOTH1[,2]),xaxt='n',yaxt='n',density=20,col="darkred",add=T)
legend("topright",pch=19,col=c("darkblue","red"),c("Initial Assessment","Secondary Assessment"),cex=0.8)
mtext(paste("(Min of ",MIN_COUNTS," individuals with the same condition in the initial assessment.)",sep=""),3,cex=0.7)
if(printtopdf){dev.off()}

#check if we have people reported as having the disease in the initial but not in the secondary assessment.
n_change=0;
n_nochange=0;
n_healthynow=0;
change_idx=NULL
nochange_idx=NULL
for(i in 1:length(primary_entries)){
	if(!is.na(disease_info[primary_entries[i],6])){ 
			if(!primafield%in%disease_info[primary_entries[i],40:68]){
			n_change=n_change+1;
			change_idx=c(change_idx,primary_entries[i])
				if(is.na(disease_info[primary_entries[i],40])){
				n_healthynow=n_healthynow+1;
				}
			}
			else{
			n_nochange=n_nochange+1;
			nochange_idx=c(nochange_idx,primary_entries[i])
			}
	}
}




people_in_both_assessments=primary_entries[which(!is.na(disease_info[primary_entries,6]))]

outfile_cons=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Init_vs_Secondary_Assessment_consistency_",tmp,".pdf",sep="")
if(printtopdf){pdf(outfile_cons,width=12)}
barplot(n_nochange/(length(people_in_both_assessments)),col="darkblue",ylim=c(0,1),main=paste("Percentage of people in both assessments who keep having ",tmp,sep=""),ylab="Percentage",xlim=c(-1,2.5))
abline(h=c(0,0.2,0.4,0.6,0.8,1),col="gray",lwd=2,lty=2)
mtext(paste("People with ",tmp," in the init assessment which were included in the secondary: ",length(people_in_both_assessments),sep=""),3,cex=0.8)
if(printtopdf){dev.off()}


#now have a look at what are the changes in the secondary conditions for these individuals who keep having inflammatory bowel disease
newdiseases=sort(table(unlist(disease_info[nochange_idx,40:68])),decreasing=T)
olddiseases=sort(table(unlist(disease_info[nochange_idx,11:39])),decreasing=T)
toplot=unique(sort(c(names(which(olddiseases>=MIN_COUNTS)),names(which(newdiseases>=MIN_COUNTS)))))
BOTH2=matrix(ncol=2,nrow=length(toplot),0)
for(i in 1:length(toplot)){
BOTH2[i,1]=olddiseases[which(names(olddiseases)==toplot[i])]
BOTH2[i,2]=newdiseases[which(names(newdiseases)==toplot[i])]
}
colnames(BOTH2)=c("Init","Secondary")
rownames(BOTH2)=toplot

BOTH2=BOTH2[order(BOTH2[,1],decreasing=T),]
disease_name=NULL;
for(i in 1:length(toplot)){
disease_name=c(disease_name,as.character(disease_coding[which(disease_coding[,1]==rownames(BOTH2)[i]),2]))}

outfile_contrast3=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Init_vs_Secondary_Assessment_min",MIN_COUNTS,"_ppl_with_",tmp,"_in_both.pdf",sep="")
if(printtopdf){pdf(outfile_contrast3)}
par(mar=c(15,4,4,2))
barplot(BOTH2[-1,1],col="darkblue",names.arg=disease_name[-1],ylim=c(0,max(BOTH2[-1,])),main=list(paste("Secondary assessment condition for individuals\nreporting ",tmp," in BOTH assessments",sep=""),cex=1.1),
las=2,cex.names=0.9,ylab=list("Counts",cex=1))
barplot(BOTH2[-1,2],col="red",border="gray",add=T,xaxt='n',yaxt='n',density=20)
mtext(paste("(Min of ",MIN_COUNTS," individuals with the same condition in the secondary assessment.)",sep=""),3,cex=0.7)
legend("topright",pch=19,col=c("darkblue","red"),c("Initial Assessment","Secondary Assessment"),cex=0.8)
if(printtopdf){dev.off()}

newdiseases=sort(table(unlist(disease_info[change_idx,40:68])),decreasing=T)
diseases_toplot=which(newdiseases>=MIN_COUNTS)
disease_id=names(diseases_toplot)
disease_name=NULL;
#also get the counts for the individuals with the same secondary condition in the initial assessment
counts_init=matrix(ncol=1,nrow=length(diseases_toplot),0)
for(i in 1:length(diseases_toplot)){
counts_init[i]=sorted_counts[which(names(sorted_counts)==names(newdiseases)[diseases_toplot[i]])]
}
for(i in 1:length(disease_id)){
disease_name=c(disease_name,as.character(disease_coding[which(disease_coding[,1]==disease_id[i]),2]))}


BOTH=cbind(c(counts_init,0),c(newdiseases[diseases_toplot],n_healthynow))
BOTH=cbind(BOTH,BOTH[,2]-BOTH[,1])
colnames(BOTH)=c("counts_init","counts_sec","diff")

outfile_contrast2=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Init_vs_Secondary_Assessment_min",MIN_COUNTS,"_ppl_initially_diagnosed_with_",tmp,".pdf",sep="")
if(printtopdf){pdf(outfile_contrast2)}
par(mar=c(15,4,4,2))
barplot(BOTH[,1],col="darkblue",names.arg=c(disease_name,"healthy"),ylim=c(0,max(BOTH)),main=list(paste("Secondary assessment condition for individuals\nreporting ",tmp," ONLY in the initial assessment",sep=""),cex=1.1),
las=2,cex.names=0.9,ylab=list("Counts",cex=1))
barplot(BOTH[,2],col="red",names.arg=c(disease_name,"healthy"),border="gray",add=T,xaxt='n',yaxt='n',density=20)
mtext(paste("(Min of ",MIN_COUNTS," individuals with the same condition in the secondary assessment.)",sep=""),3,cex=0.7)
legend("topright",pch=19,col=c("darkblue","red"),c("Initial Assessment","Secondary Assessment"),cex=0.8)
if(printtopdf){dev.off()}




#END