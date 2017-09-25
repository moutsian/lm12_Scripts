
YEAR_WINDOW=-1 #in years - but it accepts decimals (e.g. 1.5) . Leave as -1 if no limit
printtopdf=T;
if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/ukb4940.tab",sep="\t",head=T)}
disease_coding=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/disease_coding.csv",sep=",",head=T)
Disease="type 1 diabetes" #examples: "crohns disease" "irritable bowel syndrome" "sclerosing cholangitis" "malabsorption/coeliac disease" "infectious mononucleosis / glandular fever / epstein barr virus (ebv)"
primafield=disease_coding[which(disease_coding[,2]==Disease),1]
if(length(primafield)==0){print(paste("Disease ",Disease,"could not be found"))}


tmp=gsub(' ','_',Disease)
tmp=gsub('/','_',tmp)
tmp=gsub('\\(','_',tmp)
tmp=gsub('\\)','_',tmp)

primary_entries=NULL
primary_entries_aao=NULL
secondary_entries=NULL
secondary_entries_aao=NULL

for(i in 11:39){
primary_entries=c(primary_entries,which(disease_info[,i]==primafield))
primary_entries_aao=c(primary_entries_aao,disease_info[which(disease_info[,i]==primafield),i+58])
}
for(i in 40:68){
secondary_entries=c(secondary_entries,which(disease_info[,i]==primafield))
secondary_entries_aao=c(secondary_entries_aao,disease_info[which(disease_info[,i]==primafield),i+58])
}

primary_entries_aao[primary_entries_aao<0]=NA

outfile_AAO=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","AAO_",tmp,".pdf",sep="")
if(printtopdf){pdf(outfile_AAO,width=10)}
par(mar=c(4,4,4,2))
par(mfrow=c(1,2))
hist(primary_entries_aao,col="#A2A475",main=list(paste("AAO for individuals with ",tmp,"\nin initial assessment"),cex=1),xlab=list("AAO",cex=1))
primean=mean(primary_entries_aao,na.rm=T)
abline(v=primean,col="darkred",lty=2,lwd=2)
mtext(paste("Mean age:",primean,". Median age: ",median(primary_entries_aao,na.rm=T),sep=""),3,cex=0.7)
hist(secondary_entries_aao,col="#A2A475",main=list(paste("AAO for individuals with ",tmp,"\nin secondary assessment"),cex=1),xlab=list("AAO",cex=1))
secmean=mean(secondary_entries_aao,na.rm=T)
abline(v=secmean,col="darkred",lty=2,lwd=2)
mtext(paste("Mean age:",secmean,". Median age: ",median(secondary_entries_aao,na.rm=T),sep=""),3,cex=0.7)
if(printtopdf){dev.off()}

#now check what diseases did people have BEFORE getting the disease
disease_filtered=disease_info[primary_entries,11:39]
previous_diseases=c(NA,NA,NA,NA)
later_diseases=c(NA,NA,NA,NA)

yearwindow=""
yearwindowtitle=""
if(YEAR_WINDOW!=-1){
yearwindow=paste(YEAR_WINDOW,"yr",sep="")
yearwindowtitle=paste("(within ",YEAR_WINDOW,"year)",sep="")

	for(i in 1:length(primary_entries)){
	idx=which(disease_info[primary_entries[i],11:39]==primafield)
	aao=disease_info[primary_entries[i],11+idx-1+58]
	idx2=which(!is.na(disease_info[primary_entries[i],11:39]))
	idx2=idx2[-which(idx2==idx)]
		if(length(idx2)>0){
			for(j in 1:length(idx2)){
				if(disease_info[primary_entries[i],11+idx2[j]-1+58]<aao & disease_info[primary_entries[i],11+idx2[j]-1+58]>(aao-YEAR_WINDOW)){
				previous_diseases=rbind(previous_diseases,c(primary_entries[i],disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),1],as.character(disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),2]),disease_info[primary_entries[i],11+idx2[j]-1+58]))
				}else if(disease_info[primary_entries[i],11+idx2[j]-1+58]>aao & disease_info[primary_entries[i],11+idx2[j]-1+58]<aao+YEAR_WINDOW){
				later_diseases=rbind(later_diseases,c(primary_entries[i],disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),1],as.character(disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),2]),disease_info[primary_entries[i],11+idx2[j]-1+58]))
				}
			}
		}
	}
}else{
	for(i in 1:length(primary_entries)){
		idx=which(disease_info[primary_entries[i],11:39]==primafield)
		aao=disease_info[primary_entries[i],11+idx-1+58]
		idx2=which(!is.na(disease_info[primary_entries[i],11:39]))
		idx2=idx2[-which(idx2==idx)]
			if(length(idx2)>0){
				for(j in 1:length(idx2)){
					if(disease_info[primary_entries[i],11+idx2[j]-1+58]<aao){
					previous_diseases=rbind(previous_diseases,c(NA,disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),1],as.character(disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),2]),disease_info[primary_entries[i],11+idx2[j]-1+58]))
					}else if(disease_info[primary_entries[i],11+idx2[j]-1+58]>aao){
					later_diseases=rbind(later_diseases,c(NA,disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),1],as.character(disease_coding[which(disease_coding[,1]==disease_info[primary_entries[i],11+idx2[j]-1]),2]),disease_info[primary_entries[i],11+idx2[j]-1+58]))
					}
				}
			}
		}
}



#plot diseases before the main disease
MIN_COUNTS=200
sorted_previous_diseases=sort(table(previous_diseases[,3]),decreasing=T)
toplot_previous=which(sorted_previous_diseases>=MIN_COUNTS)
sorted_later_diseases=sort(table(later_diseases[,3]),decreasing=T)
toplot_later=which(sorted_later_diseases>=MIN_COUNTS)

#now to assess significance get prevalence of these diseases across all Uk Biobank individuals with the same mean age as the people who get Crohns.
percentage_early=matrix(ncol=5,nrow=length(toplot_previous),0);
colnames(percentage_early)=c("Prevalence","Pt_disease_before_maindisease_mean_age","Prevalence_before_meandisease_mean_age","people_with_disease_whogot_it_before_main_trait","people_with_disease_whogot_it_after_main_trait")
for(i in 1:length(toplot_previous)){
disease=disease_coding[which(disease_coding[,2]==names(sorted_previous_diseases)[i]),1]
disease_entries_aao=NULL
for(j in 11:39){
disease_entries_aao=c(disease_entries_aao,disease_info[which(disease_info[,j]==disease),j+58])
disease_entries_aao[disease_entries_aao<0]=NA
percentage_early[i,1]=length(disease_entries_aao)/dim(disease_info)[1]
percentage_early[i,2]=sum(disease_entries_aao<primean,na.rm=T)/length(disease_entries_aao)
percentage_early[i,3]=percentage_early[i,2]*percentage_early[i,1]
}
percentage_early[i,4]=sorted_previous_diseases[i]
percentage_early[i,5]=sorted_later_diseases[which(names(sorted_later_diseases)==names(sorted_previous_diseases[i]))]
}

percentage_late=matrix(ncol=5,nrow=length(toplot_later),0);
colnames(percentage_late)=c("Prevalence","Pt_disease_after_maindisease_mean_age","Prevalence_after_meandisease_mean_age","people_with_disease_whogot_it_after_main_trait","people_with_disease_whogot_it_before_main_trait")
for(i in 1:length(toplot_later)){
disease=disease_coding[which(disease_coding[,2]==names(sorted_later_diseases)[i]),1]
disease_entries_aao=NULL
for(j in 11:39){
disease_entries_aao=c(disease_entries_aao,disease_info[which(disease_info[,j]==disease),j+58])
disease_entries_aao[disease_entries_aao<0]=NA
percentage_late[i,1]=length(disease_entries_aao)/dim(disease_info)[1]
percentage_late[i,2]=sum(disease_entries_aao>primean,na.rm=T)/length(disease_entries_aao)
percentage_late[i,3]=percentage_late[i,2]*percentage_late[i,1]
}
percentage_late[i,4]=sorted_later_diseases[i]
percentage_late[i,5]=sorted_previous_diseases[which(names(sorted_previous_diseases)==names(sorted_later_diseases[i]))]
}

outfile_previous=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Co-occurent_diseases_with_",tmp,"_with_earlier_AAO_min",MIN_COUNTS,"_",yearwindow,".pdf",sep="")
if(printtopdf){pdf(outfile_previous,width=10,height=14)}
par(mfrow=c(2,1))
par(mar=c(15,5,5,2))
barplot(sorted_previous_diseases[toplot_previous],main=list(paste("Co-occurent diseases (excl. cancer) with\n ",Disease," that have earlier AAO",yearwindowtitle,sep=""),cex=1.1),col="#5B1A18",
ylab=list("Counts",cex=1.1),las=2)
mtext("Results from Initial Assessment",3)
TOPLOT=t(cbind(percentage_early[,4]/(percentage_early[,4]+percentage_early[,5]),percentage_early[,2]))
par(mar=c(2,5,2,2))
barplot(TOPLOT[,toplot_previous],beside=T,col=c("darkblue","red"),ylim=c(0,1),space=c(.25,1.75))
legend("topleft",pch=19,col=c("darkblue","red"),c(paste("% People who develop ",Disease," BEFORE the second trait",sep=""),paste("% People in UK Biobank who develop the second trait before the mean age for ",Disease,"(",primean,")",sep="")),cex=0.8)
if(printtopdf){dev.off()}

outfile_later=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Co-occurent_diseases_with_",tmp,"_with_later_AAO_min",MIN_COUNTS,"_",yearwindow,".pdf",sep="")
if(printtopdf){pdf(outfile_later,width=10,height=14)}
par(mfrow=c(2,1))
par(mar=c(15,5,5,2))
barplot(sorted_later_diseases[toplot_later],main=list(paste("Co-occurent diseases (excl. cancer) with\n ",Disease," that have later AAO",yearwindowtitle,sep=""),cex=1.1),col="#5B1A18",
ylab=list("Counts",cex=1.1),las=2)
mtext("Results from Initial Assessment",3)
TOPLOT=t(cbind(percentage_late[,4]/(percentage_late[,4]+percentage_late[,5]),percentage_late[,2]))
par(mar=c(2,5,2,2))
barplot(TOPLOT[,toplot_later],beside=T,col=c("darkblue","red"),ylim=c(0,1),space=c(.25,1.75))
legend("topleft",pch=19,col=c("darkblue","red"),c(paste("% People who develop ",Disease," AFTER the second trait",sep=""),paste("% People in UK Biobank who develop the second trait after the mean age for ",Disease,"(",primean,")",sep="")),cex=0.8)
if(printtopdf){dev.off()}


### A table having the before/after files for each disease is below ('FULL')
all_diseases=unique(sort(c(names(sorted_previous_diseases),names(sorted_later_diseases))))
FULL=matrix(ncol=3,nrow=length(all_diseases),0)
FULL[,1]=all_diseases
colnames(FULL)=c("Disease","CasesBefore","CasesAfter")
for(i in 1:dim(FULL)[1]){
num1=sorted_previous_diseases[which(names(sorted_previous_diseases)==FULL[i,1])]
if(length(num1)==1){FULL[i,2]=num1}
num2=sorted_later_diseases[which(names(sorted_later_diseases)==FULL[i,1])]
if(length(num2)==1){FULL[i,3]=num2}
}
 
which(abs(as.numeric(FULL[,2])-as.numeric(FULL[,3]))>=10)

#END