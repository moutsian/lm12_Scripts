
printtopdf=T;
if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/ukb4940.tab",sep="\t",head=T)}
disease_coding=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/disease_coding.csv",sep=",",head=T)
Disease="type 1 diabetes" #examples: "crohns disease" "irritable bowel syndrome" "sclerosing cholangitis" "malabsorption/coeliac disease" "infectious mononucleosis / glandular fever / epstein barr virus (ebv)"
primafield=disease_coding[which(disease_coding[,2]==Disease),1]
if(length(primafield)==0){print(paste("Disease ",Disease,"could not be found"))}

ADD=""
if(ASSESS=="I"){
ASSESSMENT="Initial Assessment (Field: 20002)"
START=11; # 11 for initial assessment, 40 for repeat assessment
END=39; #39 for initial assesssment, 68 for repeat assessment
}else{
ASSESSMENT="Secondary Assessment (Field:20002.1)"
START=40; # 11 for initial assessment, 40 for repeat assessment
END=68; #39 for initial assesssment, 68 for repeat assessment
ADD="_secondary_assessment"
}

n_individuals=dim(disease_info)[1]
if(ASSESS=="S"){
n_individuals=sum(!is.na(disease_info[,6]))
}

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

#gender is field 31. 0 is female, 1 male.
preva=dim(disease_info[primary_entries,])[1]/n_individuals
preva_fems=sum(disease_info[primary_entries,2]==0)/sum(disease_info[,2]==0)
preva_males=sum(disease_info[primary_entries,2]==1)/sum(disease_info[,2]==1)

outfile_later=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Prevalence_by_gender_",tmp,ADD,".pdf",sep="")
if(printtopdf){pdf(outfile_later,width=10,height=14)}
par(mfrow=c(1,1))
par(mar=c(2,5,5,2))
PREVA=c(preva,preva_fems,preva_males)
barplot(PREVA,main=list(paste("Prevalence of ",Disease," by gender",sep=""),cex=1.1),col=c("#5B1A18","pink","blue"),ylab=list("Prevalence",cex=1.1),names.arg=c("All","F","M"), ylim=c(0,max(PREVA)+0.05*max(PREVA)))
mtext("Results from Initial Assessment",3)
if(printtopdf){dev.off()}

#now by age
hist(primary_entries_aao)


#END