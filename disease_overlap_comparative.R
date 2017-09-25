
printtopdf=T;

MIN_COUNTS=15;
ASSESS="I" # "I" for initial, "S" for secondary
START=-1; # 11 for initial assessment, 40 for repeat assessment
END=-1; #39 for initial assesssment, 68 for repeat assessment

if(ASSESS=="I"){
ASSESSMENT="Initial Assessment (Field: 20002)"
START=11; # 11 for initial assessment, 40 for repeat assessment
END=39; #39 for initial assesssment, 68 for repeat assessment
ADD="";
}else{
ASSESSMENT="Secondary Assessment (Field:20002.1)"
START=40; # 11 for initial assessment, 40 for repeat assessment
END=68; #39 for initial assesssment, 68 for repeat assessment
ADD="_secondary_assessment"
}

if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/ukb4940.tab",sep="\t",head=T)}
disease_coding=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/disease_coding.csv",sep=",",head=T)
Disease="chronic fatigue syndrome" #examples: "crohns disease" "irritable bowel syndrome" "sclerosing cholangitis" "malabsorption/coeliac disease"
primafield=disease_coding[which(disease_coding[,2]==Disease),1]
if(length(primafield)==0){print(paste("Disease ",Disease,"could not be found"))}

ENTRIES=NULL
for(i in START:END){
ENTRIES=c(ENTRIES,which(disease_info[,i]==primafield))
}

n_individuals=dim(disease_info)[1]
if(ASSESS=="S"){
n_individuals=sum(!is.na(disease_info[,6]))
}

primary_entries=NULL
secondary_entries=NULL
for(i in 11:39){
primary_entries=c(primary_entries,which(disease_info[,i]==primafield))
}
for(i in 40:68){
secondary_entries=c(secondary_entries,which(disease_info[,i]==primafield))
}
newcases=which(!secondary_entries%in%primary_entries)

disease_filtered=disease_info[ENTRIES,START:END]
healthy_individuals=sum(is.na(disease_info[,START]))
n_diseases=table(rowSums(!is.na(disease_filtered))) #this vector contains the number of individuals with at least x diseases (e.g. n_diseases[4] is the number of ppl with 4 diseases)

disease_counts=table(unlist(disease_filtered))
sorted_counts=sort(disease_counts,decreasing=T)
print(paste("Total individuals with ",Disease,": ",sorted_counts[1],sep=""));

sorted_diseases=matrix(ncol=2,nrow=dim(sorted_counts)[1],NA)
for(i in 1:dim(sorted_diseases)[1]){
idx=which(disease_coding[,1]==as.numeric(names(sorted_counts[i])))
sorted_diseases[i,]=c(disease_coding[idx,1],as.character(disease_coding[idx,2]))
}


prevalence_in_uk_biobank=matrix(ncol=1,nrow=dim(sorted_diseases)[1],0)
for(i in START:END){
for(j in 1:dim(sorted_diseases)[1]){
prevalence_in_uk_biobank[j]=prevalence_in_uk_biobank[j]+length(which(disease_info[,i]==sorted_diseases[j,1]))
}}


prevalence_in_uk_biobank_pt=prevalence_in_uk_biobank/n_individuals #this should be divided with the number of people who had a secondary assessment
tmp=gsub(' ','_',Disease)
tmp=gsub('/','_',tmp)
tmp=gsub('\\(','_',tmp)
tmp=gsub('\\)','_',tmp)

#prepare variables to use downstream
PREV=cbind(prevalence_in_uk_biobank,sorted_counts)

prevtitle=paste("Prevalence_in_",tmp,sep="")
colnames(PREV)=c("Prevalence_UKBiobank",prevtitle) # to see what the diseases are, look at matrix sorted_diseases
PREVPT=PREV
PREVPT[,1]=PREVPT[,1]/n_individuals
PREVPT[,2]=PREVPT[,2]/PREV[1,2]
reldif=function(X){res=abs(X[1]-X[2])/max(X[1],X[2]);return(res)}
RELDIF=apply(PREVPT,1,reldif)
FULLMAT=cbind(sorted_diseases,PREV,PREVPT,RELDIF)
FULLMAT=FULLMAT[order(as.numeric(FULLMAT[,7]),decreasing=T),]
most_different=which(as.numeric(FULLMAT[,4])>MIN_COUNTS & as.numeric(FULLMAT[,6])<1)

rownames(FULLMAT)=FULLMAT[,2]

outfile_comp=paste("/lustre/scratch115/teams/anderson/ukbiobank/plots/","Diseases_min",MIN_COUNTS,"co_occuring_w_",tmp,"_excl_cancer",ADD,".pdf",sep="")
if(printtopdf){pdf(outfile_comp,width=12)}
par(mar=c(15,4,4,2))
TOPLOT=most_different[1:10]
barplot(matrix(as.numeric(t(FULLMAT[TOPLOT,6:5])),ncol=10),beside=T,space=c(.25,1), main=list(paste("Biggest difference in prevalence between individuals with ",Disease,"\n and all UK Biobank individuals (excl. cancer)",sep=""),cex=1),
col=c("darkblue","red"),ylab=list("Percentage (%) of Overlap",cex=1.1),ylim=c(0,max(as.numeric(FULLMAT[TOPLOT,6:5]))),las=2,cex.names=0.9,names.arg=rownames(FULLMAT[TOPLOT,]))
mtext(paste("Results from ",ASSESSMENT,", at least ",MIN_COUNTS," individuals with the secondary condition",sep=""),3,cex=0.7)
legend("topright",pch=19,col=c("darkblue","red"),c(paste("Overlap with ",Disease,sep=""),"Prevalence in UK Biobank"),cex=0.7)
if(printtopdf){dev.off()}




#END