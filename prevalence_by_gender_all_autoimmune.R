

library(gplots,lib.loc="~/Rpackages/")
printtopdf=T;
if(!exists("disease_info")){disease_info=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/ukb4940.tab",sep="\t",head=T)}
disease_coding=read.table("/lustre/scratch115/teams/anderson/ukbiobank/uncompressed_files/disease_coding.csv",sep=",",head=T)

Traits=c("multiple sclerosis","ulcerative colitis", "crohns disease","sclerosing cholangitis","systemic lupus erythematosis/sle","vitiligo","ankylosing spondylitis","sjogrens syndrome/sicca syndrome","psoriasis","alopecia / hair loss","malabsorption/coeliac disease","type 1 diabetes","sarcoidosis","myasthenia gravis","rheumatoid arthritis","acute infective polyneuritis/guillain-barre syndrome","primary biliary cirrhosis","graves disease","dermatopolymyositis","dermatomyositis")
"sjogren's syndrome/sicca syndrome"
#Disease="systemic lupus erythematosis/sle" #examples: "crohns disease" "irritable bowel syndrome" "sclerosing cholangitis" "malabsorption/coeliac disease" "infectious mononucleosis / glandular fever / epstein barr virus (ebv)"

DISEASES=matrix(ncol=7,nrow=length(Traits),0)
colnames(DISEASES)=c("primary_field","counts","counts_m","counts_f","mean_aao","counts_secondary","mean_aao_secondary")
rownames(DISEASES)=Traits
for(i in 1:dim(DISEASES)[1]){
primafield=disease_coding[which(disease_coding[,2]==Traits[i]),1]
if(length(primafield)==0){print(paste("Disease ",Traits[i],"could not be found"))}else{
if(length(primafield)==2){ #relevant for myasthenia gravis
DISEASES[i,1]=primafield[1]
DISEASES=rbind(DISEASES,c(primafield[2],rep(0,dim(DISEASES)[2]-1)))
}else{
DISEASES[i,1]=primafield
}
}}
rownames(DISEASES)[dim(DISEASES)[1]]="myasthenia gravis"


n_individuals=dim(disease_info)[1]
n_individuals_secondary=sum(!is.na(disease_info[,6]))


primary_entries=NULL
primary_entries_aao=NULL
secondary_entries=NULL
secondary_entries_aao=NULL

for( d in 1:dim(DISEASES)[1]){
for(i in 11:39){
DISEASES[d,2]=DISEASES[d,2]+length(which(disease_info[,i]==DISEASES[d,1]))
DISEASES[d,3]=DISEASES[d,3]+length(which(disease_info[,i]==DISEASES[d,1] & disease_info[,2]==1))
DISEASES[d,4]=DISEASES[d,4]+length(which(disease_info[,i]==DISEASES[d,1] & disease_info[,2]==0))
DISEASES[d,5]=DISEASES[d,5]+sum(disease_info[which(disease_info[,i]==DISEASES[d,1]),i+58])
}
for(i in 40:68){
DISEASES[d,6]=DISEASES[d,6]+length(which(disease_info[,i]==DISEASES[d,1]))
DISEASES[d,7]=DISEASES[d,7]+sum(disease_info[which(disease_info[,i]==DISEASES[d,1]),i+58])
}
}
DISEASES[,5]=DISEASES[,5]/DISEASES[,2]
DISEASES[,7]=DISEASES[,7]/DISEASES[,6]


#male to female ratio
colfunc=colorRampPalette(c("blue", "pink"))
COLS=colfunc(dim(DISEASES)[1])
outfile_ratio="/lustre/scratch115/teams/anderson/ukbiobank/plots/AID_gender_ratio.pdf"
if(printtopdf){pdf(outfile_ratio,width=10,height=10)}
RATIO=DISEASES[,3]/DISEASES[,4]
par(mar=c(15,5,5,2))
barplot(sort(RATIO,decreasing=T),las=2,main=list("Male to female ratio for various immune-mediated disorders",cex=1.1),ylab="male to female ratio",col=COLS)
abline(h=seq(0,1.5,by=0.25),col="gray",lty=2,lwd=2)
abline(h=1,col="red",lty=2,lwd=2)
if(printtopdf){dev.off()}


autoimmune_disease_counts=matrix(ncol=1,nrow=dim(DISEASES)[1],0)
#check total prevalence (Initial Assessment)
idx=NULL
for(i in 11:39){
idx=c(idx,which(disease_info[,i]%in%DISEASES[,1]))
}
autoimmune_disease_prevalence=length(unique(sort(idx)))/dim(disease_info)[1]
PREV=DISEASES[,2]/dim(disease_info)[1]
#now plot prevalence
outfile_prev="/lustre/scratch115/teams/anderson/ukbiobank/plots/AID_prevalence.pdf"
if(printtopdf){pdf(outfile_prev,width=10,height=10)}
par(mar=c(15,5,5,2))
barplot(sort(PREV,decreasing=T),col="darkblue",las=2,main=list("Prevalence by trait",cex=1.1))
mtext(paste("People in UKB with at least 1 trait in the list: ",autoimmune_disease_prevalence,sep=""),3,cex=0.8)
if(printtopdf){dev.off()}

cooccurence=table(table(idx))
outfile_cooc="/lustre/scratch115/teams/anderson/ukbiobank/plots/AID_cooccurence.pdf"
if(printtopdf){pdf(outfile_cooc,width=10,height=10)}
par(mar=c(5,5,5,2))
barplot(log(cooccurence),col="darkblue",las=2,main=list("Individuals with multiple AIDs",cex=1.1),ylab="log10(Counts)",xlab="# Traits per individual")
if(printtopdf){dev.off()}

#now get specific pairs
idx_multi=names(which(table(idx)>1))
COMORB=matrix(ncol=dim(DISEASES)[1],nrow=dim(DISEASES)[1],0)
colnames(COMORB)=DISEASES[,1]
rownames(COMORB)=DISEASES[,1]

for(i in 1:length(idx_multi)[1]){
	combination=which(colnames(COMORB)%in%disease_info[idx_multi[i],])
	if(length(combination)<2){
	print(paste("something went wrong, idx_multi: ",i,sep=""))
	}else{
		COMORB[combination[1],combination[2]]=COMORB[combination[1],combination[2]]+1
		COMORB[combination[2],combination[1]]=COMORB[combination[2],combination[1]]+1
		if(length(combination)>=3){ #all pairwise couples
		COMORB[combination[1],combination[3]]=COMORB[combination[1],combination[3]]+1
		COMORB[combination[3],combination[1]]=COMORB[combination[3],combination[1]]+1	
		COMORB[combination[2],combination[3]]=COMORB[combination[2],combination[3]]+1
		COMORB[combination[3],combination[2]]=COMORB[combination[3],combination[2]]+1		
		}
		if(length(combination>=4)){
		COMORB[combination[1],combination[4]]=COMORB[combination[1],combination[4]]+1
		COMORB[combination[4],combination[1]]=COMORB[combination[4],combination[1]]+1	
		COMORB[combination[2],combination[4]]=COMORB[combination[2],combination[4]]+1
		COMORB[combination[4],combination[2]]=COMORB[combination[4],combination[2]]+1			
		COMORB[combination[3],combination[4]]=COMORB[combination[3],combination[4]]+1
		COMORB[combination[4],combination[3]]=COMORB[combination[4],combination[3]]+1		
		}
		if(length(combination>=5)){
		COMORB[combination[1],combination[5]]=COMORB[combination[1],combination[5]]+1
		COMORB[combination[5],combination[1]]=COMORB[combination[5],combination[1]]+1	
		COMORB[combination[2],combination[5]]=COMORB[combination[2],combination[5]]+1
		COMORB[combination[5],combination[2]]=COMORB[combination[5],combination[2]]+1			
		COMORB[combination[3],combination[5]]=COMORB[combination[3],combination[5]]+1
		COMORB[combination[5],combination[3]]=COMORB[combination[5],combination[3]]+1		
		COMORB[combination[4],combination[5]]=COMORB[combination[4],combination[5]]+1
		COMORB[combination[5],combination[4]]=COMORB[combination[5],combination[4]]+1			
		}		
		if(length(combination>=6)){
		COMORB[combination[1],combination[6]]=COMORB[combination[1],combination[6]]+1
		COMORB[combination[6],combination[1]]=COMORB[combination[6],combination[1]]+1	
		COMORB[combination[2],combination[6]]=COMORB[combination[2],combination[6]]+1
		COMORB[combination[6],combination[2]]=COMORB[combination[6],combination[2]]+1			
		COMORB[combination[3],combination[6]]=COMORB[combination[3],combination[6]]+1
		COMORB[combination[6],combination[3]]=COMORB[combination[6],combination[3]]+1		
		COMORB[combination[4],combination[6]]=COMORB[combination[4],combination[6]]+1
		COMORB[combination[6],combination[4]]=COMORB[combination[6],combination[4]]+1			
		COMORB[combination[5],combination[6]]=COMORB[combination[5],combination[6]]+1
		COMORB[combination[6],combination[5]]=COMORB[combination[6],combination[5]]+1			
		}	
	}
}

## now plot COMORB

# i) Without accounting for prevalence

#ii) Accounting for prevalence
colnames(COMORB)=rownames(DISEASES)
rownames(COMORB)=rownames(DISEASES)
outfile_comorb="/lustre/scratch115/teams/anderson/ukbiobank/plots/AID_cooccurence_matrix.pdf"
if(printtopdf){pdf(outfile_comorb,width=15,height=15)}
par(oma=c(10,2,2,10))
MYTITLE="Co-morbidity between Traits"
heatmap.2(COMORB,Colv=NA,Rowv=NA,trace="none",dendrogram="none",
main=list(MYTITLE,cex=1.5),
density.info="none",colsep=seq(from=1,to=dim(COMORB)[2],by=1),
cellnote=COMORB,notecex=1.4,key=F,notecol="black",na.color=par("bg"),rowsep=seq(from=1,to=dim(COMORB)[1],by=1),
sepwidth=c(0.05,0.05),sepcolor="grey",col =colorRampPalette(c("white","red"))(100),margin=c(4,4))# cell labeling
if(printtopdf){dev.off()}


#ii) Accounting for prevalence
COMORBPC=COMORB/DISEASES[,2]
COMORBPC=round(COMORBPC,digits=2)
colnames(COMORBPC)=rownames(DISEASES)
rownames(COMORBPC)=rownames(DISEASES)
outfile_comorbpc="/lustre/scratch115/teams/anderson/ukbiobank/plots/AID_cooccurence_matrix_pt.pdf"
if(printtopdf){pdf(outfile_comorbpc,width=15,height=15)}
par(oma=c(10,2,2,10))
MYTITLE="Co-morbidity between Traits"
heatmap.2(COMORBPC,Colv=NA,Rowv=NA,trace="none",dendrogram="none",
main=list(MYTITLE,cex=1.5),
density.info="none",colsep=seq(from=1,to=dim(COMORBPC)[2],by=1),
cellnote=COMORBPC,notecex=1.4,key=F,notecol="black",na.color=par("bg"),rowsep=seq(from=1,to=dim(COMORBPC)[1],by=1),
sepwidth=c(0.05,0.05),sepcolor="grey",col =colorRampPalette(c("white","red"))(100),margin=c(4,4))# cell labeling
if(printtopdf){dev.off()}



#END