## June 2012 - Loukas Moutsianas

setwd("/nfs/team152/loukas/from_oxford/MS2012/")
#Example input.
#file="PP_DPB1.txt"
#allele_res="4D"

prep_alleleprobs2012 = function(file,allele_res){

INPUTFILE=read.table(file,head=T,sep="\t")
dim(INPUTFILE)
colnames(INPUTFILE)
ALLELES= c(sapply(strsplit(colnames(INPUTFILE)[-1], "_"), "[[",2 ), sapply(strsplit(colnames(INPUTFILE)[-1], "_"), "[[",1 ))
ALLELES=gsub("X","",ALLELES)
if(allele_res=="2D"){
ALLELES=apply(as.matrix(ALLELES),1,substr,1,2)}

ALLELES=ALLELES[!duplicated(ALLELES)] #now, ALLELES containg a list of all unique alleles in the original file.
NEWFILE=array(0,dim=c(dim(INPUTFILE)[1],3,length(ALLELES)),
dimnames=list(INPUTFILE[,1],c("homs(0copies)","hets(1copy)","homs'(2copies)"),ALLELES)) #three dimensions: samples, probs, alleles.

for (i in 1:length(ALLELES)){
	
	#HETS=grep(ALLELES[i],colnames(INPUTFILE))
	HETS1=grep(paste("X",ALLELES[i],sep=""),colnames(INPUTFILE))
	HETS2=grep(paste("_",ALLELES[i],sep=""),colnames(INPUTFILE))
	HETS=unique(c(HETS1,HETS2))
	#HOM= grep(paste("X",ALLELES[i],"_",ALLELES[i],sep=""),colnames(INPUTFILE))
	HOM=HETS1[which(HETS1%in%HETS2)]
	if(length(HOM)==0){ #if there is no entry for homozygotes for this allele (cld happen for rare ones), set it to 0 and ignore it.
		#HOM=0; #do nothing.
	}else{
		# correct HETS indices by removing the HOM index.
		HETS=HETS[-which(HETS==HOM)]
	}
	NEWFILE[,2,i]=apply(INPUTFILE[,HETS],1,sum)
	if(length(HOM)>0){
		if(length(HOM)>1){ #this will be needed in 2D analysis
			NEWFILE[,1,i]=apply(INPUTFILE[,HOM],1,sum)
		}else{
			NEWFILE[,1,i]=INPUTFILE[,HOM]	
		}
	}
	NEWFILE[,3,i]=1-NEWFILE[,2,i]-NEWFILE[,1,i]
}

NEWFILE[NEWFILE<0]=0; # to avoid any issues arising due to rounding.
NEWFILE[NEWFILE>1]=1;

#now prepare output file: rows equal to number of samples for which we have imputations, columns are three: P_hom, P_het, P_hom' .

#colnames(NEWFILE)=c("homs(0copies)","hets(1copy)","homs'(2copies)") #not sure how to name cols of 3d matrix....

return(NEWFILE);
}

#END



