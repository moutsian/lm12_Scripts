
trait="UC"
utrait=toupper(trait)
#Vep_results_UC_suggestive_filtered_lowfreq.txt
data=read.table(paste("~lm12/IBD_conditional/Vep_results_",utrait,"_suggestive_filtered_lowfreq.txt",sep=""),head=T)
ourdata=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.txt",sep=""),head=T)
if(utrait=="IBD"){
colnames(ourdata)=c("chr","rsid","pos","allele_A","allele_B","P_value","coded_af","beta","se","Q","P_heterogeneity","I2","P_cohort_1","P_cohort_2","P_cohort_3","P_cohort_4")
}else{
colnames(ourdata)=c("chr","rsid","pos","allele_A","allele_B","P_value","coded_af","beta","se","Q","P_heterogeneity","I2","P_cohort_1","P_cohort_2","P_cohort_3")
}
tokeep=which(as.character(data[,4])=="splice_donor_variant" | as.character(data[,4])=="start_lost" | as.character(data[,4])=="protein_altering_variant" |  as.character(data[,4])=="transcript_amplification" | as.character(data[,4])=="transcript_ablation" | as.character(data[,4])=="splice_acceptor_variant" | as.character(data[,4])=="stop_gained" | as.character(data[,4])=="stop_lost" | as.character(data[,4])=="frameshift_variant" |as.character(data[,4])=="missense_variant" | as.character(data[,4])=="inframe_deletion" | as.character(data[,4])=="inframe_insertion"      )
tokeep_inourdata=matrix(ncol=1,nrow=length(tokeep),NA)
for(i in 1:length(tokeep)){
chr_pos=unlist(strsplit(unlist(strsplit(as.character(data[tokeep[i],2]),"-"))[1],":"))
 idx=which(as.character(ourdata[,1])==chr_pos[1] & as.character(ourdata[,3])==chr_pos[2])
	 if(length(idx)>0){
	 tokeep_inourdata[i]=idx
	 }else{
		chr=unlist(strsplit(unlist(strsplit(as.character(data[tokeep[i],2]),"-"))[1],":"))[1]
		pos=unlist(strsplit(as.character(data[tokeep[i],2]),"-"))[2]
		idx=which(as.character(ourdata[,1])==chr & as.character(ourdata[,3])==pos)
		if(length(idx)>0){
		 tokeep_inourdata[i]=idx	
		}else{
			#if none of the positions correspond to our variant, then look whether our variant is actually between the two positions
			chr_pos=unlist(strsplit(unlist(strsplit(as.character(data[tokeep[i],2]),"-"))[1],":"))
			pos=unlist(strsplit(as.character(data[tokeep[i],2]),"-"))[2]
			idx=which(as.character(ourdata[,1])==chr_pos[1] & as.numeric(as.character(ourdata[,3]))>as.numeric(chr_pos[2]) & as.numeric(as.character(ourdata[,3]))>as.numeric(pos)  )
			if(length(idx)>0){
			tokeep_inourdata[i]=idx[1]	
			}
		
		}
	}
 }
#OUT=cbind(ourdata[tokeep_inourdata,],data[tokeep,1:45],as.character(data[tokeep,46]),data[tokeep,47:50]) 
OUT=cbind(ourdata[tokeep_inourdata,],data[tokeep,]) 
write.table(OUT,paste("~lm12/IBD_conditional/Vep_results_",utrait,"_suggestive_filtered_lowfreq.high_and_moderate_impact_only.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

#END