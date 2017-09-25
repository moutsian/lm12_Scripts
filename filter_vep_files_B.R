
trait="IBD"
utrait=toupper(trait)
#Vep_results_UC_suggestive_filtered_lowfreq.txt
data=read.table(paste("~lm12/IBD_conditional/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.txt.annotated",sep=""),head=T,sep="\t")
datarefseq=read.table(paste("~lm12/IBD_conditional/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.txt.annotated.refseq",sep=""),head=T,sep="\t")

ourdata=read.table(paste("/lustre/scratch115/teams/anderson/ibd_conditional/UK_only_analysis_B_results/UK_only_analysis_B.",utrait,".1e_05.filtered.lowfreq.txt",sep=""),head=T)
if(utrait=="IBD"){
colnames(ourdata)=c("chr","rsid","pos","allele_A","allele_B","P_value","coded_af","beta","se","Q","P_heterogeneity","I2","P_cohort_1","P_cohort_2","P_cohort_3","P_cohort_4")
}else{
colnames(ourdata)=c("chr","rsid","pos","allele_A","allele_B","P_value","coded_af","beta","se","Q","P_heterogeneity","I2","P_cohort_1","P_cohort_2","P_cohort_3")
}

#tokeep=which(as.character(data[,4])=="splice_donor_variant" | as.character(data[,4])=="start_lost" | as.character(data[,4])=="protein_altering_variant" |  as.character(data[,4])=="transcript_amplification" | as.character(data[,4])=="transcript_ablation" | as.character(data[,4])=="splice_acceptor_variant" | as.character(data[,4])=="stop_gained" | as.character(data[,4])=="stop_lost" | as.character(data[,4])=="frameshift_variant" |as.character(data[,4])=="missense_variant" | as.character(data[,4])=="inframe_deletion" | as.character(data[,4])=="inframe_insertion"      )

tokeep=matrix(ncol=1,nrow=dim(data)[1],0)
tmp=unlist(regexpr('splice_donor_variant', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('start_lost', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('protein_altering_variant', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('transcript_amplification', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('transcript_ablation', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('splice_acceptor_variant', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('stop_gained', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('stop_lost', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('frameshift_variant', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('missense_variant', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('inframe_deletion', datarefseq[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('inframe_insertion', datarefseq[,7]))
tokeep[tmp!=-1]=1
#same for ensembl data
tmp=unlist(regexpr('splice_donor_variant', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('start_lost', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('protein_altering_variant', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('transcript_amplification', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('transcript_ablation', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('splice_acceptor_variant', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('stop_gained', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('stop_lost', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('frameshift_variant', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('missense_variant', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('inframe_deletion', data[,7]))
tokeep[tmp!=-1]=1
tmp=unlist(regexpr('inframe_insertion', data[,7]))
tokeep[tmp!=-1]=1
idxkeep=which(tokeep==1)
tokeep_inourdata=matrix(ncol=1,nrow=length(idxkeep),NA)
for(i in 1:length(idxkeep)){
chr_pos=unlist(strsplit(as.character(data[idxkeep[i],1]),"_"))
 idx=which(as.character(ourdata[,1])==chr_pos[1] & as.character(ourdata[,3])==chr_pos[2])
	 if(length(idx)>0){
	 tokeep_inourdata[i]=idx
	 }
 }
OUT=cbind(ourdata[tokeep_inourdata,],data[idxkeep,],datarefseq[idxkeep,]) 
write.table(OUT,paste("~lm12/IBD_conditional/Vep_results_",utrait,"_suggestive_filtered_lowfreq.high_and_moderate_impact_only.with_geneinfo.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

#END