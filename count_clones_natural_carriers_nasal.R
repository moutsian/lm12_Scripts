#
directory="/lustre/scratch115/projects/paxgene/tcrseq/pneumo/"
files_minus=c("19828_2_1.clones.txt","19828_2_5.clones.txt","19828_3_5.clones.txt","19828_5_1.clones.txt","19828_2_4.clones.txt")
files_plus=c("19828_2_2.clones.txt","19828_3_1.clones.txt","19828_4_1.clones.txt","19828_4_2.clones.txt","19828_5_2.clones.txt")


clones_minus=read.table(paste(directory,files_minus[1],sep=""),head=T,sep="\t")
added=sum(clones_minus[,1])
for(i in 2:length(files_minus)){
clones_minus=rbind(clones_minus,read.table(paste(directory,files_minus[i],sep=""),head=T,sep="\t"))
added=c(added,sum(clones_minus[,1]))
}
sum(clones_minus[,1])

clones_plus=read.table(paste(directory,files_plus[1],sep=""),head=T,sep="\t")
for(i in 2:length(files_minus)){
clones_plus=rbind(clones_plus,read.table(paste(directory,files_plus[i],sep=""),head=T,sep="\t"))
}
sum(clones_plus[,1])

all_unique_clones=unique(c(as.character(clones_plus[,2]),as.character(clones_minus[,2])))
#now count occurences
clone_counter=matrix(ncol=6,nrow=length(all_unique_clones),0)
row.names(clone_counter)=all_unique_clones
colnames(clone_counter)=c("minus_5","plus_2","minus5_perc_of_total","plus2_perc_of_total","sum_of_clones","relative_difference")
for(i in 1:length(all_unique_clones)){

idx_min=which(as.character(clones_minus[,2])==all_unique_clones[i])
if(length(idx_min)>0){
clone_counter[i,1]=sum(clones_minus[idx_min,1])
}

idx_plus=which(as.character(clones_plus[,2])==all_unique_clones[i])
if(length(idx_plus)>0){
clone_counter[i,2]=sum(clones_plus[idx_plus,1])
}
}
clone_counter[,3]=clone_counter[,1]/sum(clone_counter[,1])
clone_counter[,4]=clone_counter[,2]/sum(clone_counter[,2])
clone_counter[,5]=clone_counter[,2]+clone_counter[,1]
clone_counter[,6]=abs(clone_counter[,2]-clone_counter[,1])/(clone_counter[,2]+clone_counter[,1])
clone_counter_sorted=clone_counter[order(clone_counter[,6],decreasing=T),]
ccs=clone_counter_sorted[which(clone_counter_sorted[,5]>40),]
write.table(clone_counter_sorted,paste(directory,"v_genes_counts_in_natural_carriers_nasal.txt",sep=""),row.names=T,col.names=T,quote=F,sep="\t")
#END