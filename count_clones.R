#
directory="/lustre/scratch115/projects/paxgene/tcrseq/pneumo/"
files_minus=c("19743_3_3.clones.txt","19743_4_2.clones.txt","19743_5_5.clones.txt","19743_6_2.clones.txt","19743_7_3.clones.txt","19743_4_4.clones.txt")
files_plus=c("19743_3_4.clones.txt","19743_4_3.clones.txt","19743_4_5.clones.txt","19743_6_1.clones.txt","19743_6_3.clones.txt","19743_7_4.clones.txt")

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
clone_counter=matrix(ncol=2,nrow=length(all_unique_clones),0)
row.names(clone_counter)=all_unique_clones
col.names(clone_counter)=c("minus_5","plus_2")
for(i in 1:length(all_unique_clones)){
}
#END