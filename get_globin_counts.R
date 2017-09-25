samples=read.table("/lustre/scratch115/projects/paxgene/gene_lists/mappings_ISR.txt",head=T,sep="\t") # this is used just for the samples list, so IU or ISR doesn't matter
counts=read.table("/lustre/scratch115/projects/paxgene/SalmonOutput_newruns/ISR.transcript_to_gene.counts",head=T) # this is used for the counts, so library name matters.
colnames(counts)=gsub('X2', '2', colnames(counts))
colnames(counts)=gsub('_IU.raw', '', colnames(counts))
oursamples=samples[!is.na(samples$Vol),] #removes spike-in samples (KC)
kapa_samples=as.character(oursamples[oursamples$kit=="kapa",1])
NEB_samples=as.character(oursamples[oursamples$kit=="NEB",1])
TruSeq_samples=as.character(oursamples[oursamples$kit=="TruSeq",1])

globin_genes=read.table("/lustre/scratch115/projects/paxgene/gene_lists/globin_genes_IDs.txt",skip=1)#read list of globin genes in
globin_idx=which(as.character(row.names(counts))%in%as.character(globin_genes[,1]))
globin_counts=counts[globin_idx,]

kapa_idx=which(colnames(globin_counts)%in%kapa_samples)
NEB_idx=which(colnames(globin_counts)%in%NEB_samples)
TruSeq_idx=which(colnames(globin_counts)%in%TruSeq_samples)
results=matrix(ncol=3,nrow=3,0)
rownames(results)=c("kapa","NEB","TruSeq")
colnames(results)=c("GlobinCounts","AllCounts","PercGlobin")
results[1,1:2]=c(sum(globin_counts[,kapa_idx]),sum(counts[,kapa_idx]))
results[2,1:2]=c(sum(globin_counts[,NEB_idx]),sum(counts[,NEB_idx]))
results[3,1:2]=c(sum(globin_counts[,TruSeq_idx]),sum(counts[,TruSeq_idx]))
results[,3]=results[,1]/results[,2]
results
#END