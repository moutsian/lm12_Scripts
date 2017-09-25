
#this R file merges the step2 files across traits, and checks whether the hits we find in our study are previously reported or not. If they are not, it adds I2 and frequency info and outputs them separately

myfilename="~lm12/IBD_conditional/Vep_results_datasetB.ALL_suggestive_filtered_lowfreq.high_and_moderate_impact_only.with_geneinfo.txt"
vepres=as.matrix(read.table(myfilename,sep="\t",head=T))


vep_plus=matrix(nrow=dim(vepres)[1],ncol=dim(vepres)[2]+9,"")
vep_plus[,1:dim(vepres)[2]]=as.matrix(vepres)
colnames(vep_plus)=c(colnames(vepres),"MAF_GWAS1","MAF_GWAS2","MAF_GWAS3","MAF_IBDseq","pval_ichip","pval_IIBDGC","pval_DatC","FREQ_in_1KG_CEU_GBR","I2")
for(i in 1:dim(vepres)[1]){
print(paste("*****chr:",vepres[i,1],"pos:",vepres[i,3]))
shellcommand=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_freq_and_I2_B.sh ",toupper(trait),as.character(vepres[i,1]),as.character(vepres[i,3]),sep=" ")
shellout=unlist(strsplit(system(shellcommand,intern=T)," "))
shellfreqin=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_our_freq_and_ichip_pval_B.sh ",toupper(trait2),as.character(vepres[i,1]),as.character(vepres[i,3]),sep=" ")
shellfreqout=unlist(strsplit(system(shellfreqin,intern=T),","))
print(paste("shellfreqout:",shellfreqout))
if(length(shellfreqout)==8){
vep_plus[i,dim(vepres)[2]+1]=unlist(strsplit(shellfreqout[2],":"))[2]
vep_plus[i,dim(vepres)[2]+2]=unlist(strsplit(shellfreqout[3],":"))[2]
vep_plus[i,dim(vepres)[2]+3]=unlist(strsplit(shellfreqout[4],":"))[2]
vep_plus[i,dim(vepres)[2]+4]=unlist(strsplit(shellfreqout[5],":"))[2]
vep_plus[i,dim(vepres)[2]+5]=unlist(strsplit(shellfreqout[6],":"))[2]
vep_plus[i,dim(vepres)[2]+6]=unlist(strsplit(shellfreqout[7],":"))[2]
vep_plus[i,dim(vepres)[2]+7]=unlist(strsplit(shellfreqout[8],":"))[2]
}else{
print ("Something went wrong when getting frequencies from our data")
}

if(length(shellout==8)){
vep_plus[i,dim(vepres)[2]+8]=shellout[7]
vep_plus[i,dim(vepres)[2]+9]=shellout[8]
}else{
print("Something went wrong!")
}
}

Trait=rep(trait,dim(vep_plus)[1])
vep_plus=cbind(Trait,vep_plus)
write.table(vep_plus,"~lm12/IBD_conditional/Vep_results_datasetB.round6.ALL.high_and_moderate_impact_only.txt",sep="\t",col.names=T,quote=F,row.names=F)
vepresplus_all[10,c(1,3:9,11,12:16)]




#END