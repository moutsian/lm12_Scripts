# add_freq_and_I2 to the list, as well as an indication of how many datasets the result is based on
# -- note that for dataset C the update of the function check_if_previously_reported.R does it automatically - but I haven't yet configured that for dataset B.
#data=read.table( "~lm12/IBD_conditional/IIBDGC-meta-analysis-data-loci-20160211.txt",sep="\t",head=T)
trait="IBD"
n_alldat=3
if(toupper(trait)=="IBD"){
n_alldat=4
}

data=read.table(paste("~lm12/IBD_conditional/step2.",toupper(trait),".all.with_info_on_previous_hits_from_all_PAPERS.datasetB.txt",sep=""),sep="\t",head=T)
data_plus=matrix(nrow=dim(data)[1],ncol=dim(data)[2]+3,"")
data_plus[,1:dim(data)[2]]=as.matrix(data)
colnames(data_plus)=c(colnames(data),"FREQ_in_1KG_CEU_GBR","I2","datasets_used")
for(i in 1:dim(data)[1]){
shellcommand=paste(" sh /lustre/scratch115/teams/anderson/ibd_conditional/get_freq_and_I2_datasetB.sh ",toupper(trait),as.character(data[i,1]),as.character(data[i,3]),sep=" ")
shellout=unlist(strsplit(system(shellcommand,intern=T)," ")) # the last entry is how many datasets we do NOT have data for, for this variant. 
#For IBD this can be a max of 4 and for UC,CD a max of 3, although of course we have a GW p-value already so in practice that will be 3 and respectively.
if(length(shellout==9)){
data_plus[i,dim(data)[2]+1]=shellout[7]
data_plus[i,dim(data)[2]+2]=shellout[8]
data_plus[i,dim(data)[2]+3]=n_alldat - as.numeric(shellout[9])
}else{
print("Something went wrong!")
}
}
write.table(data_plus,paste("~lm12/IBD_conditional/Step2.",toupper(trait),".all.with_info_on_previous_hits_from_all_PAPERS.datasetB.txt.plus",sep=""),quote=F,sep="\t",row.names=F,col.names=T)

#END