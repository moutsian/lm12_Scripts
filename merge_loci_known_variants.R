##updated on Aprilt 4th to account for the bug with the last entries.

#ldvars_unord=as.matrix(read.table("/nfs/users/nfs_l/lm12/IBD_conditional/gw_vars_unique_reported.final.withLD.txt",head=T) )
ldvars_unord=as.matrix(read.table("/nfs/users/nfs_l/lm12/IBD_conditional/gw_vars_unique_reported.final.updated_with_Ellinghaus_sept16.withLD.txt",head=T) )
ldvars=ldvars_unord[order(as.numeric(ldvars_unord[,1]),as.numeric(ldvars_unord[,2])),]
noldleft=which(is.na(ldvars[,4]))
noldright=which(is.na(ldvars[,5]))
ldvars[noldleft,4]=ldvars[noldleft,2]
ldvars[noldright,5]=ldvars[noldright,2]

colnames(ldvars)=c("Chr","Pos_b37","SNP_ID","LD_left","LD_right")
previous=read.table("/nfs/users/nfs_l/lm12/IBD_conditional/IBD.CD.UC.GWAS.Hits.WithPositions.hg19.25Jan2016.LM.txt",head=T,sep="\t")
previousreported=previous[previous[,1]!="Novel",]
previousgw=previousreported[as.numeric(as.character(previousreported[,7]))<=5e-08,]

dataout=NULL
n_lead_var=1;
for(i in 2:dim(ldvars)[1]){

if( (as.numeric(ldvars[i,4])<=as.numeric(ldvars[i-1,5]) | abs(as.numeric(ldvars[i,2])-as.numeric(ldvars[i-1,2]))<500000 ) & (as.numeric(ldvars[i,1])==as.numeric(ldvars[i-1,1])) ){
ldvars[i,3]=as.matrix(unique(paste(as.character(ldvars[i,3]),as.character(ldvars[i-1,3]),sep=",")))
n_lead_var=n_lead_var+1;
ldvars[i,4]=min(as.numeric(ldvars[i,4]),as.numeric(ldvars[i-1,4]),na.rm=T)
ldvars[i,5]=max(as.numeric(ldvars[i,5]),as.numeric(ldvars[i-1,5]),na.rm=T)
}else{
dataout=rbind(dataout,c(ldvars[i-1,],n_lead_var))
n_lead_var=1; #reset
}
}
colnames(dataout)=c(colnames(ldvars),"num_of_previously_reported_gw_significant_variants_in_locus")
dataout=rbind(dataout,c(ldvars[i,],n_lead_var-1))

#now add additional info required
addedinfo=matrix(nrow=dim(dataout)[1],ncol=8,"NA")
colnames(addedinfo)=c("topSNP","topSNP_pos","topSNP_pval","topSNP_PUBMED","top_trait","all_pvals","all_PMIDs","all_traits")
for(i in 1:dim(dataout)[1]){
allvars=which(as.numeric(as.character(previousgw[,4]))>=as.numeric(as.character(dataout[i,4])) & as.numeric(as.character(previousgw[,4]))<=as.numeric(as.character(dataout[i,5])) & as.numeric(as.character(previousgw[,3]))==as.numeric(as.character(dataout[i,1])))
if(length(allvars)>0){
min_idx=which.min(previousgw[allvars,7])
addedinfo[i,1]=as.character(previousgw[allvars[min_idx],5])
addedinfo[i,2]=as.character(previousgw[allvars[min_idx],4])
addedinfo[i,3]=as.character(previousgw[allvars[min_idx],7])
addedinfo[i,4]=as.character(previousgw[allvars[min_idx],1])
addedinfo[i,5]=as.character(previousgw[allvars[min_idx],2])
addedinfo[i,6]=paste(as.character(previousgw[allvars,7]),sep=",",collapse=",")
addedinfo[i,7]=paste(as.character(previousgw[allvars,1]),sep=",",collapse=",")
addedinfo[i,8]=paste(as.character(previousgw[allvars,2]),sep=",",collapse=",")
}else{
print(paste("warning: No entry found for ",dataout[i,1],":",dataout[i,2])) #these are previously annotated "novel" loci which should be removed.
}
}
write.table(cbind(dataout,addedinfo),"/nfs/users/nfs_l/lm12/IBD_conditional/table_of_known_loci.updated_with_ellinghaus.sept16.txt",quote=F,col.names=T,row.names=F,sep="\t")

#END
