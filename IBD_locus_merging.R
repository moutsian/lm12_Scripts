mydir="/lustre/scratch115/teams/anderson/ibd_conditional/"
datain= as.matrix(read.table(paste(mydir,"loci_list_before_any_merging_r2_0.60.ALL.txt",sep=""),head=T,sep="\t"))

#before starting the merging process group together the "Trait" categories into 3: IBD, UC and CD
supergroup=matrix(ncol=1,nrow=dim(datain)[1],NA)
supergroup[datain[,2]=="Chrons" | datain[,2]=="Chrons_related" | datain[,2]=="Chrons_and_celiac" | datain[,2]=="Chrons_time_to_surgery" | datain[,2]=="Chrons_need_for_surgery" | datain[,2]=="Chrons_and_psoriasis" | datain[,2]=="Chrons_combined" | datain[,2]=="Chrons_and_sarcoidosis_combined"]="CD"
supergroup[datain[,2]=="IBD" | datain[,2]=="IBD_unsaturated" | datain[,2]=="IBD_saturated" | datain[,2]=="IBD_SAT" | datain[,2]=="IBD (early onset)"| datain[,2]=="UC_or_Crohns"] ="IBD"
supergroup[datain[,2]=="UC"]="UC"

trait="IBD" # can be "UC" or "IBD" or "CD"

data_trait=cbind(datain[which(supergroup==trait),],supergroup[which(supergroup==trait)]) #do the same for the other two once this works.
idx_cnt=100;

for(i in 1:(dim(data_trait)[1])){
idx=which((as.numeric(data_trait[,11])>as.numeric(data_trait[i,10]))& (as.numeric(data_trait[,11])<as.numeric(data_trait[i,11])) & (as.numeric(data_trait[,3])==as.numeric(data_trait[i,3])))
if(length(idx)>0){
idx_cnt=idx_cnt+1
print(paste("i:",i,sep=""))
endpos=max(as.numeric(data_trait[idx,11]),as.numeric(data_trait[i,11]))
data_trait[idx,11]=endpos
if(endpos>as.numeric(data_trait[i,10])+2000000){
print(paste("Something is wrong for i: ",i,sep=""))
break}
}
}


print(paste("idx_cnt: ",idx_cnt,sep=""));

#now change the left side too.
for(i in 1:(dim(data_trait)[1])){
idx=which((as.numeric(data_trait[,11])==as.numeric(data_trait[i,11]))& (as.numeric(data_trait[,3])==as.numeric(data_trait[i,3])))
if(length(idx)>0){
print(paste("i:",i,sep=""))
startpos=min(as.numeric(data_trait[idx,10]))
data_trait[i,10]=startpos
}
}


#sort loci:
ordIBD=data_trait[order(data_trait[,3],data_trait[,10]),]
ordIBD[1:50,c(1,3,5,10,11)] #quick view of the data
novel_idx=which(ordIBD[,1]=="Novel")
is_novel=matrix(ncol=1,nrow=length(novel_idx),NA)
for(i in 1:length(novel_idx)){
dups=which(as.numeric(ordIBD[,10])==as.numeric(ordIBD[novel_idx[i],10]))
if(length(dups)>1){
print("This is not novel:")
is_novel[i]=-1 
ordIBD[dups,]
}else{
is_novel[i]=1}
}

#now merge
uloci=rep(NA,14)
uloci=rbind(uloci,ordIBD[1,])
for(i in 2:dim(ordIBD)[1]){
if(is.na(ordIBD[i,10])){
uloci=rbind(uloci,ordIBD[i,])
}else if( as.numeric(ordIBD[i,10])==as.numeric(ordIBD[i-1,10]) & as.numeric(ordIBD[i,11])==as.numeric(ordIBD[i-1,11])){
uloci[dim(uloci)[1],1]=paste(uloci[dim(uloci)[1],1],",",ordIBD[i,1],sep="")
uloci[dim(uloci)[1],2]=paste(uloci[dim(uloci)[1],2],",",ordIBD[i,2],sep="")
uloci[dim(uloci)[1],4]=paste(uloci[dim(uloci)[1],4],",",ordIBD[i,4],sep="")
uloci[dim(uloci)[1],5]=paste(uloci[dim(uloci)[1],5],",",ordIBD[i,5],sep="")
uloci[dim(uloci)[1],6]=paste(uloci[dim(uloci)[1],6],",",ordIBD[i,6],sep="")
uloci[dim(uloci)[1],7]=paste(uloci[dim(uloci)[1],7],",",ordIBD[i,7],sep="")
uloci[dim(uloci)[1],8]=paste(uloci[dim(uloci)[1],8],",",ordIBD[i,8],sep="")
uloci[dim(uloci)[1],9]=paste(uloci[dim(uloci)[1],9],",",ordIBD[i,9],sep="")
uloci[dim(uloci)[1],12]=paste(uloci[dim(uloci)[1],12],",",ordIBD[i,12],sep="")
uloci[dim(uloci)[1],13]=paste(uloci[dim(uloci)[1],13],",",ordIBD[i,13],sep="")
}else{
uloci=rbind(uloci,ordIBD[i,])
}
}

colnames(uloci)[14]="Trait_Simplified"
write.table(uloci,paste(mydir,"loci_list_unique.",trait,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)



#END