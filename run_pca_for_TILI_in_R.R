#
ped=read.table("TILI_merged.qc3.pruned.corrected.recoded.ped",stringsAsFactors=F)
popinfo=ped[,1:6]
data=ped[,-c(1:6)]

test=matrix(ncol=1,nrow=dim(data)[2],0)
for(i in 1:dim(data)[2]){
test[i]=length(table(data[,i]))
}
#now the remove the SNPs that are monomoprhic -- consider whether this has been automatically done for the EIGENSTRAT or not. Afto mporei na egine epeidi 
#evgala merika SNPs sto telos (pou ipirxan mono se cases)
datafilt=data[,-c(which(test[,1]==1))]
#now run pca function
tili.pca=prcomp(datafilt,center=TRUE,scale.=TRUE)
summary(tili.pca)

#END