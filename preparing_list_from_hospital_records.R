
if(!exists("data")){
data=read.table("ukb5922.tab",head=T,sep="\t") #note this takes forever to load
}
diagnosis=matrix(ncol=3,nrow=dim(data)[1],0)
diagnosis[,1]=data[,1]
 colnames(diagnosis)=c("ID","Diagnosis","Field")
 for(i in 349:724){
idxA=which(as.character(data[,i])=="K900")
idxB=which(as.character(data[,i])=="K901")
idxC=which(as.character(data[,i])=="K908")
idxD=which(as.character(data[,i])=="K909")
#idxE=which(as.character(data[,i])=="K514")
#idxF=which(as.character(data[,i])=="K515")
#idxG=which(as.character(data[,i])=="K518")
#idxH=which(as.character(data[,i])=="K519")	
diag_field="41202"
if(length(idxA)>0){
	diagnosis[idxA,2]="K900"
	diagnosis[idxA,3]=diag_field
	}
	if(length(idxB)>0){
	diagnosis[idxB,2]="K901"
	diagnosis[idxB,3]=diag_field
	}
	if(length(idxC)>0){
	diagnosis[idxC,2]="K908"
	diagnosis[idxC,3]=diag_field
	}	
	if(length(idxD)>0){
	diagnosis[idxD,2]="K909"
	diagnosis[idxD,3]=diag_field
	}	
	# if(length(idxE)>0){
	# diagnosis[idxE,2]="K514"
	# diagnosis[idxE,3]=diag_field
	# }
	# if(length(idxF)>0){
	# diagnosis[idxF,2]="K515"
	# diagnosis[idxF,3]=diag_field
	# }	
	# if(length(idxG)>0){
	# diagnosis[idxG,2]="K518"
	# diagnosis[idxG,3]=diag_field
	# }	
	# if(length(idxH)>0){
	# diagnosis[idxH,2]="K519"
	# diagnosis[idxH,3]=diag_field
	# }	
}
for(i in 726:1069){
idxA=which(as.character(data[,i])=="K900")
idxB=which(as.character(data[,i])=="K904")
idxC=which(as.character(data[,i])=="K908")
idxD=which(as.character(data[,i])=="K909")
#idxE=which(as.character(data[,i])=="K514")
#idxF=which(as.character(data[,i])=="K515")
#idxG=which(as.character(data[,i])=="K518")
#idxH=which(as.character(data[,i])=="K519")	
diag_field="41204"
if(length(idxA)>0){
	diagnosis[idxA,2]="K900"
	diagnosis[idxA,3]=diag_field
	}
	if(length(idxB)>0){
	diagnosis[idxB,2]="K904"
	diagnosis[idxB,3]=diag_field
	}
	if(length(idxC)>0){
	diagnosis[idxC,2]="K908"
	diagnosis[idxC,3]=diag_field
	}	
	if(length(idxD)>0){
	diagnosis[idxD,2]="K909"
	diagnosis[idxD,3]=diag_field
	}	
	# if(length(idxE)>0){
	# diagnosis[idxE,2]="K514"
	# diagnosis[idxE,3]=diag_field
	# }
	# if(length(idxF)>0){
	# diagnosis[idxF,2]="K515"
	# diagnosis[idxF,3]=diag_field
	# }	
	# if(length(idxG)>0){
	# diagnosis[idxG,2]="K518"
	# diagnosis[idxG,3]=diag_field
	# }	
	# if(length(idxH)>0){
	# diagnosis[idxH,2]="K519"
	# diagnosis[idxH,3]=diag_field
	# }	
}
head(diagnosis[which(diagnosis[,2]!="0"),])




write.table(diagnosis,"Coeliac_malabsorption.status_from_hospital_admission.csv",col.names=T,row.names=F,quote=F,sep=",")

