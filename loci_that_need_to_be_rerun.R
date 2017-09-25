
trait="UC"
RERUNALL=NULL
dirin=paste("/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/",trait,"/",sep="")
for(chr in 1:22){
myfile=paste(dirin,"step1.filtered.",trait,".",chr,".txt",sep="")
if(file.exists(myfile)){
datain=read.table(myfile,sep="\t",head=T)
rerunleft=which(abs(as.numeric(datain[,7])-as.numeric(datain[,3]))>450000)
rerunright=which(abs(as.numeric(datain[,8])-as.numeric(datain[,3]))>450000)
RERUNL=NULL
RERUNR=NULL
if(!(identical(rerunright-500000,numeric(0)) & identical(rerunleft-500000,numeric(0))) ){
if(length(rerunright)>0){
RERUNR=datain[rerunright,1:8]
RERUNR[,8]=RERUNR[,8]+1000000
}
if(length(rerunleft)>0){
RERUNL=datain[rerunleft,1:8]
RERUNL[,7]=RERUNL[,7]-1000000
}
RERUN=rbind(RERUNL,RERUNR)
RERUN[duplicated(RERUN[,3])|duplicated(RERUN[,3],fromLast=TRUE),7]=min(RERUN[duplicated(RERUN[,3])|duplicated(RERUN[,3],fromLast=TRUE),7])
RERUN[duplicated(RERUN[,3])|duplicated(RERUN[,3],fromLast=TRUE),8]=max(RERUN[duplicated(RERUN[,3])|duplicated(RERUN[,3],fromLast=TRUE),8])
RERUN=unique(RERUN)
RERUNALL=rbind(RERUNALL,RERUN)
}}
}
write.table(RERUNALL,paste(dirin,"variants_to_rerun_ldcalc_for.",trait,".filtered.txt",sep=""),quote=F,row.names=F,col.names=T)

#END
