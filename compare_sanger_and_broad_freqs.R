broad=read.table("/lustre/scratch115/projects/crohns/exome/TIH/TILI_final_ctrls.qc1d.broad.frq",head=T)
sanger=read.table("/lustre/scratch115/projects/crohns/exome/TIH/TILI_final_ctrls.qc1d.sanger.frq",head=T)
#re-align
 disc=which(sanger[,3]!=broad[,3])
 sanger[disc,3] -> tmp
 sanger[disc,3]=sanger[disc,4]
 sanger[disc,4]=tmp
 sanger[disc,5]=1-sanger[disc,5]
 
  mrg=merge(sanger,broad,by.x="SNP",by.y="SNP")
  max_frq=apply(mrg[,c(5,10)],1,max)
  differ=abs(mrg[,5]-mrg[,10])
  diffratio=differ/max_frq
  threshold=0.4
  toremove=which(diffratio>threshold & max_frq>0.05)
  write.table(mrg[toremove,c(1:5,10:11)],"/lustre/scratch115/projects/crohns/exome/TIH/snps_to_remove_because_of_diff_frq_in_sanger_vs_broad.qc1d.txt",quote=F,col.names=T,row.names=F)
 #END