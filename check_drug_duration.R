ctrls_all=read.table("/lustre/scratch115/projects/crohns/exome/TIH/aux_and_intermediate_files/carl_controls_2017_07_19_duration_edited.txt",head=T,sep="\t")
cases_all=read.table("/lustre/scratch115/projects/crohns/exome/TIH/aux_and_intermediate_files/TILI_carl_cases_20170601_duration_edited.txt",head=T,sep="\t")
fam=read.table("/lustre/scratch115/projects/crohns/exome/TIH/TILI_merged.qc4_v2.pruned.fam")
#note that the above are before any filtering (due to PCA, high relatedness, missingness etc)
mean(cases_all$Duration.weeks.,na.rm=T)
mean(ctrls_all$azaDuration.weeks.,na.rm=T)
ctrls=ctrls_all[which(ctrls[,1]%in%fam[,1]),]
cases=cases_all[which(cases[,1]%in%fam[,1]),]
mean(cases$Duration.weeks.,na.rm=T)
mean(ctrls$azaDuration.weeks.,na.rm=T)
mean(ctrls$mpDuration.weeks.,na.rm=T)
summary(cases$Duration.weeks.,na.rm=T)
summary(ctrls$azaDuration.weeks.,na.rm=T)
summary(ctrls$mpDuration.weeks.)

ulc_col_cases=which(cases$crohns_diagnosis_date=="")
crohns_cases=which(cases$ulc_col_diagnosis_date=="")


#END