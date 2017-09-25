inputpathname = "/lustre/scratch115/projects/crohns/exome/TIH/TILI_final_ctrls.qc1c"
inputpathnamecases = "/lustre/scratch115/projects/crohns/exome/TIH/TILI_final_cases.qc1c"
hetlinesd = 3 #sd from mean heterozygosity rate
missingnessline = 0.05
missingnessline2 = 0.1

print(paste("plotting with sd =", hetlinesd, "and miss =s", missingnessline))

imiss = read.table(paste(inputpathname, ".imiss", sep=""), h=T)
imiss$logF_MISS = log10(imiss[,6])

imiss_cases = read.table(paste(inputpathnamecases, ".imiss", sep=""), h=T)
imiss_cases$logF_MISS = log10(imiss_cases[,6])


het = read.table(paste(inputpathname, ".het", sep=""), h=T)
het$meanHet = (het$N.NM. - het$O.HOM.) / het$N.NM.
het_upper=mean(het$meanHet) + (hetlinesd*sd(het$meanHet))
het_lower=mean(het$meanHet) - (hetlinesd*sd(het$meanHet))

het_cases = read.table(paste(inputpathnamecases, ".het", sep=""), h=T)
het_cases$meanHet = (het_cases$N.NM. - het_cases$O.HOM.) / het_cases$N.NM.
het_cases_upper=mean(het_cases$meanHet) + (hetlinesd*sd(het_cases$meanHet))
het_cases_lower=mean(het_cases$meanHet) - (hetlinesd*sd(het_cases$meanHet))

het_all=rbind(het,het_cases)
het_all$meanHet = (het_all$N.NM. - het_all$O.HOM.) / het_all$N.NM.
het_all_upper=mean(het_all$meanHet) + (hetlinesd*sd(het_all$meanHet))
het_all_lower=mean(het_all$meanHet) - (hetlinesd*sd(het_all$meanHet))

cases_to_remove=unique(c(which(het_cases$meanHet>het_all_upper),which(het_cases$meanHet<het_all_lower),which(imiss_cases$F_MISS>missingnessline)))
ctrls_to_remove=unique(c(which(het$meanHet>het_all_upper),which(het$meanHet<het_all_lower),which(imiss$F_MISS>missingnessline)))
write.table(het_cases[cases_to_remove,c(1:2)],"/lustre/scratch115/projects/crohns/exome/TIH/cases_to_remove.qc1c.txt",quote=F,col.names=T,row.names=F)
write.table(het[ctrls_to_remove,c(1:2)],"/lustre/scratch115/projects/crohns/exome/TIH/ctrls_to_remove.qc1c.txt",quote=F,col.names=T,row.names=F)
#END