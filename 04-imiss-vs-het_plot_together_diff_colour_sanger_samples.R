inputpathnamectrls = "/lustre/scratch115/projects/crohns/exome/TIH/TILI_final_ctrls.qc1"
inputpathnamecases = "/lustre/scratch115/projects/crohns/exome/TIH/TILI_final_cases.qc1"
sangersamples      = "/lustre/scratch115/projects/crohns/exome/TIH/TILI_ctrls_from_sanger"
sangersamples_preqc      = "/lustre/scratch115/projects/crohns/exome/TIH/TILI_ctrls_from_sanger.preqc"
outputpathname = "/lustre/scratch115/projects/crohns/exome/TIH/TILI_final.qc1.multicolour.pdf"

hetlinesd = 3
hetlinesd2= 4

missingnessline = 0.05
missingnessline2 = 0.1

print(paste("plotting with sd =", hetlinesd, "and miss =s", missingnessline))

imiss = read.table(paste(inputpathnamectrls, ".imiss", sep=""), h=T)
imiss$logF_MISS = log10(imiss[,6])


imiss_sanger = read.table(paste(sangersamples, ".imiss", sep=""), h=T)
imiss_sanger$logF_MISS = log10(imiss_sanger[,6])

imiss_sanger_preqc = read.table(paste(sangersamples_preqc, ".imiss", sep=""), h=T)
imiss_sanger_preqc$logF_MISS = log10(imiss_sanger_preqc[,6])


imiss_cases = read.table(paste(inputpathnamecases, ".imiss", sep=""), h=T)
imiss_cases$logF_MISS = log10(imiss_cases[,6])


het = read.table(paste(inputpathnamectrls, ".het", sep=""), h=T)
het_cases = read.table(paste(inputpathnamecases, ".het", sep=""), h=T)
het_sanger = read.table(paste(sangersamples, ".het", sep=""), h=T)
het_sanger_preqc = read.table(paste(sangersamples_preqc, ".het", sep=""), h=T)

het_all=rbind(het,het_cases)

het$meanHet = (het$N.NM. - het$O.HOM.) / het$N.NM.
het_upper=mean(het$meanHet) + (hetlinesd*sd(het_all$meanHet))
het_lower=mean(het$meanHet) - (hetlinesd*sd(het_all$meanHet))

het_cases$meanHet = (het_cases$N.NM. - het_cases$O.HOM.) / het_cases$N.NM.
het_cases_upper=mean(het_cases$meanHet) + (hetlinesd*sd(het_cases$meanHet))
het_cases_lower=mean(het_cases$meanHet) - (hetlinesd*sd(het_cases$meanHet))

het_sanger$meanHet = (het_sanger$N.NM. - het_sanger$O.HOM.) / het_sanger$N.NM.
het_sanger_upper=mean(het_sanger$meanHet) + (hetlinesd*sd(het_sanger$meanHet))
het_sanger_lower=mean(het_sanger$meanHet) - (hetlinesd*sd(het_sanger$meanHet))


het_sanger_preqc$meanHet = (het_sanger_preqc$N.NM. - het_sanger_preqc$O.HOM.) / het_sanger_preqc$N.NM.
het_sanger_preqc_upper=mean(het_sanger_preqc$meanHet) + (hetlinesd*sd(het_sanger_preqc$meanHet))
het_sanger_preqc_lower=mean(het_sanger_preqc$meanHet) - (hetlinesd*sd(het_sanger_preqc$meanHet))

het_all$meanHet = (het_all$N.NM. - het_all$O.HOM.) / het_all$N.NM.
het_all_upper=mean(het_all$meanHet) + (hetlinesd*sd(het_all$meanHet))
het_all_lower=mean(het_all$meanHet) - (hetlinesd*sd(het_all$meanHet))

imiss_all=rbind(imiss,imiss_cases)
imiss_all$logF_MISS = log10(imiss_all[,6])


pdf(outputpathname)
plot(imiss$logF_MISS, het$meanHet, xlim=c(-3,0),
     ylim=c(0,0.5), pch=20, xlab="Proportion of missing genotypes", 
     ylab="Heterozygosity rate", axes=F,col="darkblue")
axis(2, at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), tick=T)
axis(1, at=c(-3,-2,-1,0), labels=c(0.001,0.01,0.1,1))
points(imiss_cases$logF_MISS, het_cases$meanHet,col="darkred",pch=20)
#points(imiss_sanger$logF_MISS, het_sanger$meanHet,col="orange",pch=20)
points(imiss_sanger_preqc$logF_MISS, het_sanger_preqc$meanHet,col="gray",pch=20)

abline(h=mean(het_all$meanHet) - (hetlinesd*sd(het_all$meanHet)), col="sienna1", lty=2)
abline(h=mean(het_all$meanHet) + (hetlinesd*sd(het_all$meanHet)), col="sienna1", lty=2)
abline(h=mean(het_all$meanHet) - ((hetlinesd+1)*sd(het_all$meanHet)), col="sienna2", lty=2)
abline(h=mean(het_all$meanHet) + ((hetlinesd+1)*sd(het_all$meanHet)), col="sienna2", lty=2)
abline(h=mean(het_all$meanHet) - ((hetlinesd+2)*sd(het_all$meanHet)), col="sienna4", lty=2)
abline(h=mean(het_all$meanHet) + ((hetlinesd+2)*sd(het_all$meanHet)), col="sienna4", lty=2)
abline(v=log10(missingnessline), col="sienna2", lty=2)
abline(v=log10(missingnessline2), col="sienna4", lty=2)
dev.off()

#END


