inputpathname = "tili.i2.qc"

hetlinesd = 3
hetlinesd2= 4

missingnessline = 0.05
missingnessline2 = 0.1

print(paste("plotting with sd =", hetlinesd, "and miss =s", missingnessline))

imiss = read.table(paste(inputpathname, ".imiss", sep=""), h=T)
imiss$logF_MISS = log10(imiss[,5])


het = read.table(paste(inputpathname, ".het", sep=""), h=T)
het_lowfreq = read.table(paste(inputpathname, ".lowfreq.het", sep=""), h=T)



het$meanHet = (het$N_SITES - het$O.HOM.) / het$N_SITES
het_upper=mean(het$meanHet) + (hetlinesd*sd(het$meanHet))
het_lower=mean(het$meanHet) - (hetlinesd*sd(het$meanHet))

het_lowfreq$meanHet = (het_lowfreq$N_SITES - het_lowfreq$O.HOM.) / het_lowfreq$N_SITES
het_lowfreq_upper=mean(het_lowfreq$meanHet) + (hetlinesd*sd(het_lowfreq$meanHet))
het_lowfreq_lower=mean(het_lowfreq$meanHet) - (hetlinesd*sd(het_lowfreq$meanHet))

n_singl = read.table(paste(inputpathname, ".ifreqburden", sep=""), h=T)
n_singl_upper=mean(n_singl[,3]) + (hetlinesd*sd(n_singl[,3]))
n_singl_lower=mean(n_singl[,3]) - (hetlinesd*sd(n_singl[,3]))

plot(imiss$logF_MISS, het$meanHet, xlim=c(-3,0),
     ylim=c(0,0.1), pch=20, xlab="Proportion of missing genotypes", 
     ylab="Heterozygosity rate, all sites", axes=F,col="darkblue")
axis(2, at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), tick=T)
axis(1, at=c(-3,-2,-1,0), labels=c(0.001,0.01,0.1,1))
#points(imiss_cases$logF_MISS, het_lowfreq$meanHet,col="darkred",pch=20)
#points(imiss_sanger$logF_MISS, het_sanger$meanHet,col="orange",pch=20)
#points(imiss_sanger_preqc$logF_MISS, het_sanger_preqc$meanHet,col="gray",pch=20)

abline(h=mean(het$meanHet) - (hetlinesd*sd(het$meanHet)), col="sienna1", lty=2)
abline(h=mean(het$meanHet) + (hetlinesd*sd(het$meanHet)), col="sienna1", lty=2)
abline(h=mean(het$meanHet) - ((hetlinesd+1)*sd(het$meanHet)), col="sienna2", lty=2)
abline(h=mean(het$meanHet) + ((hetlinesd+1)*sd(het$meanHet)), col="sienna2", lty=2)
abline(h=mean(het$meanHet) - ((hetlinesd+2)*sd(het$meanHet)), col="sienna4", lty=2)
abline(h=mean(het$meanHet) + ((hetlinesd+2)*sd(het$meanHet)), col="sienna4", lty=2)
abline(v=log10(missingnessline), col="sienna2", lty=2)
abline(v=log10(missingnessline2), col="sienna4", lty=2)


#same but now heterozygosity is calculated only using low freq variants (MAF<1%)
plot(imiss$logF_MISS, het_lowfreq$meanHet, xlim=c(-3,0),
     ylim=c(0,0.015), pch=20, xlab="Proportion of missing genotypes", 
     ylab="Heterozygosity rate, only low freq sites", axes=F,col="darkblue")
axis(2, at=c(0,0.01,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), tick=T)
axis(1, at=c(-3,-2,-1,0), labels=c(0.001,0.01,0.1,1))
#points(imiss_cases$logF_MISS, het_lowfreq$meanHet,col="darkred",pch=20)
#points(imiss_sanger$logF_MISS, het_sanger$meanHet,col="orange",pch=20)
#points(imiss_sanger_preqc$logF_MISS, het_sanger_preqc$meanHet,col="gray",pch=20)

abline(h=mean(het_lowfreq$meanHet) - (hetlinesd*sd(het_lowfreq$meanHet)), col="sienna1", lty=2)
abline(h=mean(het_lowfreq$meanHet) + (hetlinesd*sd(het_lowfreq$meanHet)), col="sienna1", lty=2)
abline(h=mean(het_lowfreq$meanHet) - ((hetlinesd+1)*sd(het_lowfreq$meanHet)), col="sienna2", lty=2)
abline(h=mean(het_lowfreq$meanHet) + ((hetlinesd+1)*sd(het_lowfreq$meanHet)), col="sienna2", lty=2)
abline(h=mean(het_lowfreq$meanHet) - ((hetlinesd+2)*sd(het_lowfreq$meanHet)), col="sienna4", lty=2)
abline(h=mean(het_lowfreq$meanHet) + ((hetlinesd+2)*sd(het_lowfreq$meanHet)), col="sienna4", lty=2)
abline(v=log10(missingnessline), col="sienna2", lty=2)
abline(v=log10(missingnessline2), col="sienna4", lty=2)

#now a boxplot for number of singletons
cutoff=2100
n_singl[n_singl[,3]>cutoff,3]=cutoff
which(n_singl[,3]>n_singl_upper)
boxplot(n_singl[,3],ylim=c(0,2100),ylab="Number of singletons per sample")
abline(h=n_singl_upper,col="sienna1",lty=2)

#END


