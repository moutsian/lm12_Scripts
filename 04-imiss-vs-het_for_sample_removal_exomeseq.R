inputpathname = "tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.rsID"

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

samples_to_remove=unique(c(which(n_singl[,3]>n_singl_upper),which(het_lowfreq$meanHet>het_lowfreq_upper),which(het_lowfreq$meanHet<het_lowfreq_lower),which(het$meanHet>het_upper),which(het$meanHet<het_lower),which(imiss$F_MISS>missingnessline2)))
write.table(het[samples_to_remove,c(1:2)],"/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/samples_to_remove.qc1.txt",quote=F,col.names=T,row.names=F)
#END