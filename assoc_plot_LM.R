# 
#
#              Diabetes Genetics Initiative of Broad Institute of Harvard and MIT, Lund University and 
#                                  Novartis Institutes of BioMedical Research
#        Whole-genome association analysis identifies novel loci for type 2 diabetes and triglyceride levels
#                             Science 2007 Jun 1;316(5829):1331-6. Epub 2007 Apr 26 
#

# File Modified by Loukas #####
#setwd("U:/R_workspace")
setwd("T:/MomokoPlot/")

make.fancy.locus.plot <- function(snp, locusname, chr, locus, range, best.pval) {
#example: 
# locus <- read.table("U:/data/DGI_chr9_hit.txt", header=T, row.names=1)
# make.fancy.locus.plot("rs10811661", "CDKN2A/CDKN2B region", "9", locus, 9, 5.4e-8)

hit <- locus[snp,]

#
# size of the region
#
min.pos <- min(locus$POS) - 10000
max.pos <- max(locus$POS) + 10000
size.pos <- max.pos - min.pos
center.pos <- min.pos + ( size.pos / 2 )
center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

#
# range of y-axis
#
# this dedicates 33% of the yaxis to the genes, labels, recomb rate
offset <- ( range * 4 / 3 ) - range
big.range <- range + offset 

ystart.gene <- - offset
ystart.recomb <- - offset + (big.range / 6)


#
# recombination rate (b37)
# example: genetic_map_chr9_combined_b37.txt
#recomb <- read.table(paste("../data/genetic_map_chr", chr, "_combined_b37.txt", sep=""), header=T)
recomb <- read.table(paste("T:/MomokoPlot/GeneticMaps/genetic_map_chr", chr, "_combined_b37.txt", sep=""), header=T)

keep.recomb <- subset(recomb, recomb[,1] > min.pos & recomb[,1] < max.pos)

#
# genes in the region
#
genelist <- read.table(paste("T:/MomokoPlot/Genes/known_genes_build37_chr", chr, ".transformed.txt", sep=""), header=T)
#genelist <- read.table(paste("../data/known_genes_build37_chr", chr, ".transformed.txt", sep=""), header=T)
genes.in.locus <- subset(genelist, ( genelist$START > min.pos & genelist$START < max.pos ) | ( genelist$STOP > min.pos & genelist$STOP < max.pos) )
print(genes.in.locus)

#
# genotyped markers
#
markers.in.strong.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.8 & locus$TYPE == "typed"))
markers.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == "typed"))
markers.in.weak.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == "typed"))
markers.not.in.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR<0.2 & locus$TYPE == "typed"))
#markers.typed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "typed"))

#
# imputed SNPs
#
imputed.in.strong.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.8 & locus$TYPE == "imputed"))
imputed.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == "imputed"))
imputed.in.weak.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == "imputed"))
imputed.not.in.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR<0.2 & locus$TYPE == "imputed"))
#markers.imputed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "imputed"))


par(mar=c(4,4,3,4))

#
# start plot with recombination rate (in background)
#
plot(keep.recomb[,1], ystart.recomb + ( ( keep.recomb[,2] / 60 ) * ( 6 * big.range / 8 )), type="l", col="lightblue", lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)


#
# axes, titles and legends
#
mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
axis(1, at=c(center.100kb.pos - offset.100kb.pos, center.100kb.pos, center.100kb.pos + offset.100kb.pos), labels=c((center.100kb.pos - offset.100kb.pos) / 1000, center.100kb.pos / 1000, (center.100kb.pos + offset.100kb.pos) / 1000), las=1) 

axis(2, at=seq(0,range,2), labels=seq(0,range,2), las=1) 
mtext("Observed (-logP)", side=2, at=(range/2), line=2)

axis(4, at=c( ystart.recomb, ystart.recomb + (big.range / 4), ystart.recomb + ( 2 * big.range / 4), ystart.recomb + ( 3 * big.range / 4 ) ), labels=c("0","20","40","60"), las=1)
mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range/2), line=2)


box()
lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")

#
# this is the hit
#
points(hit$POS, -(log10(hit$PVAL)), pch=23, cex=2.5, bg="red")
text(hit$POS, -(log10(hit$PVAL)), labels=c(row.names(hit)), pos=3, offset=1)
if ( -(log10(best.pval)) < range ) {
	points(hit$POS, -(log10(best.pval)), pch=23, cex=2.5, bg="blue")
	text(hit$POS, -(log10(best.pval)), labels=c(paste("P=",best.pval,sep="")), pos=4, offset=2)
}else {
	points(hit$POS, range, pch=23, cex=2.5, bg="blue")
	text(hit$POS, range, labels=c(paste("P=",best.pval,sep="")), pos=4, offset=1)
}


#
# plot the genotyped markers
#
points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), pch=23, cex=1.0, bg="white")
points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), pch=23, cex=1.25, bg="yellow")
points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)), pch=23, cex=1.25, bg="orange")
points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), pch=23, cex=1.25, bg="red")
#points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), pch=23, cex=1.0, bg=hsv(0,markers.not.in.ld$RSQR,1))
#points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), pch=23, cex=1.25, bg=hsv(0,markers.in.weak.ld$RSQR,1))
#points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)), pch=23, cex=1.25, bg=hsv(0,markers.in.moderate.ld$RSQR,1))
#points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), pch=23, cex=1.25, bg=hsv(0,markers.in.strong.ld$RSQR,1))

#
# plot the imputed SNPs
#
#points(imputed.not.in.ld$POS, -(log10(imputed.not.in.ld$PVAL)), pch=21, cex=1.0, bg="white")
#points(imputed.in.weak.ld$POS, -(log10(imputed.in.weak.ld$PVAL)), pch=21, cex=1.0, bg="grey")
#points(imputed.in.moderate.ld$POS, -(log10(imputed.in.moderate.ld$PVAL)), pch=21, cex=1.5, bg="orange")
#points(imputed.in.strong.ld$POS, -(log10(imputed.in.strong.ld$PVAL)), pch=21, cex=1.8, bg="red")
points(imputed.not.in.ld$POS, -(log10(imputed.not.in.ld$PVAL)), pch=23, cex=1.0, bg="grey")
points(imputed.in.weak.ld$POS, -(log10(imputed.in.weak.ld$PVAL)), pch=23, cex=1.0, bg="grey")
points(imputed.in.moderate.ld$POS, -(log10(imputed.in.moderate.ld$PVAL)), pch=23, cex=1.0, bg="grey")
points(imputed.in.strong.ld$POS, -(log10(imputed.in.strong.ld$PVAL)), pch=23, cex=1.0, bg="grey")

#
# plot the genes
#
offset=offset+0.7
add_to_offset=0.5
cnt= 0.5
vpos = matrix(ncol=5,nrow=1,0)
for ( i in 1:nrow(genes.in.locus) ) { 
	if ( genes.in.locus[i,]$STRAND == "+" ) {
		if (i>1){
			j=1
			while(j<i & j<min(which(vpos==0)) ){
				if(genes.in.locus[vpos[j],]$STOP+25000<genes.in.locus[i,]$START){					
					break;
				}else{
					j=j+1
				}
			}
			vpos[j]=i;
			add_to_offset=0.5*j
			arrows(max(genes.in.locus[i,]$START, min.pos), -offset+add_to_offset, min(genes.in.locus[i,]$STOP, max.pos), -offset+add_to_offset, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")	
		}else{
			arrows(max(genes.in.locus[i,]$START, min.pos), -offset+0.5, min(genes.in.locus[i,]$STOP, max.pos), -offset+0.5, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")	
			vpos[i]=i;
		}
	}else {
		if (i>1){
			j=1
			while(j<i & j<min(which(vpos==0)) ){
				if(genes.in.locus[vpos[j],]$STOP+25000<genes.in.locus[i,]$START){
					break;
				}else{
					j=j+1
				}
			}
			vpos[j]=i;
			add_to_offset=0.5*j
			arrows(max(genes.in.locus[i,]$START, min.pos), -offset+add_to_offset, min(genes.in.locus[i,]$STOP, max.pos), -offset+add_to_offset, length=0.05, lwd=2, code=1, lty="solid", col="chocolate")			
		}else{
			arrows(max(genes.in.locus[i,]$START, min.pos), -offset+0.5, min(genes.in.locus[i,]$STOP, max.pos), -offset+0.5, length=0.05, lwd=2, code=1, lty="solid", col="chocolate")			
			vpos[i]=i;
		}
	}
	if ( ! is.na(genes.in.locus[i,]$GENE) ) {
		text(genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE / 2), -offset+add_to_offset + ( big.range / 50 ), labels=genes.in.locus[i,]$GENE, cex=0.5)
	}
}

}




#locus <- read.table("U:/data/DGI_chr9_hit.txt", header=T, row.names=1)
locus <- read.table("T:/MomokoPlot/Data/DGI_chr9_hit.txt", header=T, row.names=1)


#pdf("U:/data/assocplot.pdf", width=6, height=4)
pdf("T:/MomokoPlot/Data/assocplot.pdf", width=6, height=4)


make.fancy.locus.plot("rs10811661", "CDKN2A/CDKN2B region", "9", locus, 9, 5.4e-8)

dev.off()



#quit()

