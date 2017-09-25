#Prepare cumulative plots by concentration from the gene counts, as per Carl's request.

data=read.table("C:/Academic/SANGER/PAXGENE/ISR.transcript_to_gene.counts",head=T) #note that these are only ISR counts.
data=data[,c(13,14,15,16,8,9,2,4,11,10,3,12,6,7,5,1)] # ranking the data. Make sure this holds (as it depends on the ranking of columns in the .counts file)


#For individuals 
COLS=c("black","darkblue","blue","lightblue")
START=5
plot(cumsum(data[,START]/sum(data[,START])),col=COLS[1],main=list("cumulative reads as a percentage of total for each sample - individual 2",cex=1.1))
points(cumsum(data[,START+1]/sum(data[,START+1])),col=COLS[2])
points(cumsum(data[,START+2]/sum(data[,START+2])),col=COLS[3])
points(cumsum(data[,START+3]/sum(data[,START+3])),col=COLS[4])
legend(28000,0.2,col=COLS,legend=c("2.5R","2.5M","1M","0.5M"),pch=19,cex=0.8)
#END


#For individuals - not as a percentage
COLS=c("black","darkblue","blue","lightblue")
START=1
plot(cumsum(data[,START]),col=COLS[1],ylim=c(0,max(c(sum(data[,START]),sum(data[,START+1]),sum(data[,START+2]),sum(data[,START+2])))),main=list("cumulative reads for each sample - individual 1",cex=1.1))
points(cumsum(data[,START+1]),col=COLS[2])
points(cumsum(data[,START+2]),col=COLS[3])
points(cumsum(data[,START+3]),col=COLS[4])
legend(28000,3500000,col=COLS,legend=c("2.5R","2.5M","1M","0.5M"),pch=19,cex=0.8)
#END

#For individuals - not as a percentage
COLS=c("black","darkblue","blue","lightblue")
START=1
MAX=200
plot(cumsum(data[,START])[1:MAX],col=COLS[1],ylim=c(0,180000),main=list("cumulative reads for each sample - individual 1",cex=1.1))
points(cumsum(data[,START+1])[1:MAX],col=COLS[2])
points(cumsum(data[,START+2])[1:MAX],col=COLS[3])
points(cumsum(data[,START+3])[1:MAX],col=COLS[4])
legend(28000,3500000,col=COLS,legend=c("2.5R","2.5M","1M","0.5M"),pch=19,cex=0.8)
#END
