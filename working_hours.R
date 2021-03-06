# Usage of folder by Day and Time
# Data in format: Month Day Time

data=read.table("C:/Academic/SANGER/test.data",stringsAsFactors=F)

#Assuming 2017 as year for first implementation
MM=match(data[,1],month.abb)
Date=paste("2017-",MM,"-",data[,2],sep="")
Days=weekdays(as.Date(Date))
DaysN=as.POSIXlt(Date)$wday #this gives numbers, with one being Monday.

par(mfrow=c(1,2))
# files by weekday
hist(DaysN,breaks=15,xlim=c(0.5,5.5),xaxt='n',col="darkblue",xlab="Days of the week", main=list("No of files generated by weekday",cex=1.1))
axis(side=1, at=1:5, labels=c("Mon","Tue","Wed","Thu","Fri"))

# files generated by time
Time=sapply(strsplit(data[,3],":"),
  function(x) {
    x <- as.numeric(x)
    x[1]+x[2]/60
    }
)
hist(Time,col="gray",main=list("No of files generated by time of the day",cex=1.1),breaks=15)


#END