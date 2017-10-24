setwd("~/Desktop")
read.csv("OverviewSelCoeff_BachelerFilter.csv")->HIVdata
print(HIVdata)

summary(HIVdata)

which(HIVdata[,4]=="syn")->Synonymous 
which(HIVdata[,4]== "nonsyn")->Nonsynonymous


which(HIVdata[,13]=="0")->noCpG
print(noCpG)
which(HIVdata[,13]=="1")->CpG
print(CpG)

HIVdata[noCpG,7]->loc
HIVdata[CpG,7]->loc2
HIVdata[noCpG,3]->Freq
HIVdata[CpG,3]->Freq2
print(Freq)

print(Freq2)



plot(HIVdata$num, HIVdata$MeanFreq, col=HIVdata$TypeOfSite, pch=c(HIVdata$TypeOfSite), main="HIV", xlab="Location",ylab="Mean")

plot(noCpG, Freq, col="red")
par(new = TRUE )
plot(CpG, Freq2,col="blue")



