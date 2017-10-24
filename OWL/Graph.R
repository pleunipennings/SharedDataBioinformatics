read.csv("OverviewSelCoeff_BachelerFilter.csv")->HIVdata
print(HIVdata)


#Identify noCpG sites
which(HIVdata$makesCpG=="0")->noCpG
print(noCpG)
#Identify CpG sites
which(HIVdata$makesCpG=="1")->CpG

#Identify mean frequency in association with noCpG sites
HIVdata[noCpG,"MeanFreq"]->Freq
#Identify mean frequency in association with CpG sites
HIVdata[CpG,"MeanFreq"]->Freq2



#Making a plot of mean frequency vs. CpG/nonCpG location
plot(noCpG, Freq, log="y",col="black",pch=21.25, bg="red", main="CpG/noCpG location vs. Mean Frequency", xlab = "Location", ylab = "Mean Frequency")
#Overlapping the two graphs
par(new = TRUE )
plot(CpG, Freq2, col="black",pch=21.25, bg="blue",xaxt='n',yaxt='n', ann = FALSE)


legend("topright",c("noCpG","CpG"),lty = c(1,1),lwd=4, col=c("red","blue"),bty="n")



