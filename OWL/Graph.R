read.csv("OverviewSelCoeff_BachelerFilter.csv")->HIVdata
print(HIVdata)


#Identify noCpG sites
#column might be 13 of 15, so change accordingly 
which(HIVdata[,15]=="0")->noCpG
print(noCpG)
#Identify CpG sites
which(HIVdata[,15]=="1")->CpG

#Identify mean frequency in association with noCpG sites
HIVdata[noCpG,3]->Freq
#Identify mean frequency in association with CpG sites
HIVdata[CpG,3]->Freq2


#Making a plot of mean frequency vs. CpG/nonCpG location
plot(noCpG, Freq, col="red", main="Mean Frequency vs. CpG/noCpG location", xlab = "Location", ylab = "Mean Frequency")
#Overlapping the two graphs
par(new = TRUE )
plot(CpG, Freq2,col="blue",xaxt='n',yaxt='n', ann = FALSE)


legend("topright",c("noCpG","CpG"),lty = c(1,1),lwd=4, col=c("red","blue"),bty="n")



