setwd("~/Desktop")
read.csv("OverviewSelCoeff_BachelerFilter.csv")->HIVdata
print(HIVdata)

#Identify noCpG sites
which(HIVdata[,13]=="0")->noCpG
#Identify CpG sites
which(HIVdata[,13]=="1")->CpG

#Identify mean frequency in association with noCpG sites
HIVdata[noCpG,3]->Freq
#Identify mean frequency in association with CpG sites
HIVdata[CpG,3]->Freq2



#Making a plot of mean frequency vs. CpG/nonCpG location
plot(noCpG, Freq, col="red", main="Mean Frequency vs. CpG/noCpG location", xlab = "Location", ylab = "Mean Frequency")
#Overlapping the two graphs
par(new = TRUE )
plot(CpG, Freq2,col="blue",xaxt='n',yaxt='n', ann = FALSE)



