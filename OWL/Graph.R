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
plot(noCpG, Freq+0.0001, ylim=c(0.0001,0.5), log = "y", col="black",pch=21.25, bg=rgb(1,0,0,0.5), main="CpG/noCpG location vs. Mean Frequency", xlab = "Location", ylab = "Mean Frequency")
 
#Overlapping the two graphs
points(CpG, Freq2+0.0001, col="black",pch=21.25, bg=rgb(0,0,1,0.5),xaxt='n',yaxt='n', ann = FALSE)



#legend with circle instead of line  
legend("topleft",c("noCpG","CpG"),cex=0.6, 
      col ="black",bg=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),bty="n",
       title="Legend",inset=.02,pch =21,pt.bg = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))
