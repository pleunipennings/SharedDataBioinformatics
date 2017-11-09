#This function makes a plot that shows the location of CpG and noCpG sites vs. Frequency.
#This function assumes that you have num, freq, and CPG are already in your dataframe.

#Team: Christen, Rima, Nicole, and Kellen.
#Would also like to thank Scott and Pleuni for their help with making corrections to our code. 
virus = function(VIRUS_DATA) {
#Identify noCpG sites
which(VIRUS_DATA$CPG=="0")->noCpG
print(noCpG)
#Identify CpG sites
which(VIRUS_DATA$CPG=="1")->CpG

#Identify mean frequency in association with noCpG sites
VIRUS_DATA[noCpG,"freq"]->Freq
#Identify mean frequency in association with CpG sites
VIRUS_DATA[CpG,"freq"]->Freq2



#Making a plot of mean frequency vs. CpG/nonCpG location
plot(noCpG, Freq+0.0001, ylim=c(0.0001,0.5), col="black",pch=21.25, bg=rgb(1,0,0,0.5), main="CpG/noCpG location vs. Mean Frequency", xlab = "Location", ylab = "Mean Frequency")
 
#Overlapping the two graphs
points(CpG, Freq2+0.0001, col="black",pch=21.25, bg=rgb(0,0,1,0.5))

}
virus(VIRUS_DATA)

#legend with circle instead of line  
legend("topleft",c("noCpG","CpG"),cex=0.6, 
      col ="black",bg=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),bty="n",
       title="Legend",inset=.02,pch =21,pt.bg = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))


