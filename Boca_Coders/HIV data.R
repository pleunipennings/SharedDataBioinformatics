setwd("~/Downloads")
read.csv("HIV.csv")
HIVdata = read.csv("../Downloads/HIV.csv")
data.frame(HIVdata$num, HIVdata$WTnt, HIVdata$MeanFreq) -> HIVdataframe
summary(HIVdata)
#pull out 3 column B,C and G 2,3, & 7
plot(HIVdata$num, HIVdata$MeanFreq, log="y")
     

