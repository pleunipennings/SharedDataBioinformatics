#Liana's code for chart function

#set directory and open data faile
setwd("~/Desktop/Bioinformatics")
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#examination of data
head(data)
tail(data)
summary(data)
class(data)

#Fig 2: location vs frequency syn/nonsyn/drastic/nonsense.
#   Plot num column vs MeanFreq column in a scatterplot. The points should be colored
#   depending on whether they are syn/nonsyn/drastic/nonsense.

plot(data[,2],data[,3], main = "HIV Practice Data", xlab ="Number", 
ylab ="Mean Frequency", col=data$TypeOfSite)
plot(data)
