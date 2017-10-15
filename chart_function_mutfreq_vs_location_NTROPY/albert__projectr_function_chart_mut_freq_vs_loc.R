#Albert's code for chart function

#set directory and open data faile
setwd("C:/Users/alber/Downloads")
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#examination of data
head(data)
tail(data)
summary(data)
class(data)

#Fig 2: location vs frequency syn/nonsyn/drastic/nonsense.
#   Plot num column vs MeanFreq column in a scatterplot. The points should be colored
#   depending on whether they are syn/nonsyn/drastic/nonsense.

