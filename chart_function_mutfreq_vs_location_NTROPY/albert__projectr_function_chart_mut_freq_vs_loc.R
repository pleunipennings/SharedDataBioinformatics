#Albert's code for chart function

#set directory and open data file
setwd("C:/Users/alber/Downloads")

#read csv file into variable "data"
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#examination of data
head(data)
tail(data)
summary(data)
class(data)
list(data)
#structure
str(data)

#Fig 2: location vs frequency syn/nonsyn/drastic/nonsense.
#   Plot num column vs MeanFreq column in a scatterplot. The points should be colored
#   depending on whether they are syn/nonsyn/drastic/nonsense.

#plot of data with legend

plot(data[,2],
     data[,3],
     col=data$TypeOfSite,
     pch=c(data$TypeOfSite)
)

#add legend in top right corner
legend("topright", 
       #names of each category
       legend = levels(data$TypeOfSite), 
       #symbols matching dataframe's factors of data$TypeOfSite
       pch=c(1,2,3,4,5),
       #colors matching dataframe's factors of data$TypeOfSite
       col=c(1,2,3,4,5)
)

