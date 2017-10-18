#Albert's code for chart function

#set directory and open data faile
setwd("C:/Users/alber/Downloads")
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#examination of data
head(data)
tail(data)
summary(data)
class(data)
list(data)
str(data)

#Fig 2: location vs frequency syn/nonsyn/drastic/nonsense.
#   Plot num column vs MeanFreq column in a scatterplot. The points should be colored
#   depending on whether they are syn/nonsyn/drastic/nonsense.

#removing overlap, res from data NOT SUCCESSFUL

summary(data)
summary(data$TypeOfSite != c("overlap","res"))
summary(data[data$TypeOfSite != c("overlap","res"), ,drop=FALSE)])
summary(data[ data$TypeOfSite != c("overlap","res") , ,drop=FALSE]
)
subset(data,TypeOfSite!="overlap" | "res"))

#plot of data with legend

plot(data[,2],
     data[,3],
     col=data$TypeOfSite,
     pch=c(data$TypeOfSite)
)

#add legend but not yet right
legend("topright", 
       legend = levels(data$TypeOfSite), 
       pch=c(1,2,3,4,5),
       col=c(1,2,3,4,5)
)
str(data$TypeOfSite)
data$TypeOfSite
View(data$TypeOfSite)
