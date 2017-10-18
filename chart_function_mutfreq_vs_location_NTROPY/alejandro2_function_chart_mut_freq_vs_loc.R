#Alejandro's code for chart function

#set directory and open data faile
setwd("~/Desktop")

data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")
#HIV practice file read into R as csv 
#examination of data summary,tail,and head 
head(data)
tail(data)
summary(data)
plot(data[,2],data[,3],col=data$TypeOfSite,xlab = "Number",
     ylab = "Mean Frequency")

str(data)
dim(data)
class(data)
which(data$WTAA=="A")

