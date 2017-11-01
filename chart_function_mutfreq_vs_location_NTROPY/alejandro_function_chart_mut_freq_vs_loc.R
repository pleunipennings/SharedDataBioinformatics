#Alejandro's code for chart function

#set directory and open data faile
setwd("~/desktop")
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

plot(
  #x vector
  data$num,
  #y vector
  log(data$MeanFreq),
  #colors by category factors from column "TypeOfSite"
  col=data$TypeOfSite,
  #symbols attached by category factors from column "TypeOfSite"
  pch=20,
  #Title label
  main = "HIV Practice Data",
  #x axis label
  xlab ="Location on Sequence", 
  #y axis label
  ylab ="Mean Frequency of Mutation",
  <<<<<<< HEAD
  #set y limits to expand data
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  #grid superimposes grid onto plot 
  #nx and ny describes x and y axis 
  #NA will automatically set x-axis to default plot ticks 
  
  =======
    cex=3
  >>>>>>> c3c8d1c599cf5a611cdc4de1d1ecb160612a46a6
)

#add legend in top right corner
legend("topright", 
       #names of each category
       legend = levels(data$TypeOfSite), 
       #symbols matching dataframe's factors 1:5 of data$TypeOfSite
       pch=20,
       #colors matching dataframe's factors 1:5 of data$TypeOfSite
       col=c(1:5)
)


#examination of data
head(data)
tail(data)
summary(data)

