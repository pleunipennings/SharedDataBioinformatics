#script for Team NTROPY's function graphing Frequency of synonymous/nonsynonmous/other mutation vs location

#set directory and open data file
setwd("C:/Users/alber/Downloads")

#read csv file into variable "data"
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#plot data mut freq vs location
plot(data[,2],
     data[,3],
     col=data$TypeOfSite,
     pch=c(data$TypeOfSite),
     #Title label
    main = "HIV Practice Data",
    #x axis label
    xlab ="Location on Sequence", 
    #y axis label
    ylab ="Mean Frequency of Mutation"
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
