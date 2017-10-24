#script for Team NTROPY's function graphing Frequency of synonymous/nonsynonmous/other mutation vs location

#set directory and open data file
setwd("C:/Users/alber/Downloads")

#read csv file into variable "data"
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#plot data mut freq vs location. 
#   Assumptions: plot based on column names, "num", "MeanFreq", "TypeOfSite".
#   Your data must be in these columns! 
plot(
    #x vector
    data$num,
    #y vector
    log(data$MeanFreq),
    #make black empty circles as symbol
    pch=21,
    #make outline of symbol black
    col= "black",
    #fill inside of point with color by factor category "TypeOfSite" bg=
    bg=data$TypeOfSite,
    #Title label
    main = "HIV Practice Data",
    #x axis label
    xlab ="Location on Sequence", 
    #y axis label
    ylab ="Mean Frequency of Mutation",
    cex=3
)

#add legend in top right corner
legend("topright", 
       #names of each category
       legend = levels(data$TypeOfSite), 
       #symbols matching dataframe's factors 1:5 of data$TypeOfSite
       pch=21,
       #colors matching dataframe's factors 1:5 of data$TypeOfSite
       col="black",
       pt.bg=c(1:5)
)


