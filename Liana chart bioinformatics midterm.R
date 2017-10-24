#script for Team NTROPY's function graphing Frequency of synonymous/nonsynonmous/other mutation vs location

#set directory and open data file
setwd("~/Desktop/Bioinformatics")

#read csv file into variable "data"
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#plot data mut freq vs location. 
#   Assumptions: plot based on column names, "num", "MeanFreq", "TypeOfSite".
#   Your data must be in these columns! 

plot(
  #x vector 
  data$num,
  #y vector
  data$MeanFreq,
  log = "y",
  #colors by category factors from column "TypeOfSite
  col="black",
  #symbols attached by category factors from column "TypeOfSite"
  pch=21,
  bg = data$TypeOfSite,
  #Title label
  main = "HIV Practice Data",
  #x axis label
  xlab ="Location on Sequence", 
  #y axis label
  ylab ="Mean Frequency of Mutation", 
  ylim= c(.001,.1000), 
  cex = 3
)
  
  #set y limits to expand data
x

#add legend in top right corner
legend("topright", 
       #names of each category
       legend = levels(data$TypeOfSite), 
       #symbols matching dataframe's factors of data$TypeOfSite
       pch=21,
       #colors matching dataframe's factors of data$TypeOfSite
       col="black",
       pt.bg= c(1:5)
)
