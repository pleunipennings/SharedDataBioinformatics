install.packages("gdata")

library(gdata)                   # load gdata package 
help(read.xls)                   # documentation 
mydata = read.csv("OverviewSelCoeff_BachelerFilter.csv.csv")  # read from first sheet
print(mydata)
