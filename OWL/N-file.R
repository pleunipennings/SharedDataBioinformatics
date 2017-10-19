install.packages("gdata")

library(gdata)                   # load gdata package 
help(read.xls)                   # documentation 
mydata = read.xls("OverviewSelCoeff_BachelerFilter.csv")  # read from first sheet