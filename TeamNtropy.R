#Team NTROPY. Albert, Liana, Alejandro's project

#function to : To create function that will chart
#   the location of the mutation (x-axis) vs frequency of mutations (y-axis). 

#set directory to where "OverviewSelfCoeff_BachelerFilter.csv" is
setwd("c:/Users/alber/Downloads")

#read file
file<-read.csv("OverviewSelCoeff_BachelerFilter(1).csv")

head(file)
tail(file)
summary(file)
