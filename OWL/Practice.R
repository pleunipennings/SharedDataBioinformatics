setwd("~/Desktop")
read.csv("OverviewSelCoeff_BachelerFilter.csv")->HIVdata
print(HIVdata)
summary(HIVdata)

which(HIVdata[,4]=="syn")->Synonymous 
which(HIVdata[,4]== "nonsyn")->Nonsynonymous

which(HIVdata[,13]=="0")->noCpG
which(HIVdata[,13]=="1")->CpG

HIVdata[noCpG,7]->loc
HIVdata[CpG,7]->loc2







