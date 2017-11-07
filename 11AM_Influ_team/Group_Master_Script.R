#### Master Script ####

#TEST 
#----------------------------------------------
#DUE DATES: 
# functxn's done      12:45 Thursday 10/26/2017
# Posters done        11:10 Tuesday 10/31/2017

#----------------------------------------------
#Fig 2: location vs frequency CpG non-CpG.
  #Plot num column vs MeanFreq column in a scatterplot. The points should be colored
  #depending on whether they are Cpg / not CpG.

#remember [PULL > Comment > Commit > PUSH]

#### Fuction to print Location vs frequency CpG non-CpG Graph ####

#setwd(SharedDataBioinformatics/11AM_Influ_team)   
#following packages are required >>
library(graphics)
library(seqinr)

# value "n" will represent our data.frame of use
n <- data.frame(read.csv("OverviewSelCoeff_BachelerFilter.csv"))

# look at n, how does the data.frame pan out. i.e. what are the col vs row within. 
str(n)
head(n)


which(n$makesCpG=="1")
which(n$makesCpG=="0")

YCpG <- which(n$makesCpG=="1")
YCpG #lists which variables return a "1" these make a CpG island when mutated "yes cpg or Y"

NCpG <- which(n$makesCpG=="0")
NCpG #lists which variables return a "0" these do not make a CpG island when mutated "No cpg or N"

V1 <- n$MeanFreq[YCpG] #create a Value that looks at mean frequency by yes cpg
V2 <- n$MeanFreq[NCpG] #create a Value that looks at mean frequency by no cpg

plot.default(x = c(V1, V2), #plot it!
     xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
     col = (as.integer(n$makesCpG)),
     log = "y"
)

####### CpG graphing function() #########

<<<<<<< HEAD
LvsF_CpG_Printer <- function(data_frame){
if (T) {n$makesCpG <- n$makesCpG+2}
=======
LvsF_CpG_Printer <- function(n){
if (T) {n$makesCpG <- n$makesCpG+1}
>>>>>>> 241c824fd750a2b9ba859d514bf83ce33750302a
  YCpG <- which(n$makesCpG=="2")
  NCpG <- which(n$makesCpG=="1")
  x1 <- n$MeanFreq[YCpG]
  x2 <- n$MeanFreq[NCpG]
  plot.default(x = c(x1, x2), 
               xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
               col = (as.integer(n$makesCpG)),
               log = "y"
  )#close plot.default
  
  legend("topright",legend=levels(n$makesCpG), inset= 0)
  
  return()#close return
}#close function

LvsF_CpG_Printer(n) #run function
head(n)
#notes for fucntion:
# n = must be a data.frame
# for our data:
#       HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
#       n <- data.frame(seqinr::consensus(HPIV1a))
#----------------------------------------------------------------------------------------------
#### section for the use of our data ####
#pure Fasta format Upload
HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIV1b = read.alignment("humanparainfluenzavirus1_F.fasta_pruned.mu.trim05", format = "fasta")
HPIV1c = read.alignment("humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05", format = "fasta")
HPIV3a = read.alignment("humanparainfluenzavirus3.fasta_pruned.mu.trim05", format = "fasta")
HPIV3b = read.alignment("humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05", format = "fasta")
#average size of each sample is 15,500~ approx

# I took each fasta file and converteded it into a csv tabulated format
# the T stands for Tabulated format or .csv
HPIV1aT = read.csv("HPIV1aT.csv") #from HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIV1bT = read.csv("HPIV1bT.csv") #from HPIV1b = read.alignment("humanparainfluenzavirus1_F.fasta_pruned.mu.trim05", format = "fasta")
HPIV1cT = read.csv("HPIV1cT.csv") #from HPIV1c = read.alignment("humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05", format = "fasta")
HPIV3aT = read.csv("HPIV3aT.csv") #from HPIV3a = read.alignment("humanparainfluenzavirus3.fasta_pruned.mu.trim05", format = "fasta")
HPIV3bT = read.csv("HPIV3bT.csv") #from HPIV3b = read.alignment("humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05", format = "fasta")
#       str(HPIV1a)
#       as.table(HPIV1a, stringsAsFactors=FALSE)
#       M <- as.table(read.alignment("clean_HPIV1.txt", format = "fasta"))
#       M

#       B <-read.csv("Clean_HPIV1.csv")
#       B
#       str(B)
str(HPIV1aT)
head(HPIV1aT)

n <- data.frame(Pos = c(1:length(seqinr::consensus(HPIV1a))),
                WTnt = (seqinr::consensus(HPIV1a)),
                Trans = c(1)
)

head(n)
tail(n)


LvsF_CpG_Printer <- function(data_frame){
  YCpG <- which(n$makesCpG=="1")
  NCpG <- which(n$makesCpG=="0")
  x1 <- n$MeanFreq[YCpG]
  x2 <- n$MeanFreq[NCpG]
  
  return(plot.default(x = c(x1, x2), 
                      xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
                      col =  c("blue", "red")

#### function to find CpG Islands (NOTICE: WE DO NOT HAVE TO FIND THEM) ####
#function pruned and sent to Nathan.R






