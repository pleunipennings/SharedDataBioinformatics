#### Master Script ####

#TEST 

#Fig 2: location vs frequency CpG non-CpG.
  #Plot num column vs MeanFreq column in a scatterplot. The points should be colored
  #depending on whether they are Cpg / not CpG.

#remember [PULL > Comment > Commit > PUSH]

#### Fuction to print Location vs frequency CpG non-CpG Graph ####

#setwd(SharedDataBioinformatics/11AM_Influ_team)   
#following packages are required >>
library(graphics)
library(seqinr)
# HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
# n <- data.frame(seqinr::consensus(HPIV1a))

n <- data.frame(read.csv("OverviewSelCoeff_BachelerFilter.csv"))
  
str(n)
head(n)

which(n$makesCpG=="1")

YCpG <- which(n$makesCpG=="1")
YCpG #lists which variables return a "1" these make a CpG island when mutated "yes cpg or Y"

NCpG <- which(n$makesCpG=="0")
NCpG #lists which variables return a "0" these do not make a CpG island when mutated "No cpg or N"

V1 <- n$MeanFreq[YCpG] #create a Value that looks at mean frequency by yes cpg
V2 <- n$MeanFreq[NCpG] #create a Value that looks at mean frequency by no cpg

plot.default(x = c(V1, V2), #plot it!
     xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
     col =  c("blue","red")
)

# CpG graphing function()

LvsF_CpG_Printer <- function(data_frame){
  YCpG <- which(n$makesCpG=="1")
  NCpG <- which(n$makesCpG=="0")
  x1 <- n$MeanFreq[YCpG]
  x2 <- n$MeanFreq[NCpG]
  
  return(plot.default(x = c(x1, x2), 
                      xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
                      col =  c("blue", "red")
                      )#close plot.default
         )#close return
}#close function

LvsF_CpG_Printer(n) #run function

#notes for fucntion:
# n = must be a data.frame
#

#### function to find CpG Islands (NOTICE: WE DO NOT HAVE TO FIND THEM) ####

#setwd(SharedDataBioinformatics/11AM_Influ_team)    
seq_CG <- read.alignment("class25Influ.txt", format = "fasta")
seqinr::consensus(seq_CG) 
nuc <- data.frame(x = (seqinr::consensus(seq_CG)))
as.matrix(nuc)
nuc$x
?seq.default
#slide function scrolls through variable set looking for matching pairs then outputs T/F if pair for the phrame is found.
CpG_Finder <- function(D, window, deslength){
  total <- length(D)
  lan <- deslength
  data <- D$x
  aa=lapply(x, function(data){
        x <- gregexpr("gc", data, perl = TRUE)
        res1 <- sum(attr(x[[1]], "match.length"))
        res1
    })
  bb=lapply(x, function(data){
        x <- gregexpr("g", data, perl = TRUE)
        res2 <- regmatches(data, x)
        res2
    })
  cc=lapply(x, function(data){
        x <- gregexpr("c", data, perl = TRUE)
        res3 <- regmatches(data, x)
        res3
    })
  
  result <- vector(length = length(data))
    for(i in 1:length(data)){result[i]}
    
  
    return(result)
}
CpG_Finder(nuc,2,500)


#### notes and comments section ####
#10-17-2017 nathan 12:53
    #Formatted the group master script
    # added sections
# pure Fasta format Upload
HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIV1b = read.alignment("humanparainfluenzavirus1_F.fasta_pruned.mu.trim05", format = "fasta")
HPIV1c = read.alignment("humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05", format = "fasta")
HPIV3a = read.alignment("humanparainfluenzavirus3.fasta_pruned.mu.trim05", format = "fasta")
HPIV3b = read.alignment("humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05", format = "fasta")
#average size of each sample is 15,500~ approx
# none of this worked!
#       str(HPIV1a)
#       as.table(HPIV1a, stringsAsFactors=FALSE)
#       M <- as.table(read.alignment("clean_HPIV1.txt", format = "fasta"))
#       M

#       B <-read.csv("Clean_HPIV1.csv")
#       B
#       str(B)

# I took each fasta file and converteded it into a csv tabulated format
# the T stands for Tabulated format or .csv
HPIV1aT = read.csv("HPIV1aT.csv") #from HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIV1bT = read.csv("HPIV1bT.csv") #from HPIV1b = read.alignment("humanparainfluenzavirus1_F.fasta_pruned.mu.trim05", format = "fasta")
HPIV1cT = read.csv("HPIV1cT.csv") #from HPIV1c = read.alignment("humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05", format = "fasta")
HPIV3aT = read.csv("HPIV3aT.csv") #from HPIV3a = read.alignment("humanparainfluenzavirus3.fasta_pruned.mu.trim05", format = "fasta")
HPIV3bT = read.csv("HPIV3bT.csv") #from HPIV3b = read.alignment("humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05", format = "fasta")

#average size of each sample is 15,500~ approx
str(HPIV1aT)
head(HPIV1aT)





