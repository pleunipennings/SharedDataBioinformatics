#### Master Script ####

#TEST 

#Fig 2: location vs frequency CpG non-CpG.
  #Plot num column vs MeanFreq column in a scatterplot. The points should be colored
  #depending on whether they are Cpg / not CpG.

#remember [PULL > Comment > Commit > PUSH]




#### function to find CpG Islands ####

#setwd(SharedDataBioinformatics/11AM_Influ_team)    
seq_CG <- read.fasta("class25Influ.txt")
        x= split(seq, (0:nrow(seq_CG) %/% 500))

    #slide function scrolls through variable set looking for matching pairs then outputs T/F if pair for the phrame is found.
CpG_Finder <- function(data, window, step){
    a=lapply(x, function(vec){
        x <- gregexpr("gc", vec, perl = TRUE)
        res <- sum(attr(x[[1]], "match.length"))
        res
    })
    b=lapply(x, function(vec){
        x <- gregexpr("g", vec, perl = TRUE)
        res <- regmatches(vec, x)
        res
    }) 
    c=lapply(x, function(vec){
        x <- gregexpr("c", vec, perl = TRUE)
        res <- regmatches(vec, x)
        res
    }) 
    total <- length(data)
    spots <- seq(from = 1, to = (total-window), length.out = step)
    result <- vector(length = length(spots))
    for(i in 1:length(spots)){result[i]}
    
    return(result)
}
CpG_Finder(seq_CG,4,500)

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





