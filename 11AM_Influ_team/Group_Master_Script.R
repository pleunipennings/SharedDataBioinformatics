#### Master Script ####

#TEST 

#Fig 2: location vs frequency CpG non-CpG.
  #Plot num column vs MeanFreq column in a scatterplot. The points should be colored
  #depending on whether they are Cpg / not CpG.

#remember [PULL > Comment > Commit > PUSH]




#### function to find CpG Islands ####

CpG_finder <- function(data, window, step){
        total <- length(data)
        spots <- seq(from = 1, to = (total-window), by = step)
        result <- vector(length = length(spots))
        for(i in 1:length(spots)){result[i]}
        return(result)
    }
slide_function(seqT,3,500)




#### function to plot CpG Islands ####
CpG_plotter <- function(){
    
}

#### notes and comments section ####
#10-17-2017 nathan 12:53
    #Formatted the group master script
    # added sections
    #
HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIV1b = read.alignment("humanparainfluenzavirus1_F.fasta_pruned.mu.trim05", format = "fasta")
HPIV1c = read.alignment("humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05", format = "fasta")
HPIV3a = read.alignment("humanparainfluenzavirus3.fasta_pruned.mu.trim05", format = "fasta")
HPIV3b = read.alignment("humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05", format = "fasta")
#average size of each sample is 15,500

str(HPIV1a)
