#Bioinformatics Final Project Deshawn Hopson & Jasmine Sims - CpG v. Non-CpG Significane

library(ggplot2)
library(seqinr)
setwd("~/Final Project")

#A variable that contains the virus fasta files from our wd folder (1 of each virus)
Files <- list.files(recursive = TRUE, pattern = "trim") #Includes all files that have "trim" in name, which are all the virus files.

#Loop to run through each fasta file and output the resuts of the wilcox tests for -a and -t CpG / no-CpG transition mutations
for(i in 1:length(Files)){ #Will go the length of number of files
    
    DF <- meanFreq(Files[i]) #Will create data frame with mean freq and wtnt
    
    DF <- getWTAA(DF) #Will create WTAA column
    
    DF <- getMUTAA(DF) #Will create MUTAA colum
    
    DF <- big_aa_change(DF) #Will create big_aa_change column
    
    DF <- functionSynNonSyn(DF) #Will create syn / non-syn column
    
    DF <- CpG_finder(DF) #Will create CpG / no-CpG column
    
    ###a -> g transition mutation
    
    print(Files[i]) #Print name of file to reference outputs
    
    DF2 <- subset(DF, makesCpG=="1" & wtnt=="a") #Subset data for a->g transition that cause a CpG site
    a_cpg <- c(DF2$MeanFreq) #Apply mean frequencies of subset into a vector
    
    DF2.1 <- subset(DF, makesCpG=="0" & wtnt=="a") #Subset data for a->g transition that don't cause a CpG site
    a_nocpg <- c(DF2.1$MeanFreq) #Apply mean frequencies of subset into a vector
    
    print(wilcox.test(a_cpg, a_nocpg)) #Print wilcox.test results for these two mean frequency vectors
    
    ###t -> c transition mutation
    
    DF3 <- subset(DF, makesCpG=="1" & wtnt=="t") #Subset data for t->c transition that cause a CpG site
    t_cpg <- c(DF3$MeanFreq) #Apply mean frequencies of subset into a vector
    
    DF3.1 <- subset(DF, makesCpG=="0" & wtnt=="t") #Subset data for t->c transition that don't cause a CpG site
    t_nocpg <- c(DF3.1$MeanFreq) #Apply mean frequencies of subset into a vector
    
    print(wilcox.test(t_cpg, t_nocpg)) #Print wilcox.test results for these two mean frequenciy vectors
    
}



#///// Virus of interest (Dengue Virus)

Dengue <- meanFreq("DengueVirus1.fasta_pruned.mu.trim05") #Create data frame with mean freq and wtnt 

Dengue <- getWTAA(Dengue) #Create WTAA column

Dengue <- getMUTAA(Dengue) #Create MUTAA column

Dengue <- big_aa_change(Dengue) #Create big_aa_change column

Dengue <- functionSynNonSyn(Dengue) #Create syn / non-syn column

Dengue <- CpG_finder(Dengue) #Create CpG / no-CpG column

### a -> g transition mutation

Dengue2 <- subset(Dengue, makesCpG=="1" & wtnt=="a") #Subset data for a->g transition that cause a CpG site
a_cpg_dengue <- c(Dengue2$MeanFreq) #Apply mean frequencies of subset into a vector

Dengue2.1 <- subset(Dengue, makesCpG=="0" & wtnt=="a") #Subset data for a->g transition that don't cause a CpG site
a_nocpg_dengue <- c(Dengue2.1$MeanFreq) #Apply mean frequencies of subset into a vector

wilcox.test(a_cpg, a_nocpg) #Determines significance of the difference between these variables

### t -> c transition mutation

Dengue3 <- subset(Dengue, makesCpG=="1" & wtnt=="t") #Subset data for t->c transition that cause a CpG site
t_cpg_dengue <- c(Dengue3$MeanFreq) #Apply mean frequencies of subset into a vector

Dengue3.1 <- subset(Dengue, makesCpG=="0" & wtnt=="t") #Subset data for t->c transition that don't cause a CpG site
t_nocpg_dengue <- c(Dengue3.1$MeanFreq) #Apply mean frequencies of subset into a vector

wilcox.test(t_cpg, t_nocpg) #Determines significance of the difference between these variables



#///// PLotting Dengue

#Histogram of number of sequences with specified mean frequency of a->g transition mutation that DO NOT cause a CpG site
hist(a_nocpg, breaks = 40, col = rgb(1,.7,0,0.6), xlim = c(.1,.5), ylim = c(0, 35),
     main = "A to G Transition Mutations of Dengue Virus", 
     xlab = "Mean Frequency", ylab = "Sequence Count", cex.lab = 1.5, cex.main = 2)
#Histogram of number of sequences with specified mean frequency of a->g transition mutation that DO cause a CpG site (layered over previous)
hist(a_cpg, breaks = 40, col = rgb(0,1,0,0.6), add = TRUE)
#Legend for mutations that do or do not create a CpG site
legend("topright", c("No", "Yes"), col=c(rgb(1,.7,0,0.6), rgb(0,1,0,0.6)), lwd = 3,
       title = "Makes CpG", cex = 1, pch = 2)

#Histogram of number of sequences with specified mean frequency of t->c transition mutation that DO NOT cause a CpG site
hist(t_nocpg, breaks = 40, col = rgb(1,.7,0,0.6), xlim = c(.1,.5), ylim = c(0, 30),
     main = "T to C Transition Mutations of Dengue Virus",
     xlab = "Mean Frequency", ylab = "Sequence Count", cex.lab = 1.5, cex.main = 2)
#Histogram of number of sequences with specified mean frequency of t->c transition mutation that DO cause a CpG site (layered over previous)
hist(t_cpg, breaks = 40, col = rgb(0,1,0,0.6), add = TRUE)
#Legend for mutations that do or do not create a CpG site
legend("topright", c("No", "Yes"), col=c(rgb(1,.7,0,0.6), rgb(0,1,0,0.6)), lwd = 3,
       title = "Makes CpG", cex = 1, pch = 2)

