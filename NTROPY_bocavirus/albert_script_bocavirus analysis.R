setwd("~/bioinformatics/bioinformaticsproject")
library(seqinr)
library(ape)

#bi2sfsu@gmail.com
#bi2sfsu!

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05") 
#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format="fasta")
#works when using con() instead of consensus()
seqinr::consensus(bocaNS1seqsAli)->cons

#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format = "fasta",as.character=TRUE)
# Writes new function to create transition mutations in nucleotides
transition <- function(nt){
    if(nt=="a") {return("g")}
    if(nt=="g") {return("a")}
    if(nt=="c") {return("t")}
    if(nt=="t") {return("c")}}
# For-loop to calculate mean frequency of transition mutations for each nucleotide:
MeanFreq<-c()
for (i in 1:ncol(bocaNS1seqsDNA)){
    MeanFreq<-c(MeanFreq,(length(which(bocaNS1seqsDNA[,i]==transition(WTnt[i])))/ncol(bocaNS1seqsDNA)))}

#basic data frame with mean freq
BoNS1df<-data.frame("num"=c(1:ncol(bocaNS1seqsDNA)),"WTnt" = cons,MeanFreq)

#CPG sites input
CpG_finder <- function(new_virus_data){
    #reads data into function as CSV file
    virus_data <- new_virus_data
    
    #singles out column of nucleotides -- this assumes that the new datafile has the same column headers such as WTnt as the HIV file does
    WTnt <- virus_data$WTnt
    
    #gets the length of the data file for use in a loop
    data_length <- nrow(virus_data)
    
    #creates an empty vector (like a list) of the same length as the data file to be used to record the results of the loop below
    makesCpG <- vector(mode = "numeric", length = data_length)
    
    #loop that determines if a CpG could occur due to mutation at each spot in the list of nucleotides (WTnt)
    #loops from row 1 to the the last row of the column WTnt
    for(x in 1:data_length){
        
        #assigns a name (current_nuc) to the nucleotide at row x in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
        current_nuc <- toupper(WTnt[x])
        
        #assigns a name (current_neighbor) to the nucleotide in the next row down in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
        current_neighbor <- toupper(WTnt[x + 1])
        
        #begins a loop that compares the two nucleotides that were just isolated
        #this if statement is only evaluated if the current nucleotide is a T
        if(current_nuc == "T"){
            #this if statement is only evaluated if the neighbor is a G
            if(current_neighbor == "G"){
                #if the current nucleotide is a T and the neighbor is a G (thus a mutation can cause a CpG), then the of the current nucleotide in the list of 0s we made is changed to a 1
                makesCpG[x] = 1
            }
            #otherwise, nothing is changed
            else{}
        }
        #repeat the same loop as before, but modified slightly to look for CA pairs
        if(current_nuc == "C"){
            if(current_neighbor == "A"){
                #here, instead of changing the position of the current nuc to a 1, we change the neighbor's position because that is where the mutation would be, hence x+1
                makesCpG[x+1] = 1
            }
            else{}
        }
        else{}
        #print(c("Data length is ", nrow(virus_data), "List length is ", length(makesCpG)))
    }
    
    #append the list of 0s and 1s we created to the end of the data file we imported
    virus_data$makesCpG <- makesCpG
    
    #when the entire function is done, return the amended data file
    return(virus_data)
}
CpG_finder(BoNS1df)->BoNS1df

#make variables for for-loop for getting WTaa, wild type amino acid,  and MUTaa, mutated amino acid 
WtAA<-list()
MUTAA<-list()
#for loop for transition mutation at each position and translating it 
for(x in seq(1, length(cons), 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    mutated_codon <- codon
    if(codon[1] == "a"){
        mutated_codon <- replace(x=mutated_codon, values=c("g", codon[2], codon[3]))
    }
    if(codon[1] == "g"){
        mutated_codon <- replace(x=mutated_codon, values=c("a", codon[2], codon[3]))
    }
    if(codon[1] == "c"){
        mutated_codon <- replace(x=mutated_codon, values=c("t", codon[2], codon[3]))
    }
    if(codon[1] == "t"){
        mutated_codon <- replace(x=mutated_codon, values=c("c", codon[2], codon[3]))
    }
    MUTAA[x] <- translate(mutated_codon)
    
    if(codon[2] == "a"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "g", codon[3]))
    }
    if(codon[2] == "g"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "a", codon[3]))
    }
    if(codon[2] == "c"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "t", codon[3]))
    }
    if(codon[2] == "t"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "c", codon[3]))
    }
    MUTAA[x+1] <- translate(mutated_codon)
    
    if(codon[3] == "a"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "g"))
    }
    if(codon[3] == "g"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "a"))
    }
    if(codon[3] == "c"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "t"))
    }
    if(codon[3] == "t"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "c"))
    }
    MUTAA[x+2] <- translate(mutated_codon)
}

for(x in seq(1, length(cons) - 2, 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    new_AA <- translate(codon)
    WtAA[x] <- new_AA
    WtAA[x+1] <- new_AA
    WtAA[x+2] <- new_AA
}
#Insert WTaa and MTaa into dataframe
View(BoNS1df)
WtAA->BoNS1df$WtAA
MUTAA->BoNS1df$MUTAA
