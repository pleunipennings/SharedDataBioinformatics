setwd("~/desktop/Bioinformatics")
library(seqinr)
library(ape)

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu copy.trim05") 


#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu copy.trim05", format="fasta")


#works when using con() instead of consensus()
seqinr::consensus(bocaNS1seqsAli)->WTnt

#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu copy.trim05", format = "fasta",as.character=TRUE)

# Writes new function to create transition mutations in nucleotides
transition <- function(nt){
  if(nt=="a") {return("g")}
  if(nt=="g") {return("a")}
  if(nt=="c") {return("t")}
  if(nt=="t") {return("c")}}

# For-loop to calculate mean frequency of transition mutations for each nucleotide:

MeanFreq<-c()
for (i in 1:ncol(bocaNS1seqsDNA)){
  MeanFreq<-c(MeanFreq,(length(which(bocaNS1seqsDNA[,i]==transition(WTnt[i])))/ncol(bocaNS1seqsDNA)))
}



BoNS1df<-data.frame("num"=c(1:ncol(bocaNS1seqsDNA)),
                    WTnt,
                    MeanFreq)

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
CpG_finder(BoNS1df)

#start of CpG Graph Function 
LvsF_CpG_Printer <- function(data_frame){
  if (T) {n$makesCpG <- n$makesCpG+1}
  YCpG <- which(n$makesCpG=="2")
  NCpG <- which(n$makesCpG=="1")
  x1 <- n$MeanFreq[YCpG]
  x2 <- n$MeanFreq[NCpG]
  plot.default(x = c(x1, x2), 
               xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
               col = (as.integer(n$makesCpG)),
               log = "y"
  )#close plot.default
  
  return()#close return
}#close function

LvsF_CpG_Printer(data_frame) #run function
head(data_frame)