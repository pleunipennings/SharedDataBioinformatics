setwd("~/Dropbox/2015_fallBioinformatics/RTutorials/Influenza")
library(seqinr)
library(ape)

NAseqs<-read.fasta("NA_alignment.fasta") #NA stands for neuraminidase, the gene for which we have the sequences 

#How to get the consensus (= most common nucleotide for each position)
NAseqsAli<-read.alignment("NA_alignment.fasta", format="fasta") #NA stands for neuraminidase, the gene for which we have the sequences 
cons<-consensus(NAseqsAli)
seqinr::consensus()

#use read.dna to get the data in matrix form, this makes it easier to count
NAseqsDNA<-read.dna("NA_alignment.fasta", format = "fasta",as.character=TRUE)

#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons<-length(which(NAseqsDNA[,1]==cons[1]))
#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans<-length(which(NAseqsDNA[,1]=="g"))

numTrans<-length(which(NAseqsDNA[,1]==transitionfunction(cons[1])))

 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
getwd()
setwd("/Users/shantothenel/Desktop/HumanBocaVirus")
library(seqinr)
library(ape)

#How to get the consensus (= most common nucleotide for each position)
frame1 <- read.fasta("frame1.txt")
boca_align <- read.alignment("frame1.txt", format="fasta")
boca_align

cons_boca<- seqinr::consensus(boca_align)
cons_boca

#use read.dna to get the data in matrix form, this makes it easier to count
boca_DNA <- read.dna("frame1.txt", format = "fasta", as.character=TRUE)
boca_DNA

#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons <- length(which(boca_DNA[,1]==cons_boca[1]))
numCons

#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans01 <- length(which(boca_DNA[,1]=="g"))
numTrans01

#mean frequency of the transition mutation

#making the transition function, based off the nucleotide transition functions
transition_vec = function(vector) {
  # %in% means if you find any of these things in the list, return 
  # what I want anyways
  for(i in 1:length(vector))
  if(vector[i] %in% c("A","a")) {
    vector[i] = "G" 
  } else if(vector[i] %in% c("G","g")) {
    vector[i] = "A" 
  } else if(vector[i] %in% c("T","t")) {
    vector[i] = "C"
  } else if(vector[i] %in% c("C","c")) {
    vector[i] = "T"
  } else if(vector[i] != T) {
    vector[i] = "ERROR"
  }
  return(vector)
}

#changing consensus into a vector
length(cons_boca)
consTotal = cons_boca[1:1917]
consTotal
transition_vec(consTotal)
#or just do this BECAUSE WOW IT WAS A VECTOR ALL ALONG LOL
transition_vec(cons_boca)

numTrans03 <-length(which(gene_DNA[,1]==transition_vec("G")))
numTrans03

numTrans02 <- length(which(gene_DNA[,1]== transition_vec(cons_boca[1])))
numTrans02


consTotalCol = transition_vec(consTotal)
consTotalCol

hiv_df = data.frame(hiv_csv$num, hiv_csv$WTnt, hiv_csv$MeanFreq, hiv_csv$TypeOfSite)
hiv_df

#probably don't need this
hiv_df$consTransition <- consTotalCol[1:984]
hiv_df

?consensus


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
getwd()
setwd("/Users/shantothenel/Desktop/HumanBocaVirus")
library(seqinr)
library(ape)

frame2 <- read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.txt")
boca_align2 <- read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.txt", format="fasta")
boca_align2

cons_boca2 <- seqinr::consensus(boca_align2, method = "threshold")
cons_boca2

boca_DNA2 <- read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.txt", format = "fasta", as.character=TRUE)
boca_DNA2

numCons2 <- length(which(boca_DNA2[,1]==cons_boca2[1]))
numCons2

numTrans2 <- length(which(boca_DNA2[,1]=="g"))
numTrans2
numTrans2 <-length(which(boca_DNA2[,1]==transition_nuc("g")))
numTrans2

numTrans022 <- length(which(boca_DNA2[,1]== transition_vec(cons_boca[1])))
numTrans022

#making the transition function, based off the nucleotide transition functions
transition_vec = function(vector) {
  # %in% means if you find any of these things in the list, return 
  # what I want anyways
  for(i in 1:length(vector))
    if(vector[i] %in% c("A","a")) {
      vector[i] = "G" 
    } else if(vector[i] %in% c("G","g")) {
      vector[i] = "A" 
    } else if(vector[i] %in% c("T","t")) {
      vector[i] = "C"
    } else if(vector[i] %in% c("C","c")) {
      vector[i] = "T"
    } else if(vector[i] != T) {
      vector[i] = "ERROR"
    }
  return(vector)
}

length(cons_boca)
consTotal = cons_boca[1:1917]
consTotal
transition_vec(consTotal)



#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons2<-length(which(boca_DNA2[,1]==cons_boca2[1]))
numCons2
#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans2<-length(which(boca_DNA2[,1]=="g"))
numTrans2

numTrans3<-length(which(boca_DNA2[,1]==transition_vec(cons_boca2[1])))
numTrans3

















