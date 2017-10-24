setwd("~/bioinformatics/HumanBocaVirus/HumanBocaVirus")
library(seqinr)
library(ape)

#boca NS1 stands for bocavirus NS1 gene
bocaVP1seqs<-read.fasta("HumanBocavirus1_VP1.fasta_pruned.mu.trim05") 

#How to get the consensus (= most common nucleotide for each position)
bocaVP1seqsAli<-read.alignment("HumanBocavirus1_VP1.fasta_pruned.mu.trim05",format="fasta")
bocaNS1seqsAli

#CANT MAKE IT WORK returns error "Error: $ operator is invalid for atomic vectors"
seqinr::consensus(bocaVP1Ali)
cons<-consensus(bocaVP1seqs)

#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format = "fasta",as.character=TRUE)
bocaNS1seqsDNA
class(bocaNS1seqsDNA)
summary(bocaNS1seqsDNA)

#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons<-length(which(bocaNS1seqsDNA[,1]==cons[1]))
#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans<-length(which(NAseqsDNA[,1]=="g"))

numTrans<-length(which(NAseqsDNA[,1]==transitionfunction(cons[1])))