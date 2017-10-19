setwd("~/Dropbox/2015_fallBioinformatics/RTutorials/Influenza")
library(seqinr)
library(ape)

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05") 

#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format="fasta")
bocaNS1seqsAli
class(bocaNS1seqsAli)
tail(bocaNS1seqsAli)
bocaNS1seqsAli[3]

cons<-consensus(bocaNS1seqsAli[[1:31]])

#use read.dna to get the data in matrix form, this makes it easier to count
NAseqsDNA<-read.dna("NA_alignment.fasta", format = "fasta",as.character=TRUE)

#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons<-length(which(NAseqsDNA[,1]==cons[1]))
#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans<-length(which(NAseqsDNA[,1]=="g"))

numTrans<-length(which(NAseqsDNA[,1]==transitionfunction(cons[1])))




