setwd("~/Dropbox/2015_fallBioinformatics/RTutorials/Influenza")
library(seqinr)
library(ape)

NAseqs<-read.fasta("NA_alignment.fasta") #NA stands for neuraminidase, the gene for which we have the sequences 

#How to get the consensus (= most common nucleotide for each position)
NAseqsAli<-read.alignment("NA_alignment.fasta", format="fasta") #NA stands for neuraminidase, the gene for which we have the sequences 
cons<-consensus(NAseqsAli)

#use read.dna to get the data in matrix form, this makes it easier to count
NAseqsDNA<-read.dna("NA_alignment.fasta", format = "fasta",as.character=TRUE)

#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons<-length(which(NAseqsDNA[,1]==cons[1]))
#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans<-length(which(NAseqsDNA[,1]=="g"))

numTrans<-length(which(NAseqsDNA[,1]==transitionfunction(cons[1])))




