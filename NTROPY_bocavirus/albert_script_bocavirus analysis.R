setwd("~/bioinformatics/HumanBocaVirus/HumanBocaVirus")
library(seqinr)
library(ape)

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05") 
bocaNS1seqs
class(bocaNS1seqs)
class(read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05") )
tail(bocaNS1seqs)
#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format="fasta")
class(bocaNS1seqsAli)

#works when using con() instead of consensus()
cons<-con(bocaNS1seqsAli)
cons

#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format = "fasta",as.character=TRUE)
bocaNS1seqsDNA
bocaNS1seqsDNA[[17]]
class(bocaNS1seqsDNA)
summary(bocaNS1seqsDNA)

#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons<-length(which(bocaNS1seqsDNA[,1]==cons[1]))
numCons
#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans<-length(which(NAseqsDNA[,1]=="g"))

numTrans<-length(which(NAseqsDNA[,1]==transitionfunction(cons[1])))




