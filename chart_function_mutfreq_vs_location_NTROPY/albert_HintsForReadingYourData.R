setwd("~/bioinformatics/bioinformaticsproject")
library(seqinr)
library(ape)

bvseqs<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05") 
#bv stands for neuraminidase, the gene for which we have the sequences 
bvseqs

#How to get the consensus (= most common nucleotide for each position)
bvseqsAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format="fasta") 
bvseqsAli
#consensus fxn to find the most common nucleotide
cons<-consensus(bvseqsAli)
cons
class(cons)

#use read.dna to get the data in matrix form, this makes it easier to count
bvseqsDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format = "fasta",as.character=TRUE)
bvseqsDNA
bvseqsDNA[,1]
#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide
numCons<-length(which(bvseqsDNA[,1]==cons[1]))
numCons
#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans<-length(which(bvseqsDNA[,1]=="a"))
numTrans

transition<-function(nuc){
    if(nuc %in% c("a", "A")) {mutnuc="G"}
    if(nuc %in% c("g","G")) {mutnuc="A"}
    if(nuc %in% c("c","C")) {mutnuc="T"}
    if(nuc %in% c("t","T")) {mutnuc="C"} 
    return(mutnu)
}

numTrans<-length(which(NAseqsDNA[,1]==transitionfunction(cons[1])))




