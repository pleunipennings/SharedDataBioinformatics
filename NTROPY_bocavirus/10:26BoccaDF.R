setwd("~/desktop")
library(seqinr)
library(ape)
SN1seq<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05")
SN1seqAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format="fasta")
SN1cons<-seqinr::consensus(SN1seqAli, method = "majority")
NS1seqDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format = "fasta",as.character=TRUE)
numCons<-length(which(SN1seqDNA[,1]==SN1cons[1]))
Transition<-function(nuc){
  if(nuc %in% c("a", "A")) {mutnuc="G"}
  if(nuc %in% c("g","G")) {mutnuc="A"}
  if(nuc %in% c("c","C")) {mutnuc="T"}
  if(nuc %in% c("t","T")) {mutnuc="C"} 
  
  return(mutnuc) 
}
numTrans<-length(which(SN1seqDNA[,1]==Transition(SN1cons[1])))
Transition("a")

#make data frame
class(SN1seqDNA)
dim(SN1seqDNA)
a<-c(1:5157)
bocans1df<-data.frame(a)
bocans1df
