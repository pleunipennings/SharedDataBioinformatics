setwd("~/bioinformatics/bioinformaticsproject")
library(seqinr)
library(ape)

#bi2sfsu@gmail.com
#bi2sfsu!

read.csv("OverviewSelCoeff_BachelerFilter(1).csv") -> sampledata

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05") 
bocaNS1seqs

#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format="fasta")
class(bocaNS1seqsAli)

#works when using con() instead of consensus()
seqinr::consensus(bocaNS1seqsAli)->WTnt
summary(WTnt)
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
    MeanFreq<-c(MeanFreq,(length(which(bocaNS1seqsDNA[,i]==transition(WTnt[i])))/ncol(bocaNS1seqsDNA)))
}

summary(MeanFreq)
str(MeanFreq)
length(MeanFreq)

BoNS1df<-data.frame("num"=c(1:ncol(bocaNS1seqsDNA)),
                    WTnt,
                    MeanFreq)
BoNS1df
View(BoNS1df)

