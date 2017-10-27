setwd("~/bioinformatics/bioinformaticsproject")
library(seqinr)
library(ape)

read.csv("OverviewSelCoeff_BachelerFilter(1).csv") -> sampledata
head(sampledata)

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
str(cons)

#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format = "fasta",as.character=TRUE)
str(bocaNS1seqsDNA)
bocaNS1seqsDNA[31,]
class(bocaNS1seqsDNA])
summary(bocaNS1seqsDNA)

#Now you can use length and which and the == operator to count the number of sequences with the consensus nucleotide

BoNS1df<-data.frame("num"=c(1:5157), "WTnt"=cons)
BoNS1df

#use for loop to run through each location of 1 to 5157 nt
numCons<-length(which(bocaNS1seqsDNA[,nt]==cons[,nt]))

numCons

#Here you should use a function to determine the transition from the consensus. For now, I just put "g", but it's better to use a function! 
numTrans<-length(which(bocaNS1seqsDNA[,1]=="g"))
numTrans

#function to run to see transition from consensus
for (s in 1:5157){
    BoNS1df$newcol<-cons
}
View(BoNS1df)

meanfreq<-function(x){
    if (cons[,1]!=bocaNS1seqsDNA[,1]) {
        length(which(bocaNS1seqsDNA[,1]=="g"))
    }
    else {
        if (bocaNS1seqsDNA[,1=="a"])
    }
}

numTrans<-length(which(bocaNS1seqsDNA[,1]==transitionfunction(cons[1])))

BoNS1df<-data.frame("num"=c(1:5157), "WTnt"=cons)
BoNS1df


