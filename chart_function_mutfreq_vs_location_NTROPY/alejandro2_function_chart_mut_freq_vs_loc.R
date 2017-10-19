#Alejandro's code for chart function

#set directory and open data faile
setwd("~/Desktop")
library(seqinr)
library(ape)
read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.trim05")->SN1seq
read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.trim05",format = "fasta")->SN1seqalign
cons<-consensus(SN1seqalign)
SN1data<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.trim05", format = "fasta",as.character=TRUE)


summary(HBOVdf)
plot(data[,2],data[,3],col=data$TypeOfSite,xlab = "Number",
     ylab = "Mean Frequency")

str(data)
dim(data)
class(data)
which(data$WTAA=="A")

