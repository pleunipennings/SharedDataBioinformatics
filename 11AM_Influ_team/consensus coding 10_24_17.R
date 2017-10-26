library(seqinr)
library(ape)
#-- Read HIV protease alignment
###aln <- read.fasta("BKpolyomavirus_VP1.fasta.mu.trim05.txt")

aln2 <- aln[[1]]

aln<-read.alignment("BKpolyomavirus_VP1.fasta.mu.trim05.txt", format = "fasta", forceToLower = TRUE)

#aln<-read.dna("BKpolyomavirus_VP1.fasta.mu.trim05.txt")

aln<-as.matrix.alignment(aln)
# Generate consensus
con<-seqinr::consensus(aln)
con

#store lenght of the sequence into a variable called Total_nuceleotide
Total_nucleotides<-length(aln)
Total_nucleotides

#counts the number of As,Cs,Gs and Ts in the sequence and store it in nucleotide 
count(aln,1)
nucleotide<-count(aln,1)

#calculates frequency of As, Cs, Gs, Ts. 
a_freq<-(nucleotide[1]/Total_nucleaotides)

c_freq<-(nucleotide[2]/Total_nucleaotides)

g_freq<-(nucleotide[3]/Total_nucleaotides)

t_freq<-(nucleotide[4]/Total_nucleaotides)

t_freq
g_freq
c_freq
a_freq





