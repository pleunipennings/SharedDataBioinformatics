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


#store lenght of the sequence into a variable called X
X<-length(aln)
X

#counts the number of As,Cs,Gs and Ts in the sequence 
a<-count(aln,1)

#calculates frequency of A's 
a_freq<-(a[1]/X)

c_freq<-(a[2]/X)

g_freq<-(a[3]/X)

t_freq<-(a[4]/X)

t_freq
g_freq
c_freq
a_freq





