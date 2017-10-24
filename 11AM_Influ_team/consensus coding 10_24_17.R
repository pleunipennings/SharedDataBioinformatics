library(seqinr)
library(ape)
#-- Read HIV protease alignment
###aln <- read.fasta("BKpolyomavirus_VP1.fasta.mu.trim05.txt")

aln2 <- aln[[1]]

aln<-read.alignment("BKpolyomavirus_VP1.fasta.mu.trim05.txt", format = "fasta", forceToLower = TRUE)

aln<-read.dna("BKpolyomavirus_VP1.fasta.mu.trim05.txt")

aln<-as.matrix.alignment(aln)
# Generate consensus
con<-seqinr::consensus(aln)
print(con$seq)

##### from nate for data.frame() #####
n <- data.frame(Pos = c(1:length(seqinr::consensus(HPIV1a))),
                WTnt = (seqinr::consensus(HPIV1a)),
                Trans = c("frequency()")
)#closes data.frame
