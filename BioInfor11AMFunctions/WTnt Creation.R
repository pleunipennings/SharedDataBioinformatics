### Remember to set working directory and replace your own file/dataframe names!

library(seqinr)
align <- read.alignment("BKpolyomavirus_VP1.fasta.mu.trim05", format = "fasta")
WTnt <- consensus(align)
df <- cbind(WTnt)