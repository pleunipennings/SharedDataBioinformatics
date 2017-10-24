setwd("~/Desktop")
library(seqinr)
read.alignment("HumanParaInfluenza/humanparainfluenzavirus1_F.fasta_pruned.mu.trim05",format="fasta")->Virus
print(Virus)
seqinr::consensus(Virus)
