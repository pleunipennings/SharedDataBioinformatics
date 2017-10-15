library(seqinr)
setwd("~/SharedDataBioinformatics/InfluenzaVirus/InfluenzaVirus")

read.fasta("InfluenzaAvirus_HA_H1N1.fasta.mu.fasta") -> seqData

seqData[1]

seqData[[1]][3]

seqData[[12500]]

influenzaTable <- matrix(ncol = getLength(seqData[1]), nrow = length(seqData))

for(i  in 1:length(seqData)){
    for (j in 1:getLength(seqData[i])){
        
        influenzaTable[i,j] <- seqData [[i]][j]
        
    }
}

