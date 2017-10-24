

read.fasta("DengueVirus1.fasta") -> DEN

length(names(DEN))

WOW <- data.frame(matrix(unlist(DEN),nrow=10689))


NROW(WOW)


WOW[nrow(WOW) + 1,] = which.max(WOW[,1])

WOW[nrow(WOW) + 1,] = DF[nrow(DF),]
