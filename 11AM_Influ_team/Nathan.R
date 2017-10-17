#Nathan
#course stuff!
h1 <- read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05",format="fasta")
h1
#?read.alignment
h2<- dist.alignment(h1)
h2
nj(h1)
?nj
h3<-nj(h2
)
plot(h3)


#nate stuff
G <- read.alignment("class25Influ.txt", format = "fasta")
G

G1 <- dist.alignment(G)
G1

G2 <- nj(G1)
G2

G3 <- plot(G2)
G3

?`ape-package`
library(help = "ape")