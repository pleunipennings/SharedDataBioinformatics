
#read Boca virus file 
read.csv("Bocavirus1.txt")
#save Boca file as a variable
library(seqinr)
library(ape)
read.alignment("Bocavirus1.txt", format = "fasta")-> DNABOCANS1
#read Boca file in fasta format and align file and save as a variable
DNABOCANS1 <-read.dna("Bocavirus1.txt", format = "fasta",as.character=TRUE)
#read Boca file in dna matrix format to make it easier to read/ see nucleotides

str(DNABOCANS1)

#First we'll create the data frame
Pos<-c(1:5157)
bocadata<-data.frame(Pos)
bocadata$WtNt=""
bocadata$TrNtFreq=""
bocadata$WTAA=""
bocadata$MUTAA=""
bocadata$WTAAcat=""
bocadata$MUTAAcat=""
bocadata$bigAAchange=""

#read Boca file in dna matrix format to make it easier to read/ see nucleotides
DNABOCANS1 <-read.dna("Bocavirus1.txt", format = "fasta",as.character=TRUE)

#Gets WTnt for each nt position
seqinr::consensus(DNABOCANS1)->cons

#Insert concensus into dataframe
bocadata$WtNt<-cons
WTAA<-NULL
MUTAA<-NULL

#Function to return transition mutation
transition<-function(basepair){
  #basepair<-("A", "C", "T", "G"),
  if(basepair=="a") {return("g")}
  if(basepair=="g") {return ("a")}
  if(basepair=="t") {return ("c")}
  if(basepair=="c") {return ("t")}
}

#Loop to insert the freq of transition mutations
for(i in 1:nrow(bocadata)){
  bocadata$TrNtFreq[i]<-length(which(DNABOCANS1[,i]==transition(cons[i])))/
    (nrow(DNABOCANS1))
}

#for loop for transitioning each position and translating it 
for(x in seq(1, length(cons), 3)){
  codon <- c(cons[x], cons[x+1], cons[x+2])
  mutated_codon <- codon
  if(codon[1] == "a"){
    mutated_codon <- replace(x=mutated_codon, values=c("g", codon[2], codon[3]))
  }
  if(codon[1] == "g"){
    mutated_codon <- replace(x=mutated_codon, values=c("a", codon[2], codon[3]))
  }
  if(codon[1] == "c"){
    mutated_codon <- replace(x=mutated_codon, values=c("t", codon[2], codon[3]))
  }
  if(codon[1] == "t"){
    mutated_codon <- replace(x=mutated_codon, values=c("c", codon[2], codon[3]))
  }
  MUTAA[x] <- translate(mutated_codon)
  
  if(codon[2] == "a"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "g", codon[3]))
  }
  if(codon[2] == "g"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "a", codon[3]))
  }
  if(codon[2] == "c"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "t", codon[3]))
  }
  if(codon[2] == "t"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "c", codon[3]))
  }
  MUTAA[x+1] <- translate(mutated_codon)
  
  if(codon[3] == "a"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "g"))
  }
  if(codon[3] == "g"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "a"))
  }
  if(codon[3] == "c"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "t"))
  }
  if(codon[3] == "t"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "c"))
  }
  MUTAA[x+2] <- translate(mutated_codon)
}

for(x in seq(1, length(cons) - 2, 3)){
  codon <- c(cons[x], cons[x+1], cons[x+2])
  new_AA <- translate(codon)
  WTAA[x] <- new_AA
  WTAA[x+1] <- new_AA
  WTAA[x+2] <- new_AA
}

#places respective things in columns
DNABOCANS1$WTAA<- WTAA
DNABOCANS1$MUTAA <- MUTAA

#Amino Acid Changes 
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
  if(regexpr(pos, AA) > 0){ return(0) }
  if(regexpr(neg, AA) > 0){ return(1) }
  if(regexpr(unc, AA) > 0){ return(2) }
  if(regexpr(spe, AA) > 0){ return(3) }
  if(regexpr(hyd, AA) > 0){ return(4) }
  return(5)
}

#Assign wild type AA category
for(j in 1:nrow(bocadata)){
  bocadata$WTAAcat[j]=amCat(bocadata$WTAA[j])
}

#Assign mutated AA category
for(j in 1:nrow(bocadata)){
  bocadata$MUTAAcat[j]=amCat(bocadata$MUTAA[j])
}

#Loop for drastic change or not 
for(i in 1:nrow(bocadata)){
  if (bocadata$WTAAcat[i]==bocadata$MUTAAcat[i]){
    bocadata$bigAAchange[i]= "0"
  }
  if (bocadata$WTAAcat[i]!=bocadata$MUTAAcat[i]){
    bocadata$bigAAchange[i] = "1"
  }
}

