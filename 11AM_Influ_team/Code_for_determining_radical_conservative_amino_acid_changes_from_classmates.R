setwd("~/desktop") 

library(ape)

library(pylr)

#First we'll create the data frame
Pos<-c(1:891)
enterodata<-data.frame(Pos)
enterodata$WtNt=""
enterodata$TrNtFreq=""
enterodata$WTAA=""
enterodata$MUTAA=""
enterodata$bigAAchange=""

#Read in shortened entero data
enteroseqs<-read.fasta("enteroshort.txt")

#Align entero data
enteroaligned<-read.alignment("enteroshort.txt", format="fasta")


#Gets WTnt for each nt position
cons<-seqinr::consensus(enteroaligned)

#Insert concensus into dataframe
enterodata$WtNt=cons
WTAA<-list()
MUTAA<-list()

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
  WtAA[x] <- new_AA
  WtAA[x+1] <- new_AA
  WtAA[x+2] <- new_AA
}

#places respective things in columns
enterodata$WTAA<- WTAA
enterodata$MUTAA <- MUTAA

#Loop to insert the freq of transition mutations
for(i in 1:nrow(enterodata)){
  enterodata$TrNtFreq[i]<-length(which(enteromatrix[,i]==transition(cons[i])))/
    (nrow(enteromatrix))
}

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

#for loop to categorize WTAA
for(i in 1:nrow(enterodata)){                                              
  enterodata[i,7]=amCat(enterodata[i,4])                                    
}

#for loop to categorize MUTAA
for(i in 1:nrow(enterodata)){                                              
  enterodata[i,8]=amCat(enterodata[i,5])                                 
}


#for loop to determine drastic change or not
for(i in 1:nrow(enterodata)){
  if (enterodata$WTAAcat[i]==enterodata$MUTAAcat[i]){
    enterodata$bigAAchange[i]= "0"
  }
  if (enterodata$WTAAcat[i]!=enterodata$MUTAAcat[i]){
    enterodata$bigAAchange[i] = "1"
  }
}

