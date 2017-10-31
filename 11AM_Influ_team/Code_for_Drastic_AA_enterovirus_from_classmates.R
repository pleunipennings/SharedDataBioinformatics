library(seqinr)
library(ape)

#First we'll create the data frame
Pos<-c(1:891)
enterodata<-data.frame(Pos)
enterodata$WtNt=""
enterodata$TrNtFreq=""
enterodata$WTAA=""
enterodata$MUTAA=""
enterodata$WTAAcat=""
enterodata$MUTAAcat=""
enterodata$bigAAchange=""

#Read in shortened entero data
enteroseqs<-read.fasta("enteroshort.txt")

#Align entero data
enteroaligned<-read.alignment("enteroshort.txt", format="fasta")

#Gets DNA into matrix form
enteromatrix<-read.dna("enteroshort.txt", format="fasta", as.character = TRUE)

#Gets WTnt for each nt position
cons<-seqinr::consensus(enteroaligned)

#Insert concensus into dataframe
enterodata$WtNt=cons

#Function to return transition mutation
transition<-function(basepair){
  #basepair<-("A", "C", "T", "G"),
  if(basepair=="a") {return("g")}
  if(basepair=="g") {return ("a")}
  if(basepair=="t") {return ("c")}
  if(basepair=="c") {return ("t")}
}

#Loop to insert the freq of transition mutations
for(i in 1:nrow(enterodata)){
  enterodata$TrNtFreq[i]<-length(which(enteromatrix[,i]==transition(cons[i])))/
    (nrow(enteromatrix))
}

#Translate the consensus to AAs
translate(cons) -> wildAA

#Insert the translated consensus into data frame
count=0
for(i in 1:length(wildAA)) {
  count = count+1
  enterodata[count,]$WTAA=wildAA[i]
  count=count+1
  enterodata[count,]$WTAA=wildAA[i]
  count=count+1
  enterodata[count,]$WTAA=wildAA[i]
}

#This is the loop but it doesn't run through all the data
for(i in seq(1,891,3)){
  enterodata$MUTAA[i]=translate(c(transition(enterodata$WtNt[i]), enterodata$WtNt[i+1], enterodata$WtNt[i+2]))
  enterodata$MUTAA[i+1]=translate(c(enterodata$WtNt[i], transition(enterodata$WtNt[i+1]), enterodata$WtNt[i+2]))
  enterodata$MUTAA[i+2]=translate(c(enterodata$WtNt[i], enterodata$WtNt[i+1], transition(enterodata$WtNt[i+2])))
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

#Assign wild type AA category
for(i in 1:nrow(enterodata)){
  enterodata$WTAAcat[j]=amCat(enterodata$WTAA[j])
}

#Assign mutated AA category
for(i in 1:nrow(enterodata)){
  enterodata$MUTAAcat[j]=amCat(enterodata$MUTAA[j])
}

#Loop for drastic change or not 
for(i in 1:nrow(enterodata)){
  if (enterodata$WTAAcat[i]==enterodata$MUTAAcat[i]){
    enterodata$bigAAchange[i]= "0"
  }
  if (enterodata$WTAAcat[i]!=enterodata$MUTAAcat[i]){
    enterodata$bigAAchange[i] = "1"
  }
}
