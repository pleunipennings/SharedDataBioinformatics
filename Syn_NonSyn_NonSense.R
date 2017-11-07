
##Our team utilized Deshawn's code for syn/non-syn/nonsense. THIS CODE DOES NOT COMPLETELY WORK.
##Fernando,Krystal,Jasmine,Anjani,Darryl 

#Obtains the library files
library(seqinr)
library(ape)

#Read in shortened entero data
enteroseqs<-read.fasta("InfluenzaBvirus_NA.fasta.mu.trim05")

#Create the data frame
Pos<-c(1:1401)
enterodata<-data.frame(Pos)
enterodata$WtNt=""
enterodata$TrNtFreq=""
enterodata$WTAA=""
enterodata$MUTAA=""
enterodata$WTAAcat=""
enterodata$MUTAAcat=""
enterodata$bigAAchange=""

#Aligns entero data
enteroaligned<-read.alignment("InfluenzaBvirus_NA.fasta.mu.trim05", format="fasta")

#Gets DNA into matrix form
enteromatrix<-read.dna("InfluenzaBvirus_NA.fasta.mu.trim05", format="fasta", as.character = TRUE)

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

####Determining the category of WTAA
WT_cat<-c()
for(i in 1:length(enterodata$WTAA)){
  test<-amCat(enterodata$WTAA[i])
  WT_cat<-c(WT_cat,test)
}

WT_cat

####Storing the result of the amCat function into the enterodata data frame 
for(i in 1:length(WT_cat)){
  if(WT_cat[i] == 0){enterodata$WTAAcat[i]<- "pos"}
  if(WT_cat[i] == 1){enterodata$WTAAcat[i]<- "neg"}
  if(WT_cat[i] == 2){enterodata$WTAAcat[i]<- "unc"}
  if(WT_cat[i] == 3){enterodata$WTAAcat[i]<- "spe"}
  if(WT_cat[i] == 4){enterodata$WTAAcat[i]<- "hyd"}
  if(WT_cat[i] == 5){enterodata$WTAAcat[i]<- "stop codon"}
}

class(enterodata$WTAAcat)
#Assign wild type AA category
#for(i in 1:nrow(enterodata)){
# enterodata$WTAAcat[j]=amCat(enterodata$WTAA[j])
#}

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



#  Note that this depends on the columns MUTAA and WTAA being already complete!

functionSynNonSyn<-function(DF){
  if (length(which(names(DF)=="MUTAA"))==0){
    print("Oh oh there is a problem. No MUTAA column!")
    return(0)}
  for (h in 1:nrow(DF)){
    if(DF$MUTAA[h]== DF$WTAA[h]){
      #   if(DF[h,"MUTAA"]== DF[h,"WTAA"]){
      DF[h,"TypeOfSite"] = "syn"
    }
    if(DF[h,"MUTAA"] != DF[h,"WTAA"]){
      if(DF[h,"MUTAA"]=="*"){
        DF[h,"TypeOfSite"] = "nonsense"}
      else {
        DF[h,"TypeOfSite"] = "nonsyn"
      }
    }
  }
  DF$TypeOfSite<-as.factor(DF$TypeOfSite)
}
functionSynNonSyn(enterodata)


source("BioInfor11AMFunctions/functionSynNonSyn.R")
