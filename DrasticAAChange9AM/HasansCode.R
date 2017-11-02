#setting up enviroment

  library(seqinr)
  library(Biostrings)
  library(ape)

#reading files and setting data  
  seqs = read.dna("InfluenzaAvirus_HA_H1N1.fasta", format = "fasta", as.character=TRUE)
  a=nrow(seqs)
  b=ncol(seqs)
  

# Let's establish the new data frame
# Number of columns is arbitrary. I just want to print the WT sequence down a column
  df=data.frame(matrix(nrow=b, ncol=8)) 
  names(df)=c("WTseq", "Mutseq", "WTAA", "MutAA", "WTcat", "Mutcat", "DrasticAA", "Syn/NonSyn/Nonsense")
  
#mean frequency
  

# We're defining the wild type sequence as the most common sequence
  for(i in 0:(b-1)){
    x=table(seqs[,i+1])
    y=names(x[which.max(x)])
    df[i+1,1]=y
  }

# Translate the WT DNA sequence to AA sequence

  curSeq <- translate(paste(df[,1], sep=" "), NAstring="X", ambiguous = FALSE, sens="F")
  count <- 0
  
  for(i in 1:length(curSeq)){ # incrementing down sequence by 3 (needs work)
          count <- count + 1
          df[count,]$WTAA <- curSeq[i]
          count <- count + 1
          df[count,]$WTAA <- curSeq[i]
          count <- count + 1
          df[count,]$WTAA <- curSeq[i]
  }

#Defining AA changes
  pos <- "R|H|K"
  neg <- "D|E"
  unc <- "S|T|N|Q"
  spe <- "C|U|G|P"
  hyd <- "A|I|L|F|M|W|Y|V"
  amCat <- function(AA){
    if(regexpr(pos, AA) > 0){ return("pos") }
    if(regexpr(neg, AA) > 0){ return("neg") }
    if(regexpr(unc, AA) > 0){ return("unc") }
    if(regexpr(spe, AA) > 0){ return("spe") }
    if(regexpr(hyd, AA) > 0){ return("hyd") }
    return(5)
  }

# Categorizing the WT amino acids
for(j in 1:b){
  df[j,5]=amCat(df[j,3])
}
  
#creates Mut strain of current WTseq enters it into 
  for (i in 1:b){
   curNuc <-  df[i,1]
   
   if(curNuc == "a" || curNuc == "g"){
     
     if(curNuc == "a"){
       df [i, 2] = "g"
      }else{
       df [i, 2] = "a"
        }
      
     }else{
     
       if(curNuc == "t"){
       df [i, 2] = "c"
     }else{
       df [i, 2] = "t"
     }
   
       }
   
  }

# Let's make a column where the WT AA is printed three times over their codons 
  df$MutAA <-c(0)
  count <- 1
  i = 1
  j = 1
  for(i in 1:(nrow(df)/3)){ # Workspace for MutAA filling don't run this loop
    for (j in 1:3){
      if(j == 1){                #first nucleotide in codon
           df[count,]$MutAA <- translate(firstC <- c(as.character(df[count,]$Mutseq), as.character(df[count + 1,]$WTseq), as.character(df[count +2,]$WTseq)))
      } 
      if(j == 2){                #second nucleotide in codon
           df[count,]$MutAA <- translate(secondC <- c(as.character(df[count - 1,]$WTseq), as.character(df[count,]$Mutseq), as.character(df[count +1,]$WTseq)))
      }
      if(j == 3){                #third nucleotide in codon
           df[count, ]$MutAA <- translate(thirdC <- c(as.character(df[count - 2,]$WTseq), as.character(df[count - 1,]$WTseq), as.character(df[count,]$Mutseq)))
      }
     count <- count + 1
     print(j)
     print(i)
     print(count)
      }
  }                                           


df[,4]=seqinr::translate(paste(df[,2], sep=" "),
        NAstring="X", ambiguous=FALSE, sens="F")            # Translating from nucleic acid to AA

for(j in 1:b){                                              
  df[j,6]=amCat(df[j,4])                                    # Categorizing the AA 
}

#function for syn/non/nonsence
for (h in 1:b){
  if(df[h,7]== 0){
    df[h,8] = "syn"
  }
  if(df[h,7]==1){
    if(df[h,4]=="X"){
      df[h,8] = "non-sense"
    }
    else{
      df[h,8] = "non-syn"
    }
  }
}

# Function to compare DrasticAAChanges
DrasticChange <- function(df){
    for(h in 1:nrow(df)){
      if (df[h,]$WTcat==df[h,]$Mutcat){
        df[h,]$DrasticAA= 0                      # If WT AA category = Mut AA Category, no drastic change
      }
      if (df[h,]$WTcat!=df[h,]$Mutcat){
        df[h,]$DrasticAA= 1                      # If WT AA category =/= Mut AA Category, yes drastic change
      }
    }
  return(df)
}

#function for syn/non/nonsence
for (h in 1:b){
  if(df[h,7]== 0){
    df[h,8] = "syn"
  }
  if(df[h,7]==1){
    if(df[h,4]=="X"){
      df[h,8] = "non-sense"
    }
    else{
      df[h,8] = "non-syn"
    }
  }
}

df <- DrasticChange(df)
save(df,file="df.Rda")
load("df.Rda")
