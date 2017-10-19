#setting up enviroment

  library(seqinr)
  library(Biostrings)
  library(ape)

#reading files and setting data  
  seqs = read.dna("InfluenzaAvirus_HA_H1N1.fasta.mu.fasta", format = "fasta", as.character=TRUE)
  a=nrow(seqs)
  b=ncol(seqs)

# Let's establish the new data frame
# Number of columns is arbitrary. I just want to print the WT sequence down a column
  df=data.frame(matrix(nrow=b, ncol=7)) 
  names(df)=c("WTseq", "Mutseq", "WTAA", "MutAA", "WTcat", "Mutcat", "DrasticAA")

# We're defining the wild type sequence as the most common sequence
  for(i in 0:(b-1)){
    x=table(seqs[,i+1])
    y=names(x[which.max(x)])
    df[i+1,1]=y
  }

# Translate the WT RNA sequence to AA sequence

  curSeq <- translate(paste(df[,1], sep=" "),NAstring="X", ambiguous=FALSE, sens="F")
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

# Let's pull up a different sequence for comparison
df[,2]=seqs[89,]                                            # Picking the sequence, 89 is arbitrary
df[,4]=seqinr::translate(paste(df[,2], sep=" "),
        NAstring="X", ambiguous=FALSE, sens="F")            # Translating from nucleic acid to AA
for(j in 1:b){                                              
  df[j,6]=amCat(df[j,4])                                    # Categorizing the AA 
}

# Function to compare
for(h in 1:b){
  if (df[h,5]==df[h,6]){
    df[h,7]= "0"                      # If WT AA category = Mut AA Category, no drastic change
  }
  if (df[h,5]!=df[h,6]){
    df[h,7]= "X"                      # If WT AA category =/= Mut AA Category, yes drastic change
  }
}
