#attempt 2
library(Biostrings)
library(ape)
library(seqinr)
setwd("~/Desktop/Midterm/Enteroviruses2")
bk <- read.dna("EnterovirusC_VP2.fasta_pruned.mu.fasta", format = "fasta", as.character=TRUE)
a=nrow(bk)
b=ncol(bk)

#creates Data Frame
df=data.frame(matrix(nrow=b, ncol=9)) 
names(df)=c( "WTseq","Freq", "Mutseq", "WTAA", "MutAA", "WTcat", "Mutcat", "DrasticAA", "Syn/NonSyn/Nonsense" )

#Finds the WT by finding the most common nucleotide
for(i in 0:(b-1)){
  x=table(bk[,i+1])
  y=names(x[which.max(x)])
  df[i+1,1]=y
}


#Finds the WT Amino Acid Seq
WTAminoAcids=seqinr::translate(paste(df[,1], sep=" "), frame = 0, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE)
    
df[,4] = rep(WTAminoAcids, each= 3, times = 1)

#Defining AA changes
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

# Categorizing the WT amino acids
for(j in 1:b){
  df[j,6]=amCat(df[j,4])
}

# Let's pull up a different sequence for comparison

df[,3]=bk[2,]                                            # Picking the sequence, 89 is arbitrary
MTAminoAcid=seqinr::translate(paste(df[,3], sep=" "),
                         NAstring="X", ambiguous=FALSE, sens="F")# Translating from nucleic acid to AA

df[,5] = rep(MTAminoAcid, each= 3, times = 1)

# Categorizing the AA
for(j in 1:b){                                              
  df[j,7]=amCat(df[j,5])                                     
}

# Function to compare
for(h in 1:b){
  if (df[h,4]==df[h,5]){
    df[h,8]=0                      # If WT AA category = Mut AA Category, no drastic change
  }
  if (df[h,4]!=df[h,5]){
    df[h,8]=1                      # If WT AA category =/= Mut AA Category, yes drastic change
  }
}


#function for syn/non/nonsence
for (h in 1:b){
  if(df[h,8]== 0){
    df[h,9] = "syn"
  }
  if(df[h,8]==1){
    if(df[h,5]=="X"){
      df[h,9] = "non-sense"
    }
    else{
      df[h,9] = "non-syn"
    }
  }
}

# creates dataframe containing all data
write.csv(df,"df_data.csv")
save(df,file="df_data.Rda")
load("df_data.Rda")
