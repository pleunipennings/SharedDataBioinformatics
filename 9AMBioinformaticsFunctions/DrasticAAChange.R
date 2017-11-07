# Drastic Amino Acid change function
# Contributors: Hassan, Avery, Emily, Angeline and Fran
# We would also like to give Dwayne Evans a shout-out for helping us with a for loop!

# Our function takes in FASTA files as input. It finds the most frequent occurance of nucleotides at each site and uses
# that data to determine the wild-type sequence for the virus. The wild-type sequence is then translated into its respective
# amino acids. Using the wild-type sequence, and amino acids we are able to determine how mutations at each site 
# will affect changes in amino acids; In our case drastic vs non-drastic changes.


# This function assumes that the user has installed the packages: seqinr, Biostrings and ape.
# The file pathway in line 21 needs to be changed to adjust for the data being analyzed

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
df=data.frame(matrix(nrow=b, ncol=10)) 
names(df)=c("wtnt", "MeanFreq", "Mutseq", "WTAA", "MUTAA", "WTcat", "Mutcat", "bigAAchange", "TypeOfSite", "CpG")




# We're defining the wild type sequence as the most common sequence
for(i in 0:(b-1)){
  x=table(seqs[,i+1])
  y=names(x[which.max(x)])
  z=((x[which.max(x)])/a)
  df[i+1,]$wtnt = y
  df[i+1,]$MeanFreq = z  #mean frequency
}


# Translate the WT DNA sequence to AA sequence

curSeq <- translate(paste(df[,]$wtnt, sep=" "), NAstring="X", ambiguous = FALSE, sens="F")
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
  df[j,]$WTcat=amCat(df[j,]$WTAA)
}

#creates Mut strain of current WTseq enters it into 
for (i in 1:b){
  curNuc <-  df[i,1]
  if(curNuc == "a" || curNuc == "g"){
    if(curNuc == "a"){
      df [i, 3] = "g"
    }else{
      df [i, 3] = "a"
    }
  }else{
    if(curNuc == "t"){
      df [i, 3] = "c"
    }else{
      df [i, 3] = "t"
    }
  }
}

# Let's make a column where the WT AA is printed three times over their codons 
df$MUTAA <-c(0)
count <- 1
i = 1
j = 1
for(i in 1:(nrow(df)/3)){ # Workspace for MutAA filling don't run this loop
  for (j in 1:3){
    if(j == 1){                #first nucleotide in codon
      df[count,]$MUTAA <- translate(firstC <- c(as.character(df[count,]$Mutseq), as.character(df[count + 1,]$wtnt), as.character(df[count +2,]$wtnt)))
    } 
    if(j == 2){                #second nucleotide in codon
      df[count,]$MUTAA <- translate(secondC <- c(as.character(df[count - 1,]$wtnt), as.character(df[count,]$Mutseq), as.character(df[count +1,]$wtnt)))
    }
    if(j == 3){                #third nucleotide in codon
      df[count, ]$MUTAA <- translate(thirdC <- c(as.character(df[count - 2,]$wtnt), as.character(df[count - 1,]$wtnt), as.character(df[count,]$Mutseq)))
    }
    count <- count + 1
    print(j)
    print(i)
    print(count)
  }
}                                           


df$MUTAA=seqinr::translate(paste(df$Mutseq, sep=" "),
                           NAstring="X", ambiguous=FALSE, sens="F")            # Translating from nucleic acid to AA

for(j in 1:b){                                              
  df[j,]$Mutcat=amCat(df[j,]$MUTAA)                                    # Categorizing the AA 
}

# Function to compare DrasticAAChanges
DrasticChange <- function(df){
  for(h in 1:nrow(df)){
    if (df[h,]$WTcat==df[h,]$Mutcat){
      df[h,]$bigAAchange= 0                      # If WT AA category = Mut AA Category, no drastic change
    }
    if (df[h,]$WTcat!=df[h,]$Mutcat){
      df[h,]$bigAAchange= 1                      # If WT AA category =/= Mut AA Category, yes drastic change
    }
  }
  return(df)
}

df <- DrasticChange(df) #runs the Drastic Change Function

