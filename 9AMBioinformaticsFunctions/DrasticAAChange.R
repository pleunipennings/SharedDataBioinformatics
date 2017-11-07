# Drastic Amino Acid change function
# Contributors: Hassan, Avery, Emily, Angeline, Edgar and Fran
# We would also like to give Dwayne Evans a shout-out for helping us with a for loop!

# Our function takes in FASTA files as input. It finds the most frequent occurance of nucleotides at each site and uses
# that data to determine the wild-type sequence for the virus. The wild-type sequence is then translated into its respective
# amino acids. Using the wild-type sequence, and amino acids we are able to determine how mutations at each site 
# will affect changes in amino acids; In our case drastic vs non-drastic changes.


# This function assumes that the user has installed the packages: seqinr, Biostrings and ape.
# The file pathway in line 21 needs to be changed to adjust for the data being analyzed


DrasticChange <- function(df){

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

  df$Mutcat <- 0  #setting up new column necessary to run function
  df$WTcat <- 0
  df$bigAAchange <- 0
 
  for(j in 1:b){
    df[j,]$WTcat=amCat(df[j,]$WTAA)
  }
  
for(j in 1:nrow(df)){                                              
  df[j,]$Mutcat=amCat(df[j,]$MUTAA)                                    # Categorizing the AA 
}

# Function to compare DrasticAAChanges

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
df<-DrasticChange(df)
