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
df=data.frame(matrix(nrow=b, ncol=10)) 
names(df)=c("WTseq", "Freq", "Mutseq", "WTAA", "MutAA", "WTcat", "Mutcat", "DrasticAA", "TypeofSite", "CpG")

#mean frequency


# We're defining the wild type sequence as the most common sequence
for(i in 0:(b-1)){
  x=table(seqs[,i+1])
  y=names(x[which.max(x)])
  z=((x[which.max(x)])/a)
  df[i+1,]$WTseq = y
  df[i+1,]$Freq = z
}

# Translate the WT DNA sequence to AA sequence

curSeq <- translate(paste(df[,]$WTseq, sep=" "), NAstring="X", ambiguous = FALSE, sens="F")
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


df$MutAA=seqinr::translate(paste(df$Mutseq, sep=" "),
                           NAstring="X", ambiguous=FALSE, sens="F")            # Translating from nucleic acid to AA

for(j in 1:b){                                              
  df[j,]$Mutcat=amCat(df[j,]$MutAA)                                    # Categorizing the AA 
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

df <- DrasticChange(df) #runs the Drastic Change Function

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

#Looks for pattern tg in the data STRING, sets location sites to variable TG
TG = gregexpr(pattern ='tg', paste(df$WTseq, collapse = ''))
BELL = data.frame(matrix(unlist(TG)))

# create new column name CPG with values zero
df$CPG <- 0

# Inserting value 1 in every TG site using BELL
df[BELL[,1],"CPG"] <- 1

# For CA sites: 

#Looks for pattern ca in the data STRING, sets location sites to variable CA
CA = gregexpr(pattern ='ca', paste(df$WTseq, collapse = ''))
CASITES = data.frame(matrix(unlist(CA)))

# Inserting value 1 in every CA site using CA
df[CASITES[,1]+1,"CPG"] <- 1
View(df)

plot(
  #x vector
  df[,0],
  #y vector
  df[,2],
  #make black empty circles as symbol
  pch=21,
  #make outline of symbol black
  col= "black",
  #fill inside of point with color by factor category "TypeOfSite" bg=
  bg=df[,8],
  #Title label
  main = "HIV Practice Data",
  #x axis label
  xlab ="Location on Sequence", 
  #y axis label
  ylab ="Mean Frequency of Mutation",
  #cex changes point size
  cex=2,
  #grid superimposes grid onto plot 
  #nx and ny describes x and y axis 
  #NA will automatically set x-axis to default plot ticks 
  grid(nx = NA, ny = NULL, col = "black", lty = "dotted",
       lwd = par("lwd"), equilogs = TRUE),
  #supress y axis drawing by plot fxn, put # in front to not supress
  yaxt="n",
  # or do log of y axis, delete # symbol
  log="y"
)

#axis function to write y axis with only scale by 10s. 
#USE YOUR OWN APPROPRIATE NUMBERS
axis(
  #2 is to specify left axis (aka y axis)
  2,
  at=c(0.0001,0.001,0.01,0.1)
  # dont need to define labels if same as numbers on axis
  #labels=aty
  #tck marks
  #tck=-0.01
)

#add legend in top right corner
legend("topright", 
       #inset legend off from border
       inset= 0.01,
       #names of each category based on factors, alphabetical order of category 1-5 of TypeOfSites
       legend = levels(data$TypeOfSite), 
       #symbols matching dataframe's factors 1:5 of data$TypeOfSite
       pch=21,
       #colors matching dataframe's factors 1:5 of data$TypeOfSite
       col="black",
       #fill colors of circle matching points of plot
       pt.bg=c(1:5),
       #specify scale of whole box of legend to not block data
       cex=.75,
       #specify point size in legend
       pt.cex=3,
       #remove legend border if "n" is specificed. "o" displays border
       bty="o",
       #specific box border thickness/width
       box.lwd=2,
       #specify box border type
       #box.lty=,
       #text.width change
       text.width=10,
       #justify text legend, xjust=0 is left justified, xjust=1 means right justified
       xjust=1
)

#Load/Save Workspace 

save(df,file="df.Rda")
load("df.Rda")
