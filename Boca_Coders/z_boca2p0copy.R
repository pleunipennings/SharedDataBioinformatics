setwd("/Users/shantothenel/Desktop/HumanBocaVirus")
  
library(seqinr)
read.fasta("HumanBocavirus1_VP1.fasta_pruned.mu.txt") -> VIRUS


#THIS IS NOT MY CODE -- THIS IS GORDON'S CODE

# reference sequence
ref <- VIRUS[[1]]
ref
length(ref)->Seqlength

# dataframe columns
num <- c(1:Seqlength)
wtnt <- c()
freq <- c()

# for freq calculation later
absfreq <- c(rep(0,Seqlength))
totalcount <- c(rep(0,Seqlength))

# average WT calculation
# counts number of each nucleotide in each position
acount <- c(rep(0,Seqlength))
gcount <- c(rep(0,Seqlength))
ccount <- c(rep(0,Seqlength))
tcount <- c(rep(0,Seqlength))
nuc <- c()


# same as line 20 comment
for (i in 1:length(VIRUS)) {
  sequence <- VIRUS[[i]]
  for (j in 1:length(sequence)) {
    if (sequence[j] == 'a') {acount[j] = acount[j] + 1}
    if (sequence[j] == 'g') {gcount[j] = gcount[j] + 1}
    if (sequence[j] == 'c') {ccount[j] = ccount[j] + 1}
    if (sequence[j] == 't') {tcount[j] = tcount[j] + 1}
  }
}


# assigns wtnt based on most frequent nucleotide across all sequences per position
for (j in 1:length(sequence)) {
  nuc[j] <- max(c(acount[j],gcount[j],ccount[j],tcount[j]))
  if (max(nuc[j]) == acount[j]) {wtnt[j] <- 'a'}
  if (max(nuc[j]) == gcount[j]) {wtnt[j] <- 'g'}
  if (max(nuc[j]) == ccount[j]) {wtnt[j] <- 'c'}
  if (max(nuc[j]) == tcount[j]) {wtnt[j] <- 't'}
  nuc <- c()
}


# gives absolute totals to be used for frequency calculation
for (i in 1:length(VIRUS)) {
  sequence <- VIRUS[[i]]
  for (j in 1:length(sequence)){
    if (wtnt[j] == 'a') {
      if (sequence[j] == 'g') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 'a') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
    if (wtnt[j] == 'g') {
      if (sequence[j] == 'a') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 'g') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
    if (wtnt[j] == 'c') {
      if (sequence[j] == 't') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 'c') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
    if (wtnt[j] == 't') {
      if (sequence[j] == 'c') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 't') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
  }
}




# calculates frequency as percentage

for (i in 1:length(absfreq)) {
  freq[i] <- absfreq[i] / totalcount[i]
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


#making the dataframe so far 
VIRUS_DATA <- data.frame(num,wtnt,freq)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#making the transition function, based off the nucleotide transition functions
transition_vec = function(vector) {
  # %in% means if you find any of these things in the list, return 
  # what I want anyways
  for(i in 1:length(vector))
    if(vector[i] %in% c("A","a")) {
      vector[i] = "g" 
    } else if(vector[i] %in% c("G","g")) {
      vector[i] = "a" 
    } else if(vector[i] %in% c("T","t")) {
      vector[i] = "c"
    } else if(vector[i] %in% c("C","c")) {
      vector[i] = "t"
    } else if(vector[i] != T) {
      vector[i] = "ERROR"
    }
  return(vector)
}


#mutating the wtnt and making it into a vector, then assigning it to a variable
mutnt = transition_vec(wtnt)

#making a dataframe out of all the data compiled -- position num, wtnt, mutnt, mean freq, wtAA, and MUTAA
VIRUS_DATA <- data.frame(num,wtnt,mutnt,freq)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# -- note for 10/31 --
#need to try out drastic/ nondrastic
#also need to fix the mutated AA and wt AA somehow -- get victoria's or avery's code?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#read Boca virus file 
setwd("/Users/shantothenel/Desktop/HumanBocaVirus")

library(seqinr)
read.fasta("HumanBocavirus1_VP1.fasta_pruned.mu.txt") -> VIRUS
#save Boca file as a variable
library(ape)
WTAA<-c()
MUTAA<-c()
#for loop for transition mutation at each position and translating it 
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

for(x in seq(1, length(cons)-2, 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    new_AA <- translate(codon)
    WTAA[x] <- new_AA
    WTAA[x+1] <- new_AA
    WTAA[x+2] <- new_AA
}
#Insert WTaa and MTaa into dataframe
WTAA->VIRUS_DATA$WTAA
MUTAA->VIRUS_DATA$MUTAA


#subset to only have open reading frame
VIRUS_DATA <- VIRUS_DATA[2981:4996,]

###Get Wildtype amino acid###

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


#THIS IS NOT MY CODE -- THIS IS AVERY'S CODE THANKS AVERY!!

library(seqinr)
#library(Biostrings)
library(ape)

#reading files and setting data  
seqs = read.dna("HumanBocavirus1_VP1.fasta_pruned.mu.txt", format = "fasta", as.character=TRUE)
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


VIRUS_DATA2 <- DrasticChange(VIRUS_DATA2[2981:4997,])
View(VIRUS_DATA2)

VIRUS_DATA$WTAA <- VIRUS_DATA2$WTAA
head(VIRUS_DATA,10)

VIRUS_DATA$MUTAA <- VIRUS_DATA2$MutAA
head(VIRUS_DATA,10)

VIRUS_DATA$WTcat <- VIRUS_DATA2$WTcat
head(VIRUS_DATA,10)

VIRUS_DATA$Mutcat <- VIRUS_DATA2$Mutcat
head(VIRUS_DATA,10)

#assuming we have WTAA and MUTAA already
#getting the syn/ nonsyn/ nonsense
for (h in 1:nrow(VIRUS_DATA)){
  if(VIRUS_DATA[h,"MUTAA"]== VIRUS_DATA[h,"WTAA"]){
    VIRUS_DATA[h,"TypeOfSite"] = "syn"
  }
  if(VIRUS_DATA[h,"MUTAA"] != VIRUS_DATA[h,"WTAA"]){
    if(VIRUS_DATA[h,"MUTAA"]=="*"){
      VIRUS_DATA[h,"TypeOfSite"] = "nonsense"}
    else {
      VIRUS_DATA[h,"TypeOfSite"] = "nonsyn"
    }
  }
}

View(VIRUS_DATA)

VIRUS_DATA$bigAAchange <- VIRUS_DATA2$bigAAchange
head(VIRUS_DATA,10)

safe <- VIRUS_DATA


# # # # # # # # # # # # # # CPG SECTION # # # # # # # # # # # # # # # #

setwd("/Users/shantothenel/Desktop/HumanBocaVirus")
library(seqinr)
library(ape)

#bi2sfsu@gmail.com
#bi2sfsu!

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_VP1.fasta_pruned.mu.txt") 
#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_VP1.fasta_pruned.mu.txt", format="fasta")
#works when using con() instead of consensus()
seqinr::consensus(bocaNS1seqsAli)->WTnt
#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_VP1.fasta_pruned.mu.txt", format = "fasta",as.character=TRUE)
# Writes new function to create transition mutations in nucleotides
transition <- function(nt){
  if(nt=="a") {return("g")}
  if(nt=="g") {return("a")}
  if(nt=="c") {return("t")}
  if(nt=="t") {return("c")}}
# For-loop to calculate mean frequency of transition mutations for each nucleotide:
MeanFreq<-c()
for (i in 1:ncol(bocaNS1seqsDNA)){
  MeanFreq<-c(MeanFreq,(length(which(bocaNS1seqsDNA[,i]==transition(WTnt[i])))/ncol(bocaNS1seqsDNA)))}
BoNS1df<-data.frame("num"=c(1:ncol(bocaNS1seqsDNA)),WTnt,MeanFreq)

#CPG sites input
CpG_finder <- function(new_virus_data){
  #reads data into function as CSV file
  virus_data <- new_virus_data
  
  #singles out column of nucleotides -- this assumes that the new datafile has the same column headers such as WTnt as the HIV file does
  WTnt <- virus_data$WTnt
  
  #gets the length of the data file for use in a loop
  data_length <- nrow(virus_data)
  
  #creates an empty vector (like a list) of the same length as the data file to be used to record the results of the loop below
  makesCpG <- vector(mode = "numeric", length = data_length)
  
  #loop that determines if a CpG could occur due to mutation at each spot in the list of nucleotides (WTnt)
  #loops from row 1 to the the last row of the column WTnt
  for(x in 1:data_length){
    
    #assigns a name (current_nuc) to the nucleotide at row x in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
    current_nuc <- toupper(WTnt[x])
    
    #assigns a name (current_neighbor) to the nucleotide in the next row down in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
    current_neighbor <- toupper(WTnt[x + 1])
    
    #begins a loop that compares the two nucleotides that were just isolated
    #this if statement is only evaluated if the current nucleotide is a T
    if(current_nuc == "T"){
      #this if statement is only evaluated if the neighbor is a G
      if(current_neighbor == "G"){
        #if the current nucleotide is a T and the neighbor is a G (thus a mutation can cause a CpG), then the of the current nucleotide in the list of 0s we made is changed to a 1
        makesCpG[x] = 1
      }
      #otherwise, nothing is changed
      else{}
    }
    #repeat the same loop as before, but modified slightly to look for CA pairs
    if(current_nuc == "C"){
      if(current_neighbor == "A"){
        #here, instead of changing the position of the current nuc to a 1, we change the neighbor's position because that is where the mutation would be, hence x+1
        makesCpG[x+1] = 1
      }
      else{}
    }
    else{}
    #print(c("Data length is ", nrow(virus_data), "List length is ", length(makesCpG)))
  }
  
  #append the list of 0s and 1s we created to the end of the data file we imported
  virus_data$makesCpG <- makesCpG
  
  #when the entire function is done, return the amended data file
  return(virus_data)
}

CpG_finder(BoNS1df)->BoNS1df
head(BoNS1df, 100)


VIRUS_DATA$makesCpG <- BoNS1df[2981:5024,]$makesCpG
View(VIRUS_DATA)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

final_BOCA_VP1 <- VIRUS_DATA
View(final_BOCA_VP1)

head(final_BOCA_VP1,10)
save(final_BOCA_VP1,file="final_BOCAdf_VP1.Rda")
load("final_BOCAdf_VP1.Rda")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# MAKE THE DANG THING #

#library(ggplot2) -- not sure if we need it yet; was playing around with it

#input should be dataframe?
#the command should be something like the function pulls out information
#from the dataframe and uses that info to plot
#output should be plot of Meanfreq vs. Position (num), colored by TypeofSite
virus_plot = function(df) {
  plot(df$num, df$freq+0.1, type = "p", main = "Nucleotide Positon (num) vs. Mean Frequency (freq) of Virus", 
       log = "y", xlab = "num", ylab = "log of freq", col = as.factor(df$TypeOfSite))
}

virus_plot(final_BOCA_VP1) #THIS FINALLY WORKED???

#making a legend
legend("topleft",
       c("syn", "nonsyn", "nonsense"),
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")


which(final_BOCA_VP1$TypeOfSite=="syn") #102 positions
which(final_BOCA_VP1$TypeOfSite=="nonsyn") #1897 positions -> why is it mostly red
which(final_BOCA_VP1$TypeOfSite=="nonsense") #45 positions

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#TRYING TO SAVE THE PLOT AS A PDF NOW
#this is for reference
pdf(file="Trial One's Dropout at One Allele.pdf")
plot(oneA, main = "Trial One's Dropout at One Allele", xlab="Profile Matches", ylab="Alleles", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest
dev.off()

#this is for BOCA
pdf(file="Nuc Position (num) vs. Mean Frequency (freq) of Boca VP1")
virus_plot(final_BOCA_VP1)
legend("topleft",
       c("syn", "nonsyn", "nonsense"),
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#TRIED THE CPG FIGURE

read.csv("HumanBocavirus1_VP1.fasta_pruned.mu.txt")-> BOCAdata
print(BOCAdata)

VIRUS_DATA


#This function makes a plot that shows the location of CpG and noCpG sites vs. Frequency.
#This function assumes that you have num, freq, and CPG are already in your dataframe.

#Team: Christen, Rima, Nicole, and Kellen.
#Would also like to thank Scott and Pleuni for their help with making corrections to our code. 

#Identify noCpG sites
which(VIRUS_DATA$makesCpG=="0")->noCpG
print(noCpG)
#Identify CpG sites
which(VIRUS_DATA$makesCpG=="1")->CpG
CpG

####### beginning of JACKY'S CODE #######

VIRUS_DATA[noCpG,1] -> NEWnoCpG
NEWnoCpG
VIRUS_DATA[CpG, 1] -> NEWCpG
NEWCpG

#Identify mean frequency in association with noCpG sites
VIRUS_DATA[noCpG,"freq"]->freq
#Identify mean frequency in association with CpG sites
VIRUS_DATA[CpG,"freq"]->freq2

plot(NEWnoCpG, freq+0.0001, ylim=c(0.001,0.3), col="black",pch=21.25, bg=rgb(1,0,0,0.5), main="CpG/noCpG location vs. Mean Frequency", xlab = "Location", ylab = "Mean Frequency")
points(NEWCpG, freq2+0.0001, col="black",pch=21.25, bg=rgb(0,0,1,0.5))

legend("topleft",c("noCpG","CpG"),cex = 1,
       col ="black",bg=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),bty="n",
       title="Legend",inset=.02,pch =21,pt.bg = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))

####### ending of JACKY'S CODE #######













