




library(seqinr)

read.fasta("DengueViruses/DengueVirus1.fasta_pruned.mu.trim05") -> VIRUS

# reference sequence

ref <- VIRUS[[1]]

length(ref)->Seqlenght


# dataframe columns

num <- c(1:Seqlenght)

wtnt <- c()

freq <- c()

# for freq calculation later

absfreq <- c(rep(0,Seqlenght))

totalcount <- c(rep(0,Seqlenght))

# average WT calculation

# counts number of each nucleotide in each position

acount <- c(rep(0,Seqlenght))

gcount <- c(rep(0,Seqlenght))

ccount <- c(rep(0,Seqlenght))

tcount <- c(rep(0,Seqlenght))

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


# creates dataframe containing all data

VIRUS_DATA <- data.frame(num,wtnt,freq)

View(VIRUS_DATA)


# Investigating CPG


# Collapses data into a string and sets to variable STRING
paste(VIRUS_DATA$wtnt, collapse = '') -> STRING
STRING

#Looks for pattern tg in the data STRING, sets location sites to variable TG
gregexpr(pattern ='tg',STRING ) -> TG
TG

# Insert LIST TG into a data frame BELL
BELL <- data.frame(matrix(unlist(TG)))
BELL

# create new column name CPG with values zero
VIRUS_DATA$CPG<- 0

# Inserting value 1 in every TG site using BELL
VIRUS_DATA[BELL[,1],"CPG"] <- 1

View(VIRUS_DATA)

# For CA sites: 

#Looks for pattern ca in the data STRING, sets location sites to variable CA
gregexpr(pattern ='ca',STRING ) -> CA
CA

# Insert LIST CA into a data frame CASITES
CASITES <- data.frame(matrix(unlist(CA)))
CASITES

# Inserting value 1 in every CA site using CA
VIRUS_DATA[CASITES[,1]+1,"CPG"] <- 1

View(VIRUS_DATA)

#Hasan's code was updated 1 day ago
#setting up enviroment

library(seqinr)
library(Biostrings)
library(ape)

#reading files and setting data  
seqs = read.dna("InfluenzaAvirus_HA_H1N1.fasta.mu.fasta", format = "fasta", as.character=TRUE)
a=nrow(VIRUS_DATA)
b=ncol(VIRUS_DATA)


# Let's establish the new data frame
# Number of columns is arbitrary. I just want to print the WT sequence down a column
VIRUS_DATA=data.frame(matrix(nrow=b, ncol=7)) 
names(VIRUS_DATA)=c("WTseq", "Mutseq", "WTAA", "MutAA", "WTcat", "Mutcat", "DrasticAA")

#mean frequency


# We're defining the wild type sequence as the most common sequence
for(i in 0:(b-1)){
    x=table(VIRUS_DATA[,i+1])
    y=names(x[which.max(x)])
    VIRUS_DATA[i+1,1]=y
}

# Translate the WT DNA sequence to AA sequence

curSeq <- translate(paste(VIRUS_DATA[,1], sep=" "), NAstring="X", ambiguous = FALSE, sens="F")
count <- 0

for(i in 1:length(curSeq)){ # incrementing down sequence by 3 (needs work)
    count <- count + 1
    VIRUS_DATA[count,]$WTAA <- curSeq[i]
    count <- count + 1
    VIRUS_DATA[count,]$WTAA <- curSeq[i]
    count <- count + 1
    VIRUS_DATA[count,]$WTAA <- curSeq[i]
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
    VIRUS_DATA[j,5]=amCat(VIRUS_DATA[j,3])
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

# Let's Creates 
df$MutAA <-c(0)

count <- 1
i = 1
j = 1
for(i in 1:(nrow(df)/3)){ # Workspace for MutAA filling don't run this loop
    
    for (j in 1:3){
        
        if(j == 1){                     #first codon
            
            df[count, ]$MutAA <- translate(firstC <- c(as.character(df[count,]$Mutseq), as.character(df[count + 1,]$WTseq), as.character(df[count +2,]$WTseq)))
        } 
        if(j == 2){                #second codon
            df[count, ]$MutAA <- translate(secondC <- c(as.character(df[count - 1,]$WTseq), as.character(df[count,]$Mutseq), as.character(df[count +1,]$WTseq)))
        }
        if(j == 3){                  #third codon
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
df <- DrasticChange(df)
save(df,file="df.Rda")
load("df.Rda")