setwd("/Users/shantothenel/Desktop/HumanBocaVirus")

library(seqinr)
read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.txt") -> VIRUS


#THIS IS NOT MY CODE -- THIS IS GORDON'S CODE

# reference sequence
ref <- VIRUS[[1]]
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
View(VIRUS_DATA)


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

#making a dataframe out of all the data compiled -- position num, wtnt, mutnt, mean freq, wtAA, and mutAA
VIRUS_DATA <- data.frame(num,wtnt,mutnt,freq)
View(VIRUS_DATA)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# -- note for 10/31 --
#need to try out drastic/ nondrastic
#also need to fix the mutated AA and wt AA somehow -- get victoria's or avery's code?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#THIS IS NOT MY CODE -- THIS IS AVERY'S CODE THANKS AVERY!!

#setting up enviroment
library(seqinr)
library(Biostrings) #i was not able to get this to download
library(ape)

#reading files and setting data 

seqs = read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.txt", format = "fasta", as.character=TRUE)
a=nrow(seqs)
b=ncol(seqs)


# Let's establish the new data frame
# Number of columns is arbitrary. I just want to print the WT sequence down a column
VIRUS_DATA2=data.frame(matrix(nrow=b, ncol=7)) 
names(VIRUS_DATA2)=c("WTseq", "Mutseq", "WTAA", "MutAA", "WTcat", "Mutcat", "DrasticAA")

#mean frequency
# We're defining the wild type sequence as the most common sequence
for(i in 0:(b-1)){
    x=table(seqs[,i+1])
    y=names(x[which.max(x)])
    VIRUS_DATA2[i+1,1]=y
}

# Translate the WT DNA sequence to AA sequence
curSeq <- translate(paste(VIRUS_DATA2[,1], sep=" "), NAstring="X", ambiguous = FALSE, sens="F")
count <- 0

for(i in 1:length(curSeq)){ # incrementing down sequence by 3 (needs work)
    count <- count + 1
    VIRUS_DATA2[count,]$WTAA <- curSeq[i]
    count <- count + 1
    VIRUS_DATA2[count,]$WTAA <- curSeq[i]
    count <- count + 1
    VIRUS_DATA2[count,]$WTAA <- curSeq[i]
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

# categorizing the WT amino acids
for(j in 1:b){
    VIRUS_DATA2[j,5]=amCat(VIRUS_DATA2[j,3])
}

# creates Mut strain of current WTseq enters it into 
for (i in 1:b){
    curNuc <-  VIRUS_DATA2[i,1]
    if(curNuc == "a" || curNuc == "g"){
        if(curNuc == "a"){
            VIRUS_DATA2 [i, 2] = "g"
        }else{
            VIRUS_DATA2 [i, 2] = "a"
        }
    }else{
        if(curNuc == "t"){
            VIRUS_DATA2 [i, 2] = "c"
        }else{
            VIRUS_DATA2 [i, 2] = "t"
        }
    }
}



# Let's Creates 

VIRUS_DATA2$MutAA <-c(0)

count <- 1
i = 1
j = 1
for(i in 1:(nrow(VIRUS_DATA2)/3)){ # Workspace for MutAA filling don't run this loop
    for (j in 1:3){
        if(j == 1){                     #first codon
            VIRUS_DATA2[count, ]$MutAA <- translate(firstC <- c(as.character(VIRUS_DATA2[count,]$Mutseq), as.character(VIRUS_DATA2[count + 1,]$WTseq), as.character(VIRUS_DATA2[count +2,]$WTseq)))
        } 
        if(j == 2){                #second codon
            VIRUS_DATA2[count, ]$MutAA <- translate(secondC <- c(as.character(VIRUS_DATA2[count - 1,]$WTseq), as.character(VIRUS_DATA2[count,]$Mutseq), as.character(VIRUS_DATA2[count +1,]$WTseq)))
        }
        if(j == 3){          #third codon
            VIRUS_DATA2[count, ]$MutAA <- translate(thirdC <- c(as.character(VIRUS_DATA2[count - 2,]$WTseq), as.character(VIRUS_DATA2[count - 1,]$WTseq), as.character(VIRUS_DATA2[count,]$Mutseq)))
        }
        count <- count + 1
        print(j)
        print(i)
        print(count)
    }
}


VIRUS_DATA2[,4]=seqinr::translate(paste(VIRUS_DATA2[,2], sep=" "),
                         NAstring="X", ambiguous=FALSE, sens="F")            # Translating from nucleic acid to AA

for(j in 1:b){                                              
    VIRUS_DATA2[j,6]=amCat(VIRUS_DATA2[j,4])                                    # Categorizing the AA 
}


DrasticChange <- function(VIRUS_DATA2){
    for(h in 1:nrow(VIRUS_DATA2)){
        if (VIRUS_DATA2[h,]$WTcat==VIRUS_DATA2[h,]$Mutcat){
            VIRUS_DATA2[h,]$DrasticAA= 0                      # If WT AA category = Mut AA Category, no drastic change
        }
        if (VIRUS_DATA2[h,]$WTcat!=VIRUS_DATA2[h,]$Mutcat){
            VIRUS_DATA2[h,]$DrasticAA= 1                      # If WT AA category =/= Mut AA Category, yes drastic change
        }
    }
    return(VIRUS_DATA2)
}

VIRUS_DATA2 <- DrasticChange(VIRUS_DATA2)
View(VIRUS_DATA2)

VIRUS_DATA$WTAA <- VIRUS_DATA2$WTAA
head(VIRUS_DATA,10)

VIRUS_DATA$MutAA <- VIRUS_DATA2$MutAA
head(VIRUS_DATA,10)

VIRUS_DATA$WTcat <- VIRUS_DATA2$WTcat
head(VIRUS_DATA,10)

VIRUS_DATA$Mutcat <- VIRUS_DATA2$Mutcat
head(VIRUS_DATA,10)

#assuming we have WTAA and MutAA already
#getting the syn/ nonsyn/ nonsense
for (h in 1:nrow(VIRUS_DATA)){
    if(VIRUS_DATA[h,"MutAA"]== VIRUS_DATA[h,"WTAA"]){
        VIRUS_DATA[h,"TypeOfSite"] = "synonymous"
    }
    if(VIRUS_DATA[h,"MutAA"] != VIRUS_DATA[h,"WTAA"]){
        if(VIRUS_DATA[h,"MutAA"]=="*"){
            VIRUS_DATA[h,"TypeOfSite"] = "nonsense"}
        else {
            VIRUS_DATA[h,"TypeOfSite"] = "nonsynonymous"
        }
    }
}

View(VIRUS_DATA)

VIRUS_DATA$DrasticAA <- VIRUS_DATA2$DrasticAA
head(VIRUS_DATA,10)

View(VIRUS_DATA)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

final_BOCA <- VIRUS_DATA

final_BOCA
save(final_BOCA,file="final_BOCAdf.Rda")
load("final_BOCAdf.Rda")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# MAKE THE DANG THING #


#input should be dataframe?
#the command should be something like the function pulls out information
    #from the dataframe and uses that info to plot
#output should be plot of Meanfreq vs. Position (num), colored by TypeofSite
virus_plot = function(df) {
    {
    x <- df$num
    y <- df$freq
    v_plot <- plot(x,y, type = "p", main = "Nucleotide Positon (num) vs. Mean Frequency (freq) of Virus", 
         log = "y", xlab = "num", ylab = "log of freq", col = df$TypeOfSite)
    }
    print(v_plot)
    return(v_plot)
}

virus_plot(VIRUS_DATA)


#skeleton of the function -- so far i have the colors for the syn/ nonsyn
plot(hiv_csv$num, hiv_csv$MeanFreq, main = "Plot of Mean Frequency", log = "y", xlab = "Number", ylab = "log of Mean Frequency", col = hiv_csv$TypeOfSite)

legend("topleft",
       c("Syn", "Nonsyn", "Overlap", "Res", "Stop"),
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "grey", "purple", "black"),
       bty = "n")










