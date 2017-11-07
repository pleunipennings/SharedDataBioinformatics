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

#making a dataframe out of all the data compiled -- position num, wtnt, mutnt, mean freq, wtAA, and MUTAA
VIRUS_DATA <- data.frame(num,wtnt,mutnt,freq)
View(VIRUS_DATA)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# -- note for 10/31 --
#need to try out drastic/ nondrastic
#also need to fix the mutated AA and wt AA somehow -- get victoria's or avery's code?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# THIS IS NOT MY CODE, THIS IS SCOTT'S CODE FOR THE WTAA AND MUTAA !!!

transition_nuc = function(nuc) {
    # %in% means if you find any of these things in the list, return 
    # what I want anyways
    if(nuc %in% c("A","a")) {
        mutatedNuc = "G" 
    } else if(nuc %in% c("G","g")) {
        mutatedNuc = "A" 
    } else if(nuc %in% c("T","t")) {
        mutatedNuc = "C"
    } else if(nuc %in% c("C","c")) {
        mutatedNuc = "T"
    } else if(nuc != T) {
        mutatedNuc = "ERROR"
    }
    return(mutatedNuc)
}

for (i in seq.int(1,nrow(VIRUS_DATA),3)) {
    codon <- c(VIRUS_DATA$wtnt[i],VIRUS_DATA$wtnt[i+1],VIRUS_DATA$wtnt[i+2])
    VIRUS_DATA$WTAA[i] = seqinr::translate(codon)
    VIRUS_DATA$WTAA[i+1] = seqinr::translate(codon)
    VIRUS_DATA$WTAA[i+2] = seqinr::translate(codon)
    copy = codon
    copy[1] = transition_nuc(codon[1])
    VIRUS_DATA$MUTAA[i] = seqinr::translate(copy)
    copy = codon
    copy[2] = transition_nuc(codon[2])
    VIRUS_DATA$MUTAA[i+1] = seqinr::translate(copy)
    copy = codon
    copy[3] = transition_nuc(codon[3])
    VIRUS_DATA$WTAA[i+2] = seqinr::translate(copy)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




#THIS IS NOT MY CODE -- THIS IS AVERY'S CODE THANKS AVERY!!

#setting up enviroment
library(seqinr)
#library(Biostrings) ; i was not able to get this to download
library(ape)

#reading files and setting data 

seqs = read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.txt", format = "fasta", as.character=TRUE)
a=nrow(seqs)
b=ncol(seqs)


# Let's establish the new data frame
# Number of columns is arbitrary. I just want to print the WT sequence down a column
VIRUS_DATA2=data.frame(matrix(nrow=b, ncol=7)) 
names(VIRUS_DATA2)=c("WTseq", "Mutseq", "WTAA", "MUTAA", "WTcat", "Mutcat", "bigAAchange")

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

VIRUS_DATA2$MUTAA <-c(0)

count <- 1
i = 1
j = 1
for(i in 1:(nrow(VIRUS_DATA2)/3)){ # Workspace for MUTAA filling don't run this loop
    for (j in 1:3){
        if(j == 1){                     #first codon
            VIRUS_DATA2[count, ]$MUTAA <- translate(firstC <- c(as.character(VIRUS_DATA2[count,]$Mutseq), as.character(VIRUS_DATA2[count + 1,]$WTseq), as.character(VIRUS_DATA2[count +2,]$WTseq)))
        } 
        if(j == 2){                #second codon
            VIRUS_DATA2[count, ]$MUTAA <- translate(secondC <- c(as.character(VIRUS_DATA2[count - 1,]$WTseq), as.character(VIRUS_DATA2[count,]$Mutseq), as.character(VIRUS_DATA2[count +1,]$WTseq)))
        }
        if(j == 3){          #third codon
            VIRUS_DATA2[count, ]$MUTAA <- translate(thirdC <- c(as.character(VIRUS_DATA2[count - 2,]$WTseq), as.character(VIRUS_DATA2[count - 1,]$WTseq), as.character(VIRUS_DATA2[count,]$Mutseq)))
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




VIRUS_DATA2 <- DrasticChange(VIRUS_DATA2)
View(VIRUS_DATA2)

VIRUS_DATA$WTAA <- VIRUS_DATA2$WTAA
head(VIRUS_DATA,10)

VIRUS_DATA$MUTAA <- VIRUS_DATA2$MUTAA
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


# # # # # # # # # # # # # # CPG SECTION # # # # # # # # # # # # # # # #

setwd("/Users/shantothenel/Desktop/HumanBocaVirus")
library(seqinr)
library(ape)

#bi2sfsu@gmail.com
#bi2sfsu!

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_NS1.fasta_pruned.mu.txt") 
#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_NS1.fasta_pruned.mu.txt", format="fasta")
#works when using con() instead of consensus()
seqinr::consensus(bocaNS1seqsAli)->WTnt
#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_NS1.fasta_pruned.mu.txt", format = "fasta",as.character=TRUE)
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


VIRUS_DATA$makesCpG <- BoNS1df$makesCpG
View(VIRUS_DATA)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

final_BOCA01 <- VIRUS_DATA
final_BOCA_NS1 <- final_BOCA01
final_BOCA_NS1

head(final_BOCA_NS1,10)
save(final_BOCA_NS1,file="final_BOCAdf_NS1.Rda")
load("final_BOCAdf_NS1.Rda")

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

virus_plot(final_BOCA_NS1) #THIS FINALLY WORKED???

#making a legend
legend("topleft",
       c("syn", "nonsyn", "nonsense"),
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")


which(final_BOCA_NS1$TypeOfSite=="syn") #279 positions
which(final_BOCA_NS1$TypeOfSite=="nonsyn") #4762 positions -> why is it mostly red
which(final_BOCA_NS1$TypeOfSite=="nonsense") #116 positions

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#TRYING TO SAVE THE PLOT AS A PDF NOW
  #this is for reference
pdf(file="Trial One's Dropout at One Allele.pdf")
plot(oneA, main = "Trial One's Dropout at One Allele", xlab="Profile Matches", ylab="Alleles", col="blue")
abline(v = 25, col = "red") #indicates number of alleles where the partial match was highest
dev.off()

  #this is for BOCA
pdf(file="Nucleotide Position (num) vs. Mean Frequency (freq) of Boca NS1")
virus_plot(final_BOCA_NS1)
legend("topleft",
       c("syn", "nonsyn", "nonsense"),
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#TRIED THE CPG FIGURE

read.csv("HumanBocavirus1_NS1.fasta_pruned.mu.txt")-> BOCAdata
print(BOCAdata)

#Identify noCpG sites
which(BOCAdata$makesCpG=="0")->noCpG
print(noCpG)
#Identify CpG sites
which(BOCAdata$makesCpG=="1")->CpG
CpG

#Identify mean frequency in association with noCpG sites
BOCAdata[noCpG,"MeanFreq"]->Freq
Freq
#Identify mean frequency in association with CpG sites
BOCAdata[CpG,"MeanFreq"]->Freq2
Freq2



#Making a plot of mean frequency vs. CpG/nonCpG location
plot(noCpG, Freq+0.0001, ylim=c(0.0001,0.5), xlim = c(0.1, 1), log = "y", col="black",pch=21.25, bg=rgb(1,0,0,0.5), main="CpG/noCpG location vs. Mean Frequency", xlab = "Location", ylab = "Mean Frequency")

#Overlapping the two graphs
points(CpG, Freq2+0.0001, col="black",pch=21.25, bg=rgb(0,0,1,0.5),xaxt='n',yaxt='n', ann = FALSE)



#legend with circle instead of line  
legend("topleft",c("noCpG","CpG"),cex=0.6, 
       col ="black",bg=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),bty="n",
       title="Legend",inset=.02,pch =21,pt.bg = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))












