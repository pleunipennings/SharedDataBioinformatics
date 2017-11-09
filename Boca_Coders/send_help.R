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


#subset to only have open reading frame
#VIRUS_DATA <- VIRUS_DATA[2981:4996,]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#read Boca virus file 
setwd("/Users/shantothenel/Desktop/HumanBocaVirus")

#boca NS1 stands for bocavirus NS1 gene
bocaNS1seqs<-read.fasta("HumanBocavirus1_VP1.fasta_pruned.mu.txt") 
#How to get the consensus (= most common nucleotide for each position)
bocaNS1seqsAli<-read.alignment("HumanBocavirus1_VP1.fasta_pruned.mu.txt", format="fasta")
#works when using con() instead of consensus()
seqinr::consensus(bocaNS1seqsAli)->cons
#basic data frame with mean freq
V_DATA<-data.frame("num"=c(1:ncol(bocaNS1seqsDNA)),"wtnt" = cons)
V_DATA
#use read.dna to get the data in matrix form, this makes it easier to count
bocaNS1seqsDNA<-read.dna("HumanBocavirus1_VP1.fasta_pruned.mu.txt", format = "fasta",as.character=TRUE)


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

for(x in seq(1, length(cons) - 2, 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    new_AA <- translate(codon)
    WTAA[x] <- new_AA
    WTAA[x+1] <- new_AA
    WTAA[x+2] <- new_AA
}
#Insert WTaa and MTaa into dataframe
WTAA->VIRUS_DATA$WTAA
MUTAA->VIRUS_DATA$MUTAA

head(VIRUS_DATA)

vp1_safe <- VIRUS_DATA
head(VIRUS_DATA)
#vp1_safe -> VIRUS_DATA

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


#  Note that this depends on the columns MUTAA and WTAA being already complete!

functionSynNonSyn<-function(DF){
    DF$TypeofSite<-c()
    for (h in 1:nrow(DF)){
        if(DF[h,"MUTAA"]== DF[h,"WTAA"]){
            DF[h,"TypeOfSite"] = "syn"
        }
        if(DF[h,"MUTAA"] != DF[h,"WTAA"]){
            if(DF[h,"MUTAA"]=="*"){
                DF[h,"TypeOfSite"] = "nonsense"}
            else {
                DF[h,"TypeOfSite"] = "nonsyn"
            }
        }
    }
    DF$TypeOfSite<-as.factor(DF$TypeOfSite)
}

functionSynNonSyn(VIRUS_DATA) -> VIRUS_DATA$TypeOfSite
head(VIRUS_DATA)

vp1_safe <- VIRUS_DATA
head(VIRUS_DATA)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#THIS IS NOT MY CODE -- THIS IS AVERY'S CODE THANKS AVERY!!

VIRUS_DATA$WTcat <- 0
VIRUS_DATA$MUTcat <- 0

drasticAA<-function(DF){
    #Amino Acid Changes 
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
    
    #for loop to categorize WTAA
    for(i in 1:nrow(DF)){                                              
        DF[i,8]=amCat(DF[i,5])                                    
    }
    
    #for loop to categorize MUTAA
    for(i in 1:nrow(DF)){                                              
        DF[i,9]=amCat(DF[i,6])                                 
    }
    
    
    #for loop to determine drastic change or not
    for(i in 1:nrow(DF)){
        if (DF$WTAAcat[i]==DF$MUTAAcat[i]){
            DF$bigAAchange[i]= "conserved"
        }
        if (DF$WTAAcat[i]!=DF$MUTAAcat[i]){
            DF$bigAAchange[i] = "drastic"
        }
    }
}
drasticAA(VIRUS_DATA) ->VIRUS_DATA$bigAAchange







