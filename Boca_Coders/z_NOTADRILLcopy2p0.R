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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#subset to only have open reading frame
VIRUS_DATA <- VIRUS_DATA[178:2094,]

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




VIRUS_DATA2 <- DrasticChange(VIRUS_DATA2[178:2094,])
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
for (h in 1:nrow(VIRUS_DATA)){ #looks at each row in the dataframe df
    if(VIRUS_DATA[h,"MUTAA"]== VIRUS_DATA[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
        VIRUS_DATA[h,"TypeOfSite"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(VIRUS_DATA[h,"MUTAA"] != VIRUS_DATA[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
        if(VIRUS_DATA[h,"MUTAA"]=="*"){ #and it the MUTAA value equal "*"
            VIRUS_DATA[h,"TypeOfSite"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
        else {
            VIRUS_DATA[h,"TypeOfSite"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
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


VIRUS_DATA$makesCpG <- BoNS1df[178:2094,]$makesCpG
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


which(final_BOCA_NS1$TypeOfSite=="syn") #118 positions
which(final_BOCA_NS1$TypeOfSite=="nonsyn") #1755 positions -> why is it mostly red
which(final_BOCA_NS1$TypeOfSite=="nonsense") #44 positions

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


plot(NEWnoCpG, freq+0.0001, ylim=c(0.001,0.15), col="black",pch=21.25, bg=rgb(1,0,0,0.5), main="CpG/noCpG location vs. Mean Frequency", xlab = "Location", ylab = "Mean Frequency")
points(NEWCpG, freq2+0.0001, col="black",pch=21.25, bg=rgb(0,0,1,0.5))

legend("topleft",c("noCpG","CpG"),cex = 1,
       col ="black",bg=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),bty="n",
       title="Legend",inset=.02,pch =21,pt.bg = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))

####### ending of JACKY'S CODE #######


#Night_Crew Team
#Members Victoria, Sam, Gordon, Mordy 
#setwd("~/Desktop/Git")

#install.packages('ggplot2')
#install.packages('scales')
#scales is used to log the y axis in the graoh and ggplot2 is used to create the graph
library(scales)
library(ggplot2)

#this function graphs transition mutations separated by Synonymous Sites and Nonsynonymous sites
# this graph show the frequencies of the changes and also seperates combinations of drastic AA and CpG sites by color
#Start of the graph function, here data is your dataset. Make sure you have these columns and the column names are the same.bigAAChange, makesCpG/CPG, TypeOfSite, num, MeanFreqor/freq, wtnt/WTnt 
#In column TypeOfSite make sure that syn, nonsyn are present. No - or full names
comparing_mutation_graph = function(data){
    #building the combo lines to help sort (used later when info is subsetted) 
    #This will change your colomn name to what we orginally used to make our graph
#    colnames(data)[which(names(data) == "wtnt")] <- "WTnt"
#    colnames(data)[which(names(data) == "freq")] <- "MeanFreq"
#    colnames(data)[which(names(data) == "CPG")] <- "makesCpG"
#    data$combo<- (data$bigAAChange*3) + (data$makesCpG*2)
#    data$combo<- as.factor(data$combo)
#    levels(data$combo) <- gsub("0", "noAA noCPG", levels(data$combo))
#    levels(data$combo) <- gsub("2", "noAA yesCPG", levels(data$combo))
#    levels(data$combo) <- gsub("3", "yesAA noCPG", levels(data$combo))
#    levels(data$combo) <- gsub("5", "yesAA yesCPG", levels(data$combo))
#    #This will change your colomn name to what we orginally used to make our graph
#    colnames(data)[which(names(data) == "wtnt")] <- "WTnt"
#    colnames(data)[which(names(data) == "freq")] <- "MeanFreq"
    #Subsetting for the data that we really want just Synonymous and nonynonymous
    datatww <- subset(data, TypeOfSite=="syn" | TypeOfSite=="nonsyn")
    #if your dataframe has "-" try using the following line of code
    #datatww<-datatww[!(which(datatww$WTnt=="-")),]
    datatww$TypeOfSite <- ifelse(datatww$TypeOfSite == "syn", "Synonymous Sites", "Non-synonymous Sites")
    
    #gives values to each nucecotide and if they have a cetain combo
    #first need to introduce new columes for collecting the xvalue and the color
    datatww$xvalue<- 0 #giving each nucecotide a value depending on whether or not there is an AA change a CpG site and the nucecotide
    datatww$color <- 0 # four colors for the four combos: Blue = no bigAAChange/no makesCpG, Red = no bigAAChange/yes makesCpG, Green = yes bigAAChange/no makesCpG, Purple = yes bigAAChange/yes makesCpG 
    #for loop adds color and a xvalue to the dataset, for all nucecotides
    for (i in 1:length(datatww$num)) {
        if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "a") {
            datatww$xvalue[i] <- 1
            datatww$color[i] <- "blue"
        } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "a") {
            datatww$xvalue[i] <- 2
            datatww$color[i] <- "red"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "a") {
            datatww$xvalue[i] <- 3
            datatww$color[i] <- "green"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "a") {
            datatww$xvalue[i] <- 4
            datatww$color[i] <- "purple"
        } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "g") {
            datatww$xvalue[i] <- 5
            datatww$color[i] <- "blue"
        } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "g") {
            datatww$xvalue[i] <- 6
            datatww$color[i] <- "red"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "g") {
            datatww$xvalue[i] <- 7
            datatww$color[i] <- "green"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "g") {
            datatww$xvalue[i] <- 8
            datatww$color[i] <- "purple"
        } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "c") {
            datatww$xvalue[i] <- 9
            datatww$color[i] <- "blue"
        } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "c") {
            datatww$xvalue[i] <- 10 
            datatww$color[i] <- "red"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "c") {
            datatww$xvalue[i] <- 11
            datatww$color[i] <- "green"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "c") {
            datatww$xvalue[i] <- 12
            datatww$color[i] <- "purple"
        } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "t") {
            datatww$xvalue[i] <- 13
            datatww$color[i] <- "blue"
        } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "t") {
            datatww$xvalue[i] <- 14
            datatww$color[i] <- "red"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "t") {
            datatww$xvalue[i] <- 15
            datatww$color[i] <- "green"
        } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "t") {
            datatww$xvalue[i] <- 16
            datatww$color[i] <- "purple"
        }
    }
    #subsetiing all the data by amino acid, drastic AA, and CpG forming: Helps with analyst 
    #syn subset 
    syndata <- subset(datatww, TypeOfSite=="Synonymous Sites")
    #syn subset for no AA and no CPG
    synNNdata <- subset(syndata, combo=="noAA noCPG")
    #syn subset for no AA and no CPG for a, c, g, t (for HIV all 4 should be here)
    synNNa <- subset(synNNdata, WTnt=="a")
    synNNc <- subset(synNNdata, WTnt=="c")
    synNNg  <- subset(synNNdata, WTnt=="g")
    synNNt  <- subset(synNNdata, WTnt=="t")
    #syn subset for no AA and yes CPG
    synNYdata<- subset(syndata, combo=="noAA yesCPG")
    #syn subset for no AA and yes CPG for a, c, g, t (for HIV a, t )
    synNYa <- subset(synNYdata, WTnt=="a")
    synNYc <- subset(synNYdata, WTnt=="c")
    synNYg  <- subset(synNYdata, WTnt=="g")
    synNYt  <- subset(synNYdata, WTnt=="t")
    #syn subset for yes AA and no CPG 
    synYNdata <- subset(syndata, combo=="yesAA noCPG")
    #syn subset for yes AA and no CPG for a, c, g, t (none for HIV)
    synYNa <- subset(synYNdata, WTnt=="a")
    synYNc <- subset(synYNdata, WTnt=="c")
    synYNg  <- subset(synYNdata, WTnt=="g")
    synYNt  <- subset(synYNdata, WTnt=="t")
    #syn subset for yes AA and yes CPG
    synYYdata <- subset(syndata, combo=="yesAA yesCPG")
    #syn subset for yes AA and yes CPG for a, c, g, t (none for HIV)
    synYYa <- subset(synYYdata, WTnt=="a")
    synYYc <- subset(synYYdata, WTnt=="c")
    synYYg  <- subset(synYYdata, WTnt=="g")
    synYYt  <- subset(synYYdata, WTnt=="t")
    
    #nonsyn sub set 
    nonsyndata <- subset(datatww, TypeOfSite=="Non-synonymous Sites")
    nonyescpGdata<- subset(nonsyndata, combo=="noAA yesCPG")
    #nonsyn subset for no AA and no CPG
    nonNNdata <- subset(nonsyndata, combo=="noAA noCPG")
    #nonsyn subset for no AA and no CPG for a, c, g, t (for HIV all 4 should be here)
    nonsynNNa <- subset(nonNNdata, WTnt=="a")
    nonsynNNc <- subset(nonNNdata, WTnt=="c")
    nonsynNNg  <- subset(nonNNdata, WTnt=="g")
    nonsynNNt  <- subset(nonNNdata, WTnt=="t")
    #nonsyn subset for no AA and yes CPG
    nonNYdata<- subset(nonsyndata, combo=="noAA yesCPG")
    #syn subset for no AA and yes CPG for a, c, g, t (for HIV a, t )
    nonsynNYa <- subset(nonNYdata, WTnt=="a")
    nonsynNYc <- subset(nonNYdata, WTnt=="c")
    nonsynNYg  <- subset(nonNYdata, WTnt=="g")
    nonsynNYt  <- subset(nonNYdata, WTnt=="t")
    #syn subset for yes AA and no CPG 
    nonYNdata <- subset(nonsyndata, combo=="yesAA noCPG")
    #syn subset for yes AA and no CPG for a, c, g, t (for HIV all 4 should be here)
    nonsynYNa <- subset(nonYNdata, WTnt=="a")
    nonsynYNc <- subset(nonYNdata, WTnt=="c")
    nonsynYNg  <- subset(nonYNdata, WTnt=="g")
    nonsynYNt  <- subset(nonYNdata, WTnt=="t")
    #syn subset for yes AA and yes CPG
    nonYYdata <- subset(nonsyndata, combo=="yesAA yesCPG")
    #syn subset for yes AA and yes CPG for a, c, g, t (for HIV a, t)
    nonsynYYa <- subset(nonYYdata, WTnt=="a")
    nonsynYYc <- subset(nonYYdata, WTnt=="c")
    nonsynYYg  <- subset(nonYYdata, WTnt=="g")
    nonsynYYt  <- subset(nonYYdata, WTnt=="t")
    
    
    #the actual graph is here.
    graph <-ggplot(aes(factor(xvalue), MeanFreq), data = datatww)+
        #log scale to make the data eaisier to see
        scale_y_log10() +
        #scale_x_discrete makes a point at each number, breaks lets us seperate into sections for each nucecotide
        scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),breaks=c("2","6","10", "14"), labels=c("a -> g", "g -> a", "c -> t","t -> c"))+
        #syn data, jitter seperates the points
        geom_jitter(data= syndata,aes(colour = syndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
        #facet_wrap splits graph between syn and nonsyn
        facet_wrap(~ TypeOfSite)+
        #nonsyn data, jitter seperates the points
        geom_jitter(data= nonsyndata,aes(colour = nonsyndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
        #draws lines between the nucelotides 
        geom_vline(xintercept = 4.5, linetype="solid", color = "black", size=0.9) +
        geom_vline(xintercept = 8.1, linetype="solid", color = "black", size=0.9) +
        geom_vline(xintercept = 12, linetype="solid", color = "black", size=0.9)+
        #give points new colors and lables the colors
        scale_color_manual(labels = c("No drastic AA change (non-Cpg-forming)","Drastic AA change (non-Cpg-forming)","Drastic AA change (Cpg-forming)", "No drastic AA change (Cpg-forming)"), values = c("firebrick", "darkolivegreen","goldenrod3", "royalblue3")) +
        #labels X and Y axis
        labs(x="Mutation Type", y="Mutation Frquency",col=" ")
    
    #graphing other info basied on conditions, used as a check to see if data is being read corectly
    #if statments checks to see if the dataset exist. If things that should exist show up on the graph you would be able to determine what subset of data needs to be checked
    if (nrow(synNNa)!=0) {
        synNNaconf <- t.test(synNNa$MeanFreq)$conf.int
        graph<- graph + geom_errorbar(data = synNNa, aes(ymin = synNNaconf[[1]], ymax = synNNaconf[[2]], width = 0.2))
        graph<- graph + geom_point(data =synNNa, aes('1',mean(synNNa$MeanFreq)))
    }
    if (nrow(synNNc)!=0) {
        synNNcconf <- t.test(synNNc$MeanFreq)$conf.int
        graph<- graph + geom_errorbar(data = synNNc, aes(ymin = synNNcconf[[1]], ymax = synNNcconf[[2]], width = 0.2))
        graph<- graph + geom_point(data =synNNc, aes('9',mean(synNNc$MeanFreq)))
    } 
    if (nrow(synNNg)!=0) {
        synNNgconf <- t.test(synNNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNNg, aes(ymin = synNNgconf[[1]], ymax = synNNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNNg, aes('5',mean(synNNg$MeanFreq)))
    } 
    if (nrow(synNNt)!=0) {
        synNNtconf <- t.test(synNNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNNt, aes(ymin = synNNtconf[[1]], ymax = synNNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNNt, aes('13',mean(synNNt$MeanFreq)))
    } 
    #synNYdata
    if (nrow(synNYa)!=0) {
        synNYaconf <- t.test(synNYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYa, aes(ymin = synNYaconf[[1]], ymax = synNYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYa, aes('2',mean(synNYa$MeanFreq)))
    } 
    if (nrow(synNYc)!=0) {
        synNYcconf <- t.test(synNYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYc, aes(ymin = synNYcconf[[1]], ymax = synNYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYc, aes('10',mean(synNYc$MeanFreq)))
    } 
    if (nrow(synNYg)!=0) {
        synNYgconf <- t.test(synNYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYg, aes(ymin = synNYgconf[[2]], ymax = synNYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYg, aes('6',mean(synNYg$MeanFreq)))
    } 
    if (nrow(synNYt)!=0) {
        synNYtconf <- t.test(synNYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYt, aes(ymin = synNYtconf[[1]], ymax = synNYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYt, aes('14',mean(synNYt$MeanFreq)))
    } 
    if (nrow(synYNa)!=0) {
        synYNaconf <- t.test(synYNa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNa, aes(ymin = synYNaconf[[1]], ymax = synYNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNa, aes('3',mean(synYNa$MeanFreq)))
    } 
    if (nrow(synYNc)!=0) {
        synYNcconf <- t.test(synYNc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNc, aes(ymin = synYNcconf[[1]], ymax = synYNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNc, aes('11',mean(synYNc$MeanFreq)))
    } 
    if (nrow(synYNg)!=0) {
        synYNgconf <- t.test(synYNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNg, aes(ymin = synYNgconf[[1]], ymax = synYNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNg, aes('7',mean(synYNg$MeanFreq)))
    } 
    if (nrow(synYNt)!=0) {
        synYNtconf <- t.test(synYNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNt, aes(ymin = synYNtconf[[1]], ymax = synYNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNt, aes('15',mean(synYNt$MeanFreq)))
    } 
    if (nrow(synYYa)!=0) {
        synYYaconf <- t.test(synYYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYa, aes(ymin = synYYaconf[[1]], ymax = synYYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYa, aes('4',mean(synYYa$MeanFreq)))
    } 
    if (nrow(synYYc)!=0) {
        synYYcconf <- t.test(synYYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYc, aes(ymin = synYYcconf[[1]], ymax = synYYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYc, aes('12', mean(synYYc$MeanFreq)))
    } 
    if (nrow(synYYg)!=0) {
        synYYgconf <- t.test(synYYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYg, aes(ymin = synYYgconf[[1]], ymax = synYYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYg, aes('8',mean(synYYg$MeanFreq)))
    } 
    if (nrow(synYYt)!=0) {
        synYYtconf <- t.test(synYYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYt, aes(ymin = synYYtconf[[1]], ymax = synYYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYt, aes('16',mean(synYYt$MeanFreq)))
    }
    #nonsyn data
    if (nrow(nonsynNNa)!=0) {
        nonsynNNaconf <- t.test(nonsynNNa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNa, aes(ymin = nonsynNNaconf[[1]], ymax = nonsynNNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNa, aes('1',mean(nonsynNNa$MeanFreq)))
    } 
    if (nrow(nonsynNNc)!=0) {
        nonsynNNcconf <- t.test(nonsynNNc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNc, aes(ymin = nonsynNNcconf[[1]], ymax = nonsynNNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNc, aes('9',mean(nonsynNNc$MeanFreq)))
    } 
    if (nrow(nonsynNNg)!=0) {
        nonsynNNgconf <- t.test(nonsynNNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNg, aes(ymin = nonsynNNgconf[[1]], ymax = nonsynNNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNg, aes('5',mean(nonsynNNg$MeanFreq)))
    } 
    if (nrow(nonsynNNt)!=0) {
        nonsynNNtconf <- t.test(nonsynNNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNt, aes(ymin = nonsynNNtconf[[1]], ymax = nonsynNNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNt, aes('13',mean(nonsynNNt$MeanFreq)))
    } 
    if (nrow(nonsynNYa)!=0) {
        nonsynNYaconf <- t.test(nonsynNYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYa, aes(ymin = nonsynNYaconf[[1]], ymax = nonsynNYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYa, aes('2',mean(nonsynNYa$MeanFreq)))
    } 
    if (nrow(nonsynNYc)!=0) {
        nonsynNYcconf <- t.test(nonsynNYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYc, aes(ymin = nonsynNYcconf[[1]], ymax = nonsynNYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYc, aes('10',mean(nonsynNYc$MeanFreq)))
    } 
    if (nrow(nonsynNYg)!=0) {
        nonsynNYgconf <- t.test(nonsynNYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYg, aes(ymin = nonsynNYgconf[[1]], ymax = nonsynNYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYg, aes('6',mean(nonsynNYg$MeanFreq)))
    }  
    if (nrow(nonsynNYt)!=0) {
        nonsynNYtconf <- t.test(nonsynNYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYt, aes(ymin = nonsynNYtconf[[1]], ymax = nonsynNYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYt, aes('14',mean(nonsynNYt$MeanFreq)))
    } 
    if (nrow(nonsynYNa)!=0) {
        nonsynYNaconf <- t.test(nonsynYNa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNa, aes(ymin = nonsynYNaconf[[1]], ymax = nonsynYNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNa, aes('3',mean(nonsynYNa$MeanFreq)))
    } 
    if (nrow(nonsynYNc)!=0) {
        nonsynYNcconf <- t.test(nonsynYNc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNc, aes(ymin = nonsynYNcconf[[1]], ymax = nonsynYNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNc, aes('11',mean(nonsynYNc$MeanFreq)))
    } 
    if (nrow(nonsynYNg)!=0) {
        nonsynYNgconf <- t.test(nonsynYNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNg, aes(ymin = nonsynYNgconf[[1]], ymax = nonsynYNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNg, aes('7',mean(nonsynYNg$MeanFreq)))
    } 
    if (nrow(nonsynYNt)!=0) {
        nonsynYNtconf <- t.test(nonsynYNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNt, aes(ymin = nonsynYNgconf[[1]], ymax = nonsynYNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNt, aes('15',mean(nonsynYNt$MeanFreq)))
    } 
    if (nrow(nonsynYYa)!=0) {
        nonsynYYaconf <- t.test(nonsynYYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYa, aes(ymin = nonsynYYaconf[[1]], ymax = nonsynYYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYa, aes('4',mean(nonsynYYa$MeanFreq)))
    } 
    if (nrow(nonsynYYc)!=0) {
        nonsynYYcconf <- t.test(nonsynYYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYc, aes(ymin = nonsynYYcconf[[1]], ymax = nonsynYYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYc, aes('12',mean(nonsynYYc$MeanFreq)))
    } 
    if (nrow(nonsynYYg)!=0) {
        nonsynYYgconf <- t.test(nonsynYYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYg, aes(ymin = nonsynYYgconf[[1]], ymax = nonsynYYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYg, aes('8',mean(nonsynYYg$MeanFreq)))
    } 
    if (nrow(nonsynYYt)!=0) {
        nonsynYYtconf <- t.test(nonsynYYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYt, aes(ymin = nonsynYYtconf[[1]], ymax = nonsynYYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYt, aes('16',mean(nonsynYYt$MeanFreq)))
    } 
    #this get rid of some of the back ground lines
    graph <- graph +  theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              #panel.border = element_blank(),
              panel.background = element_blank()) 
    
    
    print(graph)
}


comparing_mutation_graph(final_BOCA_NS1)




