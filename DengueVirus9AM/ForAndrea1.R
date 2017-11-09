

library(seqinr)

#reading the data
read.fasta("DengueVirus1.fasta_pruned.mu.trim05") -> VIRUS

# reference sequence

ref <- VIRUS[[1]]

length(ref)->Seqlenght



VIRUS_DATA <- VIRUS_DATA[3621:4124,]

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



VIRUS_DATA <- VIRUS_DATA[3621:4124,]


View(VIRUS_DATA)





####### Investigating CPG


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








Drastic ##########
###################################################### 
####################################################
####################################################
####################################################

library(seqinr)
#library(Biostrings)
install.packages("ape")
library(ape)

# names(VIRUS_DATA)=c("WTseq", "Mutseq", "WTAA", "MutAA", "WTcat", "Mutcat", "DrasticAA")

VIRUS_DATA["WTseq"] <- 0
VIRUS_DATA["Mutseq"] <- 0
VIRUS_DATA["WTAA"] <- 0
VIRUS_DATA["MutAA"] <- 0
VIRUS_DATA["WTcat"] <- 0
VIRUS_DATA["Mutcat"] <- 0
VIRUS_DATA["DrasticAA"] <- 0


curSeq <- translate(paste(VIRUS_DATA[,2], sep=" "), NAstring="X", ambiguous = FALSE, sens="F")

count <- 0

curSeq





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



#reading files and setting data  

seqs = read.dna("DengueVirus1.fasta", format = "fasta", as.character=TRUE)
a=nrow(seqs)
b=ncol(seqs)


# Categorizing the WT amino acids

for(j in 1:b){
    
    VIRUS_DATA[j,9]=amCat(VIRUS_DATA[j,7])
    
}

#creates Mut strain of current WTseq enters it into 
for (i in 1:b){
    curNuc <-  VIRUS_DATA[i,2]
    
    if(curNuc == "a" || curNuc == "g"){
        
        if(curNuc == "a"){
            VIRUS_DATA [i, 6] = "g"
        }else{
            VIRUS_DATA [i, 6] = "a"
        }
        
    }else{
        
        if(curNuc == "t"){
            VIRUS_DATA [i, 6] = "c"
        }else{
            VIRUS_DATA [i, 6] = "t"
        }
        
    }
    
}


########## WORKED UNTIL HERE ####


# Let's Creates 

VIRUS_DATA$WTseq = VIRUS_DATA$wtnt

VIRUS_DATA$MutAA <-c(0)

count <- 1
i = 1
j = 1
for(i in 1:(nrow(VIRUS_DATA)/3)){ # Workspace for MutAA filling don't run this loop
    
    for (j in 1:3){
        
        if(j == 1){                     #first codon
            
            VIRUS_DATA[count, ]$MutAA <- translate(firstC <- c(as.character(VIRUS_DATA[count,]$Mutseq), as.character(VIRUS_DATA[count + 1,]$WTseq), as.character(VIRUS_DATA[count +2,]$WTseq)))
        } 
        if(j == 2){                #second codon
            VIRUS_DATA[count, ]$MutAA <- translate(secondC <- c(as.character(VIRUS_DATA[count - 1,]$WTseq), as.character(VIRUS_DATA[count,]$Mutseq), as.character(VIRUS_DATA[count +1,]$WTseq)))
        }
        if(j == 3){                  #third codon
            VIRUS_DATA[count, ]$MutAA <- translate(thirdC <- c(as.character(VIRUS_DATA[count - 2,]$WTseq), as.character(VIRUS_DATA[count - 1,]$WTseq), as.character(VIRUS_DATA[count,]$Mutseq)))
        }
        
        count <- count + 1
        print(j)
        print(i)
        print(count)
    }
}                                           


VIRUS_DATA[,8]=seqinr::translate(paste(VIRUS_DATA[,6], sep=" "),
                                 NAstring="X", ambiguous=FALSE, sens="F")            # Translating from nucleic acid to AA

for(j in 1:b){                                              
    VIRUS_DATA[j,10]=amCat(VIRUS_DATA[j,8])                                    # Categorizing the AA 
}

# Function to compare DrasticAAChanges
DrasticChange <- function(VIRUS_DATA){
    for(h in 1:nrow(VIRUS_DATA)){
        if (VIRUS_DATA[h,]$WTcat==VIRUS_DATA[h,]$Mutcat){
            VIRUS_DATA[h,]$DrasticAA= 0                      # If WT AA category = Mut AA Category, no drastic change
        }
        if (VIRUS_DATA[h,]$WTcat!=VIRUS_DATA[h,]$Mutcat){
            VIRUS_DATA[h,]$DrasticAA= 1                      # If WT AA category =/= Mut AA Category, yes drastic change
        }
    }
    return(VIRUS_DATA)
}
VIRUS_DATA <- DrasticChange(VIRUS_DATA)
save(VIRUS_DATA,file="VIRUS_DATA.Rda")
load("VIRUS_DATA.Rda")
View(VIRUS_DATA)
##############
######################
#########################
#######################
for (h in 1:nrow(VIRUS_DATA)){ #looks at each row in the dataframe df
    if(VIRUS_DATA[h,"MutAA"]== VIRUS_DATA[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
        VIRUS_DATA[h,"TypeOfSite"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(VIRUS_DATA[h,"MutAA"] != VIRUS_DATA[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
        if(VIRUS_DATA[h,"MutAA"]=="*"){ #and it the MUTAA value equal "*"
            VIRUS_DATA[h,"TypeOfSite"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
        else {
            VIRUS_DATA[h,"TypeOfSite"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
        }
    }
}
View(VIRUS_DATA)

##################################
###################################
##############################
#input should be dataframe
#the command should be something like the function pulls out information
#from the dataframe and uses that info to plot
#output should be plot of freq vs. Position (num), colored by TypeofSite
plot_NonSynNon = function(VIRUS_DATA) {
    plot(VIRUS_DATA$num, VIRUS_DATA$freq+0.1, type = "p", main = "Nucleotide Positon (num) vs. Mean Frequency (freq) of Virus", 
         log = "y", xlab = "num", ylab = "log of freq", col = as.factor(VIRUS_DATA$TypeOfSite))
}


#explanation of some plot arguments
#df$freq+0.1 we had to do +0.1 because our graph looked better that way, yours may or may not need more zeros
#type = "p" makes a scatterplot
#col = as.factor(df$TypeOfSite) was so that we could color everything based on TypeOfSite, but first we had to 
#ensure that R read it as factors; otherwise it wouldn't color it b/c "nonsyn" was not a color (tl;dr R is picky)


plot_NonSynNon(VIRUS_DATA) #this is what using the function should look like, our df was titled final_BOCA


#making a legend
legend("topleft",
       c("syn", "nonsyn", "nonsense"),
       lty = c(1,1), lwd=4, cex=0.70,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")

#who is in the team: Shannel, Milo, Ricky, Natalie
#who helped make the function: mostly Shannel, sprinkles of Milo and Ricky here and there

View(VIRUS_DATA)


#Night_Crew Team
#Members Victoria, Sam, Gordon, Mordy 
#setwd("~/Desktop/Git")

#install.packages('ggplot2')
install.packages('scales')
#scales is used to log the y axis in the graoh and ggplot2 is used to create the graph
library(scales)
library(ggplot2)

#this function graphs transition mutations separated by Synonymous Sites and Nonsynonymous sites
# this graph show the frequencies of the changes and also seperates combinations of drastic AA and CpG sites by color
#Start of the graph function, here VIRUS_DATA is your VIRUS_DATAset. Make sure you have these columns and the column names are the same.DrasticAA, CPG/CPG, TypeOfSite, num, freqor/freq, wtnt/wtnt 
#In column TypeOfSite make sure that syn, nonsyn are present. No - or full names
comparing_mutation_graph = function(VIRUS_DATA){
    #building the combo lines to help sort (used later when info is subsetted) 
    #This will change your colomn name to what we orginally used to make our graph
    colnames(VIRUS_DATA)[which(names(VIRUS_DATA) == "wtnt")] <- "wtnt"
    colnames(VIRUS_DATA)[which(names(VIRUS_DATA) == "freq")] <- "freq"
    colnames(VIRUS_DATA)[which(names(VIRUS_DATA) == "CPG")] <- "CPG"
    VIRUS_DATA$combo<- (VIRUS_DATA$DrasticAA*3) + (VIRUS_DATA$CPG*2)
    VIRUS_DATA$combo<- as.factor(VIRUS_DATA$combo)
    levels(VIRUS_DATA$combo) <- gsub("0", "noAA noCPG", levels(VIRUS_DATA$combo))
    levels(VIRUS_DATA$combo) <- gsub("2", "noAA yesCPG", levels(VIRUS_DATA$combo))
    levels(VIRUS_DATA$combo) <- gsub("3", "yesAA noCPG", levels(VIRUS_DATA$combo))
    levels(VIRUS_DATA$combo) <- gsub("5", "yesAA yesCPG", levels(VIRUS_DATA$combo))
    #This will change your colomn name to what we orginally used to make our graph
    colnames(VIRUS_DATA)[which(names(VIRUS_DATA) == "wtnt")] <- "wtnt"
    colnames(VIRUS_DATA)[which(names(VIRUS_DATA) == "freq")] <- "freq"
    #Subsetting for the VIRUS_DATA that we really want just Synonymous and nonynonymous
    VIRUS_DATAtww <- subset(VIRUS_DATA, TypeOfSite=="syn" | TypeOfSite=="nonsyn")
    #if your dataframe has "-" try using the following line of code
    #datatww<-datatww[!(which(datatww$wtnt=="-")),]
    VIRUS_DATAtww$TypeOfSite <- ifelse(VIRUS_DATAtww$TypeOfSite == "syn", "Synonymous Sites", "Non-synonymous Sites")
    
    #gives values to each nucecotide and if they have a cetain combo
    #first need to introduce new columes for collecting the xvalue and the color
    VIRUS_DATAtww$xvalue<- 0 #giving each nucecotide a value depending on whether or not there is an AA change a CpG site and the nucecotide
    VIRUS_DATAtww$color <- 0 # four colors for the four combos: Blue = no DrasticAA/no CPG, Red = no DrasticAA/yes CPG, Green = yes DrasticAA/no CPG, Purple = yes DrasticAA/yes CPG 
    #for loop adds color and a xvalue to the dataset, for all nucecotides
    for (i in 1:length(VIRUS_DATAtww$num)) {
        if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "a") {
            VIRUS_DATAtww$xvalue[i] <- 1
            VIRUS_DATAtww$color[i] <- "blue"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "a") {
            VIRUS_DATAtww$xvalue[i] <- 2
            VIRUS_DATAtww$color[i] <- "red"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "a") {
            VIRUS_DATAtww$xvalue[i] <- 3
            VIRUS_DATAtww$color[i] <- "green"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "a") {
            VIRUS_DATAtww$xvalue[i] <- 4
            VIRUS_DATAtww$color[i] <- "purple"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "g") {
            VIRUS_DATAtww$xvalue[i] <- 5
            VIRUS_DATAtww$color[i] <- "blue"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "g") {
            VIRUS_DATAtww$xvalue[i] <- 6
            VIRUS_DATAtww$color[i] <- "red"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "g") {
            VIRUS_DATAtww$xvalue[i] <- 7
            VIRUS_DATAtww$color[i] <- "green"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "g") {
            VIRUS_DATAtww$xvalue[i] <- 8
            VIRUS_DATAtww$color[i] <- "purple"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "c") {
            VIRUS_DATAtww$xvalue[i] <- 9
            VIRUS_DATAtww$color[i] <- "blue"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "c") {
            VIRUS_DATAtww$xvalue[i] <- 10 
            VIRUS_DATAtww$color[i] <- "red"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "c") {
            VIRUS_DATAtww$xvalue[i] <- 11
            VIRUS_DATAtww$color[i] <- "green"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "c") {
            VIRUS_DATAtww$xvalue[i] <- 12
            VIRUS_DATAtww$color[i] <- "purple"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "t") {
            VIRUS_DATAtww$xvalue[i] <- 13
            VIRUS_DATAtww$color[i] <- "blue"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 0 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "t") {
            VIRUS_DATAtww$xvalue[i] <- 14
            VIRUS_DATAtww$color[i] <- "red"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 0 && VIRUS_DATAtww$wtnt[i] == "t") {
            VIRUS_DATAtww$xvalue[i] <- 15
            VIRUS_DATAtww$color[i] <- "green"
        } else if (VIRUS_DATAtww$DrasticAA[i] == 1 && VIRUS_DATAtww$CPG[i] == 1 && VIRUS_DATAtww$wtnt[i] == "t") {
            VIRUS_DATAtww$xvalue[i] <- 16
            VIRUS_DATAtww$color[i] <- "purple"
        }
    }
    #subsetiing all the VIRUS_DATA by amino acid, drastic AA, and CpG forming: Helps with analyst 
    #syn subset 
    synVIRUS_DATA <- subset(VIRUS_DATAtww, TypeOfSite=="Synonymous Sites")
    #syn subset for no AA and no CPG
    synNNVIRUS_DATA <- subset(synVIRUS_DATA, combo=="noAA noCPG")
    #syn subset for no AA and no CPG for a, c, g, t (for HIV all 4 should be here)
    synNNa <- subset(synNNVIRUS_DATA, wtnt=="a")
    synNNc <- subset(synNNVIRUS_DATA, wtnt=="c")
    synNNg  <- subset(synNNVIRUS_DATA, wtnt=="g")
    synNNt  <- subset(synNNVIRUS_DATA, wtnt=="t")
    #syn subset for no AA and yes CPG
    synNYVIRUS_DATA<- subset(synVIRUS_DATA, combo=="noAA yesCPG")
    #syn subset for no AA and yes CPG for a, c, g, t (for HIV a, t )
    synNYa <- subset(synNYVIRUS_DATA, wtnt=="a")
    synNYc <- subset(synNYVIRUS_DATA, wtnt=="c")
    synNYg  <- subset(synNYVIRUS_DATA, wtnt=="g")
    synNYt  <- subset(synNYVIRUS_DATA, wtnt=="t")
    #syn subset for yes AA and no CPG 
    synYNVIRUS_DATA <- subset(synVIRUS_DATA, combo=="yesAA noCPG")
    #syn subset for yes AA and no CPG for a, c, g, t (none for HIV)
    synYNa <- subset(synYNVIRUS_DATA, wtnt=="a")
    synYNc <- subset(synYNVIRUS_DATA, wtnt=="c")
    synYNg  <- subset(synYNVIRUS_DATA, wtnt=="g")
    synYNt  <- subset(synYNVIRUS_DATA, wtnt=="t")
    #syn subset for yes AA and yes CPG
    synYYVIRUS_DATA <- subset(synVIRUS_DATA, combo=="yesAA yesCPG")
    #syn subset for yes AA and yes CPG for a, c, g, t (none for HIV)
    synYYa <- subset(synYYVIRUS_DATA, wtnt=="a")
    synYYc <- subset(synYYVIRUS_DATA, wtnt=="c")
    synYYg  <- subset(synYYVIRUS_DATA, wtnt=="g")
    synYYt  <- subset(synYYVIRUS_DATA, wtnt=="t")
    
    #nonsyn sub set 
    nonsynVIRUS_DATA <- subset(VIRUS_DATAtww, TypeOfSite=="Non-synonymous Sites")
    nonyescpGVIRUS_DATA<- subset(nonsynVIRUS_DATA, combo=="noAA yesCPG")
    #nonsyn subset for no AA and no CPG
    nonNNVIRUS_DATA <- subset(nonsynVIRUS_DATA, combo=="noAA noCPG")
    #nonsyn subset for no AA and no CPG for a, c, g, t (for HIV all 4 should be here)
    nonsynNNa <- subset(nonNNVIRUS_DATA, wtnt=="a")
    nonsynNNc <- subset(nonNNVIRUS_DATA, wtnt=="c")
    nonsynNNg  <- subset(nonNNVIRUS_DATA, wtnt=="g")
    nonsynNNt  <- subset(nonNNVIRUS_DATA, wtnt=="t")
    #nonsyn subset for no AA and yes CPG
    nonNYVIRUS_DATA<- subset(nonsynVIRUS_DATA, combo=="noAA yesCPG")
    #syn subset for no AA and yes CPG for a, c, g, t (for HIV a, t )
    nonsynNYa <- subset(nonNYVIRUS_DATA, wtnt=="a")
    nonsynNYc <- subset(nonNYVIRUS_DATA, wtnt=="c")
    nonsynNYg  <- subset(nonNYVIRUS_DATA, wtnt=="g")
    nonsynNYt  <- subset(nonNYVIRUS_DATA, wtnt=="t")
    #syn subset for yes AA and no CPG 
    nonYNVIRUS_DATA <- subset(nonsynVIRUS_DATA, combo=="yesAA noCPG")
    #syn subset for yes AA and no CPG for a, c, g, t (for HIV all 4 should be here)
    nonsynYNa <- subset(nonYNVIRUS_DATA, wtnt=="a")
    nonsynYNc <- subset(nonYNVIRUS_DATA, wtnt=="c")
    nonsynYNg  <- subset(nonYNVIRUS_DATA, wtnt=="g")
    nonsynYNt  <- subset(nonYNVIRUS_DATA, wtnt=="t")
    #syn subset for yes AA and yes CPG
    nonYYVIRUS_DATA <- subset(nonsynVIRUS_DATA, combo=="yesAA yesCPG")
    #syn subset for yes AA and yes CPG for a, c, g, t (for HIV a, t)
    nonsynYYa <- subset(nonYYVIRUS_DATA, wtnt=="a")
    nonsynYYc <- subset(nonYYVIRUS_DATA, wtnt=="c")
    nonsynYYg  <- subset(nonYYVIRUS_DATA, wtnt=="g")
    nonsynYYt  <- subset(nonYYVIRUS_DATA, wtnt=="t")
    
    
    #the actual graph is here.
    graph <-ggplot(aes(factor(xvalue), freq), data = VIRUS_DATAtww)+
        #log scale to make the data eaisier to see
        scale_y_log10() +
        #scale_x_discrete makes a point at each number, breaks lets us seperate into sections for each nucecotide
        scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),breaks=c("2","6","10", "14"), labels=c("a -> g", "g -> a", "c -> t","t -> c"))+
        #syn data, jitter seperates the points
        geom_jitter(data= synVIRUS_DATA,aes(colour = synVIRUS_DATA$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
        #facet_wrap splits graph between syn and nonsyn
        facet_wrap(~ TypeOfSite)+
        #nonsyn data, jitter seperates the points
        geom_jitter(data= nonsynVIRUS_DATA,aes(colour = nonsynVIRUS_DATA$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
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
        synNNaconf <- t.test(synNNa$freq)$conf.int
        graph<- graph + geom_errorbar(data = synNNa, aes(ymin = synNNaconf[[1]], ymax = synNNaconf[[2]], width = 0.2))
        graph<- graph + geom_point(data =synNNa, aes('1',mean(synNNa$freq)))
    }
    if (nrow(synNNc)!=0) {
        synNNcconf <- t.test(synNNc$freq)$conf.int
        graph<- graph + geom_errorbar(data = synNNc, aes(ymin = synNNcconf[[1]], ymax = synNNcconf[[2]], width = 0.2))
        graph<- graph + geom_point(data =synNNc, aes('9',mean(synNNc$freq)))
    } 
    if (nrow(synNNg)!=0) {
        synNNgconf <- t.test(synNNg$freq)$conf.int
        graph<- graph +geom_errorbar(data = synNNg, aes(ymin = synNNgconf[[1]], ymax = synNNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNNg, aes('5',mean(synNNg$freq)))
    } 
    if (nrow(synNNt)!=0) {
        synNNtconf <- t.test(synNNt$freq)$conf.int
        graph<- graph +geom_errorbar(data = synNNt, aes(ymin = synNNtconf[[1]], ymax = synNNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNNt, aes('13',mean(synNNt$freq)))
    } 
    #synNYdata
    if (nrow(synNYa)!=0) {
        synNYaconf <- t.test(synNYa$freq)$conf.int
        graph<- graph +geom_errorbar(data = synNYa, aes(ymin = synNYaconf[[1]], ymax = synNYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYa, aes('2',mean(synNYa$freq)))
    } 
    if (nrow(synNYc)!=0) {
        synNYcconf <- t.test(synNYc$freq)$conf.int
        graph<- graph +geom_errorbar(data = synNYc, aes(ymin = synNYcconf[[1]], ymax = synNYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYc, aes('10',mean(synNYc$freq)))
    } 
    if (nrow(synNYg)!=0) {
        synNYgconf <- t.test(synNYg$freq)$conf.int
        graph<- graph +geom_errorbar(data = synNYg, aes(ymin = synNYgconf[[2]], ymax = synNYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYg, aes('6',mean(synNYg$freq)))
    } 
    if (nrow(synNYt)!=0) {
        synNYtconf <- t.test(synNYt$freq)$conf.int
        graph<- graph +geom_errorbar(data = synNYt, aes(ymin = synNYtconf[[1]], ymax = synNYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYt, aes('14',mean(synNYt$freq)))
    } 
    if (nrow(synYNa)!=0) {
        synYNaconf <- t.test(synYNa$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYNa, aes(ymin = synYNaconf[[1]], ymax = synYNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNa, aes('3',mean(synYNa$freq)))
    } 
    if (nrow(synYNc)!=0) {
        synYNcconf <- t.test(synYNc$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYNc, aes(ymin = synYNcconf[[1]], ymax = synYNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNc, aes('11',mean(synYNc$freq)))
    } 
    if (nrow(synYNg)!=0) {
        synYNgconf <- t.test(synYNg$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYNg, aes(ymin = synYNgconf[[1]], ymax = synYNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNg, aes('7',mean(synYNg$freq)))
    } 
    if (nrow(synYNt)!=0) {
        synYNtconf <- t.test(synYNt$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYNt, aes(ymin = synYNtconf[[1]], ymax = synYNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNt, aes('15',mean(synYNt$freq)))
    } 
    if (nrow(synYYa)!=0) {
        synYYaconf <- t.test(synYYa$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYYa, aes(ymin = synYYaconf[[1]], ymax = synYYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYa, aes('4',mean(synYYa$freq)))
    } 
    if (nrow(synYYc)!=0) {
        synYYcconf <- t.test(synYYc$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYYc, aes(ymin = synYYcconf[[1]], ymax = synYYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYc, aes('12', mean(synYYc$freq)))
    } 
    if (nrow(synYYg)!=0) {
        synYYgconf <- t.test(synYYg$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYYg, aes(ymin = synYYgconf[[1]], ymax = synYYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYg, aes('8',mean(synYYg$freq)))
    } 
    if (nrow(synYYt)!=0) {
        synYYtconf <- t.test(synYYt$freq)$conf.int
        graph<- graph +geom_errorbar(data = synYYt, aes(ymin = synYYtconf[[1]], ymax = synYYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYt, aes('16',mean(synYYt$freq)))
    }
    #nonsyn data
    if (nrow(nonsynNNa)!=0) {
        nonsynNNaconf <- t.test(nonsynNNa$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNa, aes(ymin = nonsynNNaconf[[1]], ymax = nonsynNNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNa, aes('1',mean(nonsynNNa$freq)))
    } 
    if (nrow(nonsynNNc)!=0) {
        nonsynNNcconf <- t.test(nonsynNNc$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNc, aes(ymin = nonsynNNcconf[[1]], ymax = nonsynNNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNc, aes('9',mean(nonsynNNc$freq)))
    } 
    if (nrow(nonsynNNg)!=0) {
        nonsynNNgconf <- t.test(nonsynNNg$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNg, aes(ymin = nonsynNNgconf[[1]], ymax = nonsynNNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNg, aes('5',mean(nonsynNNg$freq)))
    } 
    if (nrow(nonsynNNt)!=0) {
        nonsynNNtconf <- t.test(nonsynNNt$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNt, aes(ymin = nonsynNNtconf[[1]], ymax = nonsynNNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNt, aes('13',mean(nonsynNNt$freq)))
    } 
    if (nrow(nonsynNYa)!=0) {
        nonsynNYaconf <- t.test(nonsynNYa$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYa, aes(ymin = nonsynNYaconf[[1]], ymax = nonsynNYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYa, aes('2',mean(nonsynNYa$freq)))
    } 
    if (nrow(nonsynNYc)!=0) {
        nonsynNYcconf <- t.test(nonsynNYc$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYc, aes(ymin = nonsynNYcconf[[1]], ymax = nonsynNYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYc, aes('10',mean(nonsynNYc$freq)))
    } 
    if (nrow(nonsynNYg)!=0) {
        nonsynNYgconf <- t.test(nonsynNYg$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYg, aes(ymin = nonsynNYgconf[[1]], ymax = nonsynNYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYg, aes('6',mean(nonsynNYg$freq)))
    }  
    if (nrow(nonsynNYt)!=0) {
        nonsynNYtconf <- t.test(nonsynNYt$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYt, aes(ymin = nonsynNYtconf[[1]], ymax = nonsynNYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYt, aes('14',mean(nonsynNYt$freq)))
    } 
    if (nrow(nonsynYNa)!=0) {
        nonsynYNaconf <- t.test(nonsynYNa$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNa, aes(ymin = nonsynYNaconf[[1]], ymax = nonsynYNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNa, aes('3',mean(nonsynYNa$freq)))
    } 
    if (nrow(nonsynYNc)!=0) {
        nonsynYNcconf <- t.test(nonsynYNc$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNc, aes(ymin = nonsynYNcconf[[1]], ymax = nonsynYNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNc, aes('11',mean(nonsynYNc$freq)))
    } 
    if (nrow(nonsynYNg)!=0) {
        nonsynYNgconf <- t.test(nonsynYNg$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNg, aes(ymin = nonsynYNgconf[[1]], ymax = nonsynYNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNg, aes('7',mean(nonsynYNg$freq)))
    } 
    if (nrow(nonsynYNt)!=0) {
        nonsynYNtconf <- t.test(nonsynYNt$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNt, aes(ymin = nonsynYNgconf[[1]], ymax = nonsynYNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNt, aes('15',mean(nonsynYNt$freq)))
    } 
    if (nrow(nonsynYYa)!=0) {
        nonsynYYaconf <- t.test(nonsynYYa$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYa, aes(ymin = nonsynYYaconf[[1]], ymax = nonsynYYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYa, aes('4',mean(nonsynYYa$freq)))
    } 
    if (nrow(nonsynYYc)!=0) {
        nonsynYYcconf <- t.test(nonsynYYc$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYc, aes(ymin = nonsynYYcconf[[1]], ymax = nonsynYYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYc, aes('12',mean(nonsynYYc$freq)))
    } 
    if (nrow(nonsynYYg)!=0) {
        nonsynYYgconf <- t.test(nonsynYYg$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYg, aes(ymin = nonsynYYgconf[[1]], ymax = nonsynYYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYg, aes('8',mean(nonsynYYg$freq)))
    } 
    if (nrow(nonsynYYt)!=0) {
        nonsynYYtconf <- t.test(nonsynYYt$freq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYt, aes(ymin = nonsynYYtconf[[1]], ymax = nonsynYYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYt, aes('16',mean(nonsynYYt$freq)))
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

comparing_mutation_graph(VIRUS_DATA)
#data is your dataframe



