

library(seqinr)

#reading the data
read.fasta("DengueVirus1.fasta") -> VIRUS

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
#output should be plot of Meanfreq vs. Position (num), colored by TypeofSite
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
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")

#who is in the team: Shannel, Milo, Ricky, Natalie
#who helped make the function: mostly Shannel, sprinkles of Milo and Ricky here and there

View(VIRUS_DATA)
