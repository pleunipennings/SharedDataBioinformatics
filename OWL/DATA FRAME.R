





library(seqinr)

read.fasta("HumanParaInfluenza/humanparainfluenzavirus1_F.fasta_pruned.mu.trim05")->VIRUS

# reference sequence

ref <- VIRUS[[1]]
print(ref)

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



which(VIRUS_DATA$CPG=="0")->noCpG
print(noCpG)
#Identify CpG sites
which(VIRUS_DATA$CPG=="1")->CpG
print(CpG)
#Identify mean frequency in association with noCpG sites
VIRUS_DATA[noCpG,"freq"]->Freq
print(Freq)
#Identify mean frequency in association with CpG sites
VIRUS_DATA[CpG,"freq"]->Freq2
print(Freq2)


#Making a plot of mean frequency vs. CpG/nonCpG location
plot(noCpG, Freq+0.0001,ylim=c(0.0001,0.5), col="black",pch=21.25, bg=rgb(1,0,0,0.5), main="CpG/noCpG location vs. Mean Frequency", xlab = "Location", ylab = "Mean Frequency")

#Overlapping the two graphs
points(CpG, Freq2+0.0001, col="black",pch=21.25, bg=rgb(0,0,1,0.5))


#legend with circle instead of line  
legend("topleft",c("noCpG","CpG"),cex=0.6, 
       col ="black",bg=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),bty="n",
       title="Legend",inset=.02,pch =21,pt.bg = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))


#  need mutaa and wt Note that this depends on the columns MUTAA and WTAA being already complete!
VIRUS_DATA$MUTAA<-("")
VIRUS_DATA$WTAA<-("")
VIRUS_DATA$TypeOfSite<-("")
transition<-function(basepair){
    #basepair<-("A", "C", "T", "G"),
    if(basepair=="a") {return("g")}
    if(basepair=="g") {return ("a")}
    if(basepair=="t") {return ("c")}
    if(basepair=="c") {return ("t")}
}

VIRUS_DATA$wtnt<-as.character(VIRUS_DATA$wtnt)

for (i in seq.int(1,nrow(VIRUS_DATA),3)) {
    codon <- c(VIRUS_DATA$wtnt[i],VIRUS_DATA$wtnt[i+1],VIRUS_DATA$wtnt[i+2])
    VIRUS_DATA$WTAA[i] = seqinr::translate(codon)
    VIRUS_DATA$WTAA[i+1] = seqinr::translate(codon)
    VIRUS_DATA$WTAA[i+2] = seqinr::translate(codon)
    copy = codon
    copy[1] = transition(codon[1])
    VIRUS_DATA$MUTAA[i] = seqinr::translate(copy)
    copy = codon
    copy[2] = transition(codon[2])
    VIRUS_DATA$MUTAA[i+1] = seqinr::translate(copy)
    copy = codon
    copy[3] = transition(codon[3])
    VIRUS_DATA$MUTAA[i+2] = seqinr::translate(copy)
}

for (h in 1:15466){
    
    if(VIRUS_DATA[h,"MUTAA"]== VIRUS_DATA[h,"WTAA"]){
        
        VIRUS_DATA[h,"TypeOfSite"] = "syn"
        
    }
    
    if(VIRUS_DATA[h,"MUTAA"] != VIRUS_DATA[h,"WTAA"]){
        
        if(VIRUS_DATA[h,"MUTAA"]=="*"){
            
            VIRUS_DATA[h,"TypeOfSite"] = "non-sense"}
        
        else {
            
            VIRUS_DATA[h,"TypeOfSite"] = "non-syn"
            
        }
        
    }
    
}

