




library(seqinr)

DEN<-read.fasta("DengueViruses/DengueViruses/DengueVirus1.fasta_pruned.mu.trim05")

read.fasta("DengueVirus1.fasta") -> DEN

# reference sequence

ref <- DEN[[1]]

# dataframe columns

num <- c(1:10689)

wtnt <- c()

freq <- c()

# for freq calculation later

absfreq <- c(rep(0,10689))

totalcount <- c(rep(0,10689))

# average WT calculation

# counts number of each nucleotide in each position

acount <- c(rep(0,10689))

gcount <- c(rep(0,10689))

ccount <- c(rep(0,10689))

tcount <- c(rep(0,10689))

nuc <- c()

# same as line 20 comment

for (i in 1:length(DEN)) {
    
    sequence <- DEN[[i]]
    
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

for (i in 1:length(DEN)) {
    
    sequence <- DEN[[i]]
    
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

Dengue <- data.frame(num,wtnt,freq)

View(Dengue)


# Investigating CPG

read.fasta("DengueVirus1.fasta") -> Dengue


# Collapses data into a string and sets to variable STRING
paste(Dengue$wtnt, collapse = '') -> STRING
STRING

#Looks for pattern tg in the data STRING, sets location sites to variable TG
gregexpr(pattern ='tg',STRING ) -> TG
TG

# Insert LIST TG into a data frame BELL
BELL <- data.frame(matrix(unlist(TG)))
BELL

# create new column name CPG with values zero
Dengue$CPG<- 0

# Inserting value 1 in every TG site using BELL
Dengue[BELL[,1],"CPG"] <- 1

View(Dengue)

# For CA sites: 

#Looks for pattern ca in the data STRING, sets location sites to variable CA
gregexpr(pattern ='ca',STRING ) -> CA
CA

# Insert LIST CA into a data frame CASITES
CASITES <- data.frame(matrix(unlist(CA)))
CASITES

# Inserting value 1 in every CA site using CA
Dengue[CASITES[,1]+1,"CPG"] <- 1



